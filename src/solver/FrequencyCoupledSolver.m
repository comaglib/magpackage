classdef FrequencyCoupledSolver < handle
    % FREQUENCYCOUPLEDSOLVER 场路耦合非线性时谐求解器 (v1.0)
    %
    % 功能:
    %   1. 求解 A-V-I 频域耦合系统 (Vector Potential + Modified V + Circuit Current).
    %   2. 支持线圈与实心导体(涡流)共存。
    %   3. 自动线性/非线性检测与自适应松弛迭代。
    %
    % 输出:
    %   V_sol 为物理电压 (已自动转换 V_phys = 1j*w*V_mod)。
    
    properties
        Assembler
        LinearSolver
        Frequency
        
        % --- 迭代配置 ---
        MaxIterations = 100
        Tolerance = 1e-5
        
        % 自适应松弛 (Line Search 变体 for Fixed Point)
        Relaxation      = 0.5
        MinRelaxation   = 0.05
        MaxRelaxation   = 1.0
        RelaxationGrow  = 1.05
        RelaxationShrink= 0.6
        
        % --- 正则化 ---
        AutoRegularization = true
        RegScaleRel = 1e-6
    end
    
    methods
        function obj = FrequencyCoupledSolver(assembler, freq)
            obj.Assembler = assembler;
            obj.Frequency = freq;
            obj.LinearSolver = LinearSolver('Auto');
            
            % 场路耦合矩阵通常是非对称的 (Unsymmetric)
            % FEM块对称, 但 A-I 耦合块 (-W vs -jwW') 不对称
            obj.LinearSolver.MumpsSymmetry = 0; 
            obj.LinearSolver.MumpsICNTL.i14 = 100; % 增强稳定性
        end
        
        function [A_sol, V_sol, I_sol, info] = solve(obj, space_A, space_V, matLib, sigmaMap, ...
                                                     circuitProps, windingObj, ...
                                                     fixedDofs_A, fixedDofs_V, init_A)
            if nargin < 10, init_A = []; end
            if nargin < 9, fixedDofs_V = []; end
            
            omega = 2 * pi * obj.Frequency;
            if omega == 0, omega = 1e-9; end
            
            fprintf('======================================================\n');
            fprintf('   Frequency Coupled Solver (%.1f Hz) - A-V-I Form    \n', obj.Frequency);
            fprintf('======================================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            
            % --- 1. 自由度分配 ---
            [~, offset_A] = dofHandler.distributeDofs(space_A);
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            
            % V 场 (仅导体)
            conductingTags = obj.getConductingRegions(sigmaMap);
            enable_V = ~isempty(conductingTags) && ~isempty(space_V);
            offset_V = -1;
            auto_ground_idx = [];
            
            if enable_V
                if ~dofHandler.DofMaps.isKey(space_V.toString())
                    [~, offset_V] = dofHandler.distributeDofs(space_V, conductingTags);
                else
                    offset_V = dofHandler.SpaceOffsets(space_V.toString());
                end
                
                % [Auto-Grounding] 防止 V 场悬浮
                hasFixedV = ~isempty(fixedDofs_V) && any(fixedDofs_V);
                if ~hasFixedV
                    fprintf('   [Auto-Fix] V-field floating. Grounding 1st node.\n');
                    auto_ground_idx = offset_V + 1;
                end
            else
                fprintf('   [Info] No eddy current regions. Pure Coil mode.\n');
            end
            
            numFemDofs = dofHandler.NumGlobalDofs;
            offset_I = numFemDofs; 
            numTotalDofs = numFemDofs + 1; % 增加电流自由度
            
            % --- 2. 预计算常数矩阵 ---
            fprintf('   [Init] Assembling invariant matrices...\n');
            
            % 扩展矩阵尺寸到 TotalDofs
            M_sigma = sparse(numTotalDofs, numTotalDofs);
            C_AV    = sparse(numTotalDofs, numTotalDofs);
            L_VV    = sparse(numTotalDofs, numTotalDofs);
            M_reg   = sparse(numTotalDofs, numTotalDofs);
            
            % 组装 FEM 部分
            if enable_V
                M_sigma_loc = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
                C_AV_loc    = obj.Assembler.assembleCoupling(space_A, space_V, sigmaMap);
                L_VV_loc    = obj.Assembler.assembleScalarLaplacian(space_V, sigmaMap);
                
                [i,j,v] = find(M_sigma_loc); M_sigma = sparse(i,j,v, numTotalDofs, numTotalDofs);
                [i,j,v] = find(C_AV_loc);    C_AV    = sparse(i,j,v, numTotalDofs, numTotalDofs);
                [i,j,v] = find(L_VV_loc);    L_VV    = sparse(i,j,v, numTotalDofs, numTotalDofs);
            end
            
            % 正则化矩阵 (仅针对 A 场)
            if obj.AutoRegularization
                M_reg_loc = obj.Assembler.assembleMass(space_A);
                [i,j,v] = find(M_reg_loc);
                M_reg = sparse(i,j,v, numTotalDofs, numTotalDofs);
            end
            
            % 绕组耦合向量 W (扩展)
            W_loc = obj.Assembler.assembleWinding(space_A, windingObj);
            W_vec = sparse(numTotalDofs, 1);
            W_vec(1:length(W_loc)) = W_loc;
            
            % --- 3. 边界条件合并 ---
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            all_fixed = idx_A;
            if enable_V
                if isempty(fixedDofs_V), idx_V = [];
                elseif islogical(fixedDofs_V), idx_V = find(fixedDofs_V);
                else, idx_V = fixedDofs_V(:); 
                end
                if ~isempty(auto_ground_idx), idx_V = unique([idx_V; auto_ground_idx]); end
                idx_V = idx_V(idx_V <= numFemDofs);
                all_fixed = [all_fixed; idx_V];
            end
            % 电流自由度 I 通常不固定
            
            % --- 4. 线性/非线性检测 ---
            isLinear = obj.checkLinearity(matLib);
            if isLinear
                max_iter = 1;
                fprintf('   [Info] Linear material detected. One-shot solve.\n');
            else
                max_iter = obj.MaxIterations;
            end
            
            % --- 5. 迭代求解 ---
            Sol_global = zeros(numTotalDofs, 1);
            if ~isempty(init_A), Sol_global(1:min(length(init_A),num_A)) = init_A; end
            
            curr_alpha = obj.Relaxation;
            prev_diff = 1e9;
            
            for iter = 1:max_iter
                % 5.1 更新材料属性 (Effective Reluctivity)
                A_curr = Sol_global(1:num_A);
                nu_vec = obj.updateEffectiveNu(space_A, A_curr, matLib);
                
                % 5.2 组装刚度矩阵
                K_AA_loc = obj.Assembler.assembleStiffness(space_A, nu_vec);
                [ki, kj, kv] = find(K_AA_loc);
                K_AA = sparse(ki, kj, kv, numTotalDofs, numTotalDofs);
                
                % 5.3 构建系统矩阵
                % FEM Block: K + jw*M_sigma + eps*M_reg
                ref_scale = max(abs(kv)); if ref_scale==0, ref_scale=1.0; end
                eps_reg = ref_scale * obj.RegScaleRel;
                
                Term_FEM = K_AA + 1j * omega * M_sigma + eps_reg * M_reg;
                
                % Coupling Terms
                % Row A: -W (Circuit I to Field A)
                % Row I: -jw*W' (Field A to Circuit I, EMF)
                [w_idx, ~, w_val] = find(W_vec);
                nz_w = length(w_idx);
                Mat_AI = sparse(w_idx, repmat(numTotalDofs, nz_w, 1), -w_val, numTotalDofs, numTotalDofs);
                Mat_IA = sparse(repmat(numTotalDofs, nz_w, 1), w_idx, -1j * omega * w_val, numTotalDofs, numTotalDofs);
                
                % Circuit Block (I-I): -(R + jw*L)
                Z_ext = circuitProps.R + 1j * omega * circuitProps.L;
                Mat_II = sparse(numTotalDofs, numTotalDofs, -Z_ext, numTotalDofs, numTotalDofs);
                
                % AV Block (Symmetric Modified Potential Form)
                % Row A: + jw*C_AV
                % Row V: + jw*C_AV' + jw*L_VV
                Term_AV = 1j * omega * C_AV;
                Term_VV = 1j * omega * L_VV;
                
                % Total Matrix
                SysMat = Term_FEM + Term_AV + Term_AV.' + Term_VV + Mat_AI + Mat_IA + Mat_II;
                
                % 5.4 右端项 (电压源)
                RHS = sparse(numTotalDofs, 1);
                RHS(end) = -circuitProps.V_source; % Circuit Eq: ... = -V
                
                % 5.5 求解
                [S_bc, R_bc] = BoundaryCondition.applyDirichlet(SysMat, RHS, all_fixed);
                Sol_target = obj.LinearSolver.solve(S_bc, R_bc);
                
                if isempty(Sol_target), warning('Solver failed.'); break; end
                
                % 5.6 自适应松弛更新
                if iter == 1
                    Sol_new = Sol_target;
                    diff_norm = 1.0;
                    curr_alpha = obj.Relaxation;
                else
                    delta = Sol_target - Sol_global;
                    diff_norm = norm(delta) / (norm(Sol_target) + 1e-10);
                    
                    if diff_norm < prev_diff
                        curr_alpha = min(curr_alpha * obj.RelaxationGrow, obj.MaxRelaxation);
                        tag = '>>';
                    else
                        curr_alpha = max(curr_alpha * obj.RelaxationShrink, obj.MinRelaxation);
                        tag = '<<';
                    end
                    Sol_new = Sol_global + curr_alpha * delta;
                    
                    I_mag = abs(Sol_new(end));
                    fprintf('   Iter %2d: RelDiff=%.4e | Alpha=%.2f %s | I=%.4f A\n', ...
                        iter, diff_norm, curr_alpha, tag, I_mag);
                end
                
                if iter > 1 && diff_norm < obj.Tolerance
                    Sol_global = Sol_new;
                    fprintf('   -> Converged.\n');
                    break;
                end
                
                prev_diff = diff_norm;
                Sol_global = Sol_new;
            end
            
            % --- 6. 结果还原 ---
            A_sol = Sol_global(1:num_A);
            I_sol = Sol_global(end);
            if enable_V
                % 还原物理电压 V_phys = jw * V_mod
                V_mod = Sol_global(offset_V+1 : numFemDofs);
                V_sol = 1j * omega * V_mod;
            else
                V_sol = [];
            end
            
            info.Converged = (iter < max_iter) || (diff_norm < obj.Tolerance);
            info.FinalDiff = diff_norm;
        end
        
        function tags = getConductingRegions(~, sigmaMap)
             tags = [];
            if isa(sigmaMap, 'containers.Map')
                keys = sigmaMap.keys;
                for i = 1:length(keys)
                    k = keys{i}; val = sigmaMap(k);
                    if max(abs(val)) > 1e-12
                        if ischar(k), k = str2double(k); end
                        tags = [tags, k]; 
                    end
                end
            elseif isnumeric(sigmaMap)
                tags = find(sigmaMap > 1e-12);
            end
        end
        
        function isLin = checkLinearity(~, matLib)
            isLin = true;
            if isa(matLib, 'containers.Map')
               k = matLib.keys;
               for i=1:length(k), m = matLib(k{i}); if strcmpi(m.Type, 'Nonlinear'), isLin=false; return; end; end
            end
        end
    end
    
    methods (Access = private)
        function nu_vec = updateEffectiveNu(obj, space, A_sol, matLib)
            % 复用 FrequencySolver 的逻辑: 计算 B^2 并查表更新 Nu
            mesh = obj.Assembler.Mesh; T = mesh.T; P = mesh.P;
            numElems = size(T, 2); Tags = mesh.RegionTags;
            nu_vec = zeros(numElems, 1);
            packed = obj.Assembler.preparePackedData(space);
            CellDofs = packed.CellDofs; Signs = packed.Signs;
            [q_pts, ~] = get_quadrature_data('tet', 1);
            [~, curl_ref_raw] = nedelec_tet_p1(q_pts); RefCurl = squeeze(curl_ref_raw)';
            
            for e = 1:numElems
                tag = Tags(e);
                if matLib.isKey(tag), mat = matLib(tag);
                elseif matLib.isKey(num2str(tag)), mat = matLib(num2str(tag));
                else, mat = []; 
                end
                
                if isempty(mat) || strcmpi(mat.Type, 'Linear')
                    if isempty(mat), val = 1/(4*pi*1e-7); else, val = mat.Nu_Linear; end
                    nu_vec(e) = val; continue;
                end
                
                dofs = CellDofs(:, e);
                a_local = A_sol(dofs);
                if isempty(a_local), a_local = zeros(6,1); end
                s_vec = Signs(:, e);
                a_local = a_local .* double(s_vec); 
                
                % 计算 B (RMS 近似)
                % B_real = real(B), B_imag = imag(B)
                % Effective B for BH curve: B_eff = sqrt(|B_x|^2 + |B_y|^2 + |B_z|^2)
                % 或者 Time-Average B^2: (|B|^2)/2 ?
                % 标准做法: Harmonic Balance 用 |B| 峰值, 这里简单用模值
                
                node_ids = T(:, e); p_elem = P(:, node_ids);
                [J_mat, detJ] = compute_jacobian_tet(p_elem);
                b_phasor = (J_mat * (RefCurl * a_local)) / detJ;
                b_mag_sq = norm(b_phasor)^2; % |B|^2
                
                if isfield(mat, 'SplineNu')
                    nu_val = ppval(mat.SplineNu, b_mag_sq);
                else
                    nu_val = mat.nu0;
                end
                nu_vec(e) = nu_val;
            end
        end
    end
end