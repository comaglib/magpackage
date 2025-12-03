classdef FrequencySolver < handle
    % FREQUENCYSOLVER 频域时谐磁场求解器 (v10.4 - Symmetric Optimized)
    %
    % 改进:
    %   1. [Optimization] 启用线性求解器对称模式 (MumpsSymmetry = 2)。
    %      针对复数对称矩阵 (Complex Symmetric) 优化，减少内存和计算时间。
    
    properties
        Assembler
        LinearSolver
        Frequency
        MaxIterations = 200
        Tolerance = 1e-5
        
        Relaxation      = 0.5
        MinRelaxation   = 0.05
        MaxRelaxation   = 1.0
        RelaxationGrow  = 1.05
        RelaxationShrink= 0.6
    end
    
    methods
        function obj = FrequencySolver(assembler, freq)
            obj.Assembler = assembler;
            obj.Frequency = freq;
            obj.LinearSolver = LinearSolver('Auto');
            
            warning('off', 'MATLAB:rankDeficientMatrix'); 
            warning('off', 'MATLAB:nearlySingularMatrix');
            
            % [MUMPS 设置]
            % i14 = 100: 增加 Pivot 阈值以提高数值稳定性
            obj.LinearSolver.MumpsICNTL.i14 = 100; 
            
            % [New] 启用对称模式
            % SYM = 2: General Symmetric (适用于复数对称矩阵 A.' == A)
            % 注意: 不能用 1 (SPD)，因为复数系统不是 Hermitian 正定的
            obj.LinearSolver.MumpsSymmetry = 2;
        end
        
        function [A_sol, V_sol] = solve(obj, space_A, space_V, matLib, sigmaMap, sourceMap, fixedDofs_A, fixedDofs_V)
            
            if nargin < 8, fixedDofs_V = []; end
            
            omega = 2 * pi * obj.Frequency;
            if omega == 0, omega = 1e-9; end 
            
            fprintf('==============================================\n');
            fprintf('   Frequency Solver (%.1f Hz) - Symmetric Mode\n', obj.Frequency);
            fprintf('==============================================\n');
            
            dofHandler = obj.Assembler.DofHandler;
            [~, offset_A] = dofHandler.distributeDofs(space_A);
            num_A = dofHandler.SpaceLocalSizes(space_A.toString());
            
            enable_V = ~isempty(space_V);
            offset_V = -1;
            auto_ground_idx = [];
            
            if enable_V
                conductingTags = obj.getConductingRegions(sigmaMap);
                if isempty(conductingTags)
                    enable_V = false;
                else
                    if ~dofHandler.DofMaps.isKey(space_V.toString())
                        [~, offset_V] = dofHandler.distributeDofs(space_V, conductingTags);
                    else
                        offset_V = dofHandler.SpaceOffsets(space_V.toString());
                    end
                    
                    hasFixedV = ~isempty(fixedDofs_V) && any(fixedDofs_V);
                    if ~hasFixedV
                        fprintf('   [Auto-Fix] V-field floating. Grounding 1st node.\n');
                        auto_ground_idx = offset_V + 1;
                    end
                end
            end
            numTotalDofs = dofHandler.NumGlobalDofs;
            
            fprintf('   [Assembly] Pre-computing matrices...\n');
            M_sigma = sparse(numTotalDofs, numTotalDofs);
            C_AV    = sparse(numTotalDofs, numTotalDofs);
            L_VV    = sparse(numTotalDofs, numTotalDofs);
            M_reg   = obj.Assembler.assembleMass(space_A);
            
            if enable_V
                M_sigma = obj.Assembler.assembleMassWeighted(space_A, sigmaMap);
                C_AV    = obj.Assembler.assembleCoupling(space_A, space_V, sigmaMap);
                L_VV    = obj.Assembler.assembleScalarLaplacian(space_V, sigmaMap);
            end
            
            F_A = obj.Assembler.assembleSource(space_A, sourceMap);
            F_sys = sparse(numTotalDofs, 1);
            if ~isempty(F_A), F_sys(1:length(F_A)) = F_A; end
            
            if islogical(fixedDofs_A), idx_A = find(fixedDofs_A); else, idx_A = fixedDofs_A(:); end
            all_fixed = idx_A;
            if enable_V
                 if isempty(fixedDofs_V), idx_V = [];
                 elseif islogical(fixedDofs_V), idx_V = find(fixedDofs_V);
                 else, idx_V = fixedDofs_V(:); 
                 end
                 if ~isempty(auto_ground_idx), idx_V = unique([idx_V; auto_ground_idx]); end
                 idx_V = idx_V(idx_V <= numTotalDofs);
                 all_fixed = [all_fixed; idx_V];
            end
            
            hasNonlinear = false;
            if isa(matLib, 'containers.Map')
                keys = matLib.keys;
                for k = 1:length(keys)
                    mat = matLib(keys{k});
                    if isfield(mat, 'Type') && strcmpi(mat.Type, 'Nonlinear'), hasNonlinear = true; break; end
                end
            end
            max_iter = hasNonlinear * (obj.MaxIterations - 1) + 1;
            
            Sol_global = zeros(numTotalDofs, 1);
            A_curr = zeros(num_A, 1);
            prev_diff = 1e9;
            curr_alpha = obj.Relaxation;
            
            for iter = 1:max_iter
                nu_vec = obj.updateEffectiveNu(space_A, A_curr, matLib);
                K_AA = obj.Assembler.assembleStiffness(space_A, nu_vec);
                
                ref_scale = max(abs(diag(K_AA)));
                if ref_scale == 0, ref_scale = 1.0; end
                eps_reg = ref_scale * 1e-6; 
                
                Term_A = K_AA + 1j * omega * M_sigma + eps_reg * M_reg;
                Term_C = 1j * omega * C_AV; 
                Term_L = 1j * omega * L_VV;
                
                % [Note] 这里构造的是严格对称的矩阵
                SysMat = Term_A + Term_C + Term_C.' + Term_L;
                
                [S_sys, R_sys] = BoundaryCondition.applyDirichlet(SysMat, F_sys, all_fixed);
                
                % LinearSolver 将利用 MumpsSymmetry = 2 进行加速
                Sol_target = obj.LinearSolver.solve(S_sys, R_sys);
                
                if isempty(Sol_target), warning('Solver failed.'); break; end
                
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
                    fprintf('   Iter %2d: RelDiff = %.4e | Alpha = %.3f %s\n', iter, diff_norm, curr_alpha, tag);
                end
                
                if iter > 1 && diff_norm < obj.Tolerance
                    fprintf('   -> Converged.\n');
                    Sol_global = Sol_new;
                    break;
                end
                
                prev_diff = diff_norm;
                Sol_global = Sol_new;
                A_curr = Sol_global(1:num_A);
            end
            
            A_sol = Sol_global(1:num_A);
            if enable_V
                V_mod = Sol_global(offset_V+1 : end);
                V_sol = 1j * omega * V_mod; 
            else
                V_sol = [];
            end
            fprintf('==============================================\n');
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
    end
    
    methods (Access = private)
        function nu_vec = updateEffectiveNu(obj, space, A_sol, matLib)
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
                node_ids = T(:, e); p_elem = P(:, node_ids);
                [J_mat, detJ] = compute_jacobian_tet(p_elem);
                b_vec = (J_mat * (RefCurl * a_local)) / detJ;
                b_mag = norm(b_vec); 
                
                if isfield(mat, 'SplineNu')
                    nu_val = ppval(mat.SplineNu, b_mag^2);
                else
                    nu_val = mat.nu0;
                end
                nu_vec(e) = nu_val;
            end
        end
    end
end