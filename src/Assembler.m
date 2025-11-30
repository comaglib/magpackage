classdef Assembler < handle
    % ASSEMBLER 有限元组装器 (v6 - MEX Interface)
    
    properties
        Mesh
        DofHandler
        Config
    end
    
    methods
        function obj = Assembler(mesh, dofHandler, config)
            obj.Mesh = mesh;
            obj.DofHandler = dofHandler;
            if nargin < 3
                obj.Config = FemConfig.Default();
            else
                obj.Config = config;
            end
        end
        
        function K = assembleStiffness(obj, space, materialMap)
            packedData = obj.preparePackedData(space);
            elemTags = obj.Mesh.RegionTags;
            try
                Nu_vec = materialMap(elemTags);
            catch
                error('MaterialMap Error.');
            end
            packedData.Nu = Nu_vec(:);
            
            if strcmpi(space.Type, 'Nedelec')
                % ---------------------------------------------------------
                %  自动检测并使用 MEX 加速内核
                % ---------------------------------------------------------
                if exist('assemble_curl_curl_kernel_mex', 'file') == 3
                    % === 路径 A: C++/MEX 加速 ===
                    
                    if isprop(obj.Config, 'DefaultQuadratureOrder')
                        order = obj.Config.DefaultQuadratureOrder;
                    else
                        order = 2; 
                    end
                    
                    [q_pts, q_w] = get_quadrature_data('tet', order); 
                    [~, curl_ref] = nedelec_tet_p1(q_pts); 
                    
                    % [CRITICAL FIX] 显式转换为 double
                    % MEX 接口 mxGetPr 仅支持 double 类型指针。
                    % 实际运行中 Signs 通常为 int8，若不转换会被解析为 NaN。
                    [I, J, V] = assemble_curl_curl_kernel_mex(...
                        packedData.P, ...
                        double(packedData.T), ...        % int32 -> double
                        double(packedData.CellDofs), ... % int32 -> double
                        double(packedData.Signs), ...    % int8 -> double (修复 NaN 源头)
                        packedData.Nu, ...
                        q_w, ...
                        curl_ref);
                        
                else
                    % === 路径 B: MATLAB 并行回退 ===
                    [I, J, V] = assemble_curl_curl_kernel(packedData, obj.Config);
                end
                
                K = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space.');
            end
        end
        
        function M = assembleMass(obj, space)
            packedData = obj.preparePackedData(space);
            if strcmpi(space.Type, 'Nedelec')
                % --- MEX Call ---
                if exist('assemble_mass_kernel_mex', 'file') == 3
                    if isprop(obj.Config, 'DefaultQuadratureOrder')
                        order = obj.Config.DefaultQuadratureOrder;
                    else
                        order = 2; 
                    end
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts); % Note: Mass uses val_ref, not curl_ref
                    
                    [I, J, V] = assemble_mass_kernel_mex(...
                        packedData.P, ...
                        double(packedData.T), ...
                        double(packedData.CellDofs), ...
                        double(packedData.Signs), ...
                        [], ... % Coeff empty
                        q_w, val_ref);
                else
                    [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                end
                
                M = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space.');
            end
        end
        
        function M = assembleMassWeighted(obj, space, coeffMap)
            packedData = obj.preparePackedData(space);
            elemTags = obj.Mesh.RegionTags;
            try
                Coeff_vec = coeffMap(elemTags);
            catch
                error('MaterialMap Error (Mass).');
            end
            packedData.Coeff = Coeff_vec(:);
            
            if strcmpi(space.Type, 'Nedelec')
                % --- MEX Call ---
                if exist('assemble_mass_kernel_mex', 'file') == 3
                    order = obj.Config.DefaultQuadratureOrder;
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts);
                    
                    [I, J, V] = assemble_mass_kernel_mex(...
                        packedData.P, ...
                        double(packedData.T), ...
                        double(packedData.CellDofs), ...
                        double(packedData.Signs), ...
                        packedData.Coeff, ... % Pass Coeff
                        q_w, val_ref);
                else
                    [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                end
                M = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space.');
            end
        end
        
        function F = assembleSource(obj, space, sourceMap)
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            
            % --- MEX Call ---
            if exist('assemble_source_kernel_mex', 'file') == 3
                order = obj.Config.DefaultQuadratureOrder;
                [q_pts, q_w] = get_quadrature_data('tet', order);
                [val_ref, ~] = nedelec_tet_p1(q_pts);
                
                % 注意: assemble_source_kernel_mex 返回的是 I 和 V
                [I_idx, V_val] = assemble_source_kernel_mex(...
                    packedData.P, ...
                    double(packedData.T), ...
                    double(packedData.CellDofs), ...
                    double(packedData.Signs), ...
                    double(packedData.RegionTags), ... % Cast Tags!
                    sourceMap, ...
                    q_w, val_ref);
                
                F = sparse(I_idx, 1, V_val, obj.DofHandler.NumGlobalDofs, 1);
            else
                F = assemble_source_kernel(packedData, sourceMap, obj.Config);
            end
        end
        
        function C = assembleWinding(obj, space, windingObj)
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            
            % --- MEX Call ---
            if exist('assemble_winding_kernel_mex', 'file') == 3
                order = obj.Config.DefaultQuadratureOrder;
                [q_pts, q_w] = get_quadrature_data('tet', order);
                [val_ref, ~] = nedelec_tet_p1(q_pts);
                
                % 准备数据
                if isempty(windingObj.DirectionField)
                    dirField = [];
                else
                    dirField = windingObj.DirectionField;
                end
                
                scale = windingObj.Turns / windingObj.CrossSectionArea;
                
                [I, V] = assemble_winding_kernel_mex(...
                    packedData.P, ...
                    double(packedData.T), ...
                    double(packedData.CellDofs), ...
                    double(packedData.Signs), ...
                    double(packedData.RegionTags), ...
                    double(windingObj.RegionID), ...
                    dirField, ...              % [3 x Ne] or empty
                    windingObj.Direction, ...  % [3 x 1]
                    scale, ...
                    q_w, val_ref);
                    
                C = sparse(I, 1, V, obj.DofHandler.NumGlobalDofs, 1);
            else
                C = assemble_winding_kernel(packedData, windingObj, obj.Config);
            end
        end
        
        function K = assembleScalarLaplacian(obj, space, coeffMap)
            % ASSEMBLESCALARLAPLACIAN 组装标量拉普拉斯矩阵 (e.g., 热传导, 静电场)
            % 支持自动 MEX 加速
            
            packedData = obj.preparePackedData(space);
            
            % 处理材料系数 (Coeff)
            if nargin >= 3 && ~isempty(coeffMap)
                elemTags = obj.Mesh.RegionTags;
                try
                    Coeff_vec = coeffMap(elemTags);
                    packedData.Coeff = Coeff_vec(:);
                catch
                    error('MaterialMap Error: Failed to map region tags to coefficients.');
                end
            else
                packedData.Coeff = [];
            end
            
            % ---------------------------------------------------------
            %  自动检测并使用 MEX 加速内核
            % ---------------------------------------------------------
            if exist('assemble_scalar_laplacian_kernel_mex', 'file') == 3
                % === 路径 A: C++/MEX 加速 ===
                
                % 1. 准备积分数据 (Lagrange P1 梯度为常数，但使用统一接口)
                order = 1; 
                [q_pts, q_w] = get_quadrature_data('tet', order);
                [~, grad_ref] = lagrange_tet_p1(q_pts); % 获取参考梯度 [4 x 3 x Nq]
                
                % 2. 调用 MEX 函数
                % 注意: 必须显式转为 double，特别是 T 和 CellDofs
                [I, J, V] = assemble_scalar_laplacian_kernel_mex(...
                    packedData.P, ...
                    double(packedData.T), ...
                    double(packedData.CellDofs), ...
                    packedData.Coeff, ...
                    q_w, ...
                    grad_ref);
                    
            else
                % === 路径 B: MATLAB 并行回退 ===
                % 调用纯 MATLAB 内核
                [I, J, V] = assemble_scalar_laplacian_kernel(packedData, obj.Config);
            end
            
            % ---------------------------------------------------------
            %  组装稀疏矩阵
            % ---------------------------------------------------------
            K = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
        end
        
        function [J_mat, R_mat] = assembleJacobian(obj, space, solutionA, matLibData, calcJ)
            if nargin < 5, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            
            if strcmpi(space.Type, 'Nedelec')
                
                if exist('assemble_jacobian_kernel_mex', 'file') == 3
                    % === MEX Path ===
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [~, curl_ref] = nedelec_tet_p1(q_pts);
                    
                    % Pack Materials
                    [m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt] = ...
                        obj.packMaterialData(matLibData);
                    
                    % [FIX] 添加额外的 0 作为 NumTags 占位符，满足 MEX 的 17 个输入要求
                    [I, J, V, R_idx, R_val] = assemble_jacobian_kernel_mex(...
                        packedData.P, ...
                        double(packedData.T), ...
                        double(packedData.CellDofs), ...
                        double(packedData.Signs), ...
                        double(packedData.RegionTags), ...
                        solutionA, ...
                        q_w, curl_ref, ...
                        m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt, 0, ... % <--- Added 0 here
                        double(calcJ));
                        
                    if calcJ
                        J_mat = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                    else
                        J_mat = [];
                    end
                    
                    % Residual vector
                    if isempty(R_idx)
                         R_mat = sparse(obj.DofHandler.NumGlobalDofs, 1);
                    else
                         R_mat = sparse(R_idx, 1, R_val, obj.DofHandler.NumGlobalDofs, 1);
                    end
                    
                else
                    % === MATLAB Path ===
                    [I, J, V, R_vec] = assemble_jacobian_kernel(packedData, solutionA, matLibData, obj.Config, calcJ);
                    if calcJ
                        J_mat = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                    else
                        J_mat = [];
                    end
                    R_mat = R_vec; 
                end
            else
                error('Unsupported space.');
            end
        end
        
        function [J_triplets, Res_mat] = assembleHBFEM(obj, space, solHarmonics, aftObj, matLibData, calcJ)
            if nargin < 6, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            
            if strcmpi(space.Type, 'Nedelec')
                if exist('assemble_hbfem_kernel_mex', 'file') == 3
                    % === MEX Path (Complex Fixed) ===
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [~, curl_ref] = nedelec_tet_p1(q_pts);
                    
                    H_vec = aftObj.Harmonics;
                    Scalings = ones(length(H_vec), 1);
                    Scalings(H_vec > 0) = 2.0; 
                    
                    [m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt] = ...
                        obj.packMaterialData(matLibData);
                    
                    % 拆分复数输入
                    solH_full = full(solHarmonics);
                    SolH_R = real(solH_full); SolH_I = imag(solH_full);
                    
                    Dmat_R = real(aftObj.D_matrix); Dmat_I = imag(aftObj.D_matrix);
                    Pmat_R = real(aftObj.P_matrix); Pmat_I = imag(aftObj.P_matrix);

                    % 调用 MEX (23 输入)
                    [I, J, V, R_rows, R_cols, R_val_R, R_val_I] = assemble_hbfem_kernel_mex(...
                        packedData.P, ...
                        double(packedData.T), ...
                        double(packedData.CellDofs), ...
                        double(packedData.Signs), ...
                        double(packedData.RegionTags), ...
                        SolH_R, SolH_I, ... % Split Sol
                        Dmat_R, Dmat_I, ... % Split Dmat
                        Pmat_R, Pmat_I, ... % Split Pmat
                        q_w, curl_ref, ...
                        m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt, 0, ...
                        double(calcJ), ...
                        Scalings);
                    
                    if calcJ
                        % J 矩阵通常是实部近似 (Block Diagonal)，但如果是 full Newton 可能也是复数
                        % 当前 MEX 实现只返回实部 J (nu_dc 近似)，所以不需要 imag
                        J_triplets.I = I; J_triplets.J = J; J_triplets.V = V;
                    else
                        J_triplets = [];
                    end
                    
                    if isempty(R_rows)
                        Res_mat = sparse(obj.DofHandler.NumGlobalDofs, aftObj.NumHarmonics);
                    else
                        % 重建复数残差
                        R_val_Complex = complex(R_val_R, R_val_I);
                        Res_mat = sparse(R_rows, R_cols, R_val_Complex, obj.DofHandler.NumGlobalDofs, aftObj.NumHarmonics);
                    end
                    
                else
                    % === MATLAB Path ===
                    [I, J, V, R_mat] = assemble_hbfem_kernel(packedData, solHarmonics, aftObj, matLibData, obj.Config, calcJ);
                    if calcJ
                        J_triplets.I = I; J_triplets.J = J; J_triplets.V = V;
                    else
                        J_triplets = [];
                    end
                    Res_mat = R_mat;
                end
            else
                error('Unsupported space.');
            end
        end
        
        function packedData = preparePackedData(obj, space)
            packedData = obj.DofHandler.packForKernel(space);
            packedData.P = obj.Mesh.P;
            packedData.T = obj.Mesh.T;
        end
    end
    
    methods (Access = private)
        function [linNu, isNon, maxB, breaks, coefs, starts, counts] = packMaterialData(~, matLibData)
            % PACKMATERIALDATA 将材料数据打包为 MEX 友好的扁平数组 (优化版)
            % 1. 移除了冗余的 if-else。
            % 2. 增加了 breaks 和 coefs 的预分配逻辑。
            
            % --- 1. 确定遍历范围 ---
            if isa(matLibData, 'containers.Map')
                loop_keys = cell2mat(matLibData.keys);
                maxTag = max(loop_keys);
            else
                % 假设结构体数组的索引即为 Tag
                loop_keys = 1:length(matLibData);
                maxTag = length(matLibData);
            end
            
            % --- 2. 预分配固定大小数组 ---
            % 索引直接对应 Tag，因此大小为 maxTag + 1
            linNu  = zeros(maxTag+1, 1);
            isNon  = zeros(maxTag+1, 1);
            maxB   = zeros(maxTag+1, 1);
            starts = zeros(maxTag+1, 1);
            counts = zeros(maxTag+1, 1);
            
            % --- 3. 第一遍遍历：计算动态数组所需的总大小 ---
            total_breaks_len = 0;
            total_coefs_len = 0;
            
            for tag = loop_keys
                if tag <= 0, continue; end
                
                % 通用获取方式
                mat = matLibData(tag);
                
                if strcmp(mat.Type, 'Nonlinear')
                    pp = mat.SplineNu;
                    % breaks 长度
                    total_breaks_len = total_breaks_len + length(pp.breaks);
                    % coefs 总元素数 (Pieces * Order)
                    total_coefs_len = total_coefs_len + numel(pp.coefs);
                end
            end
            
            % --- 4. 预分配动态数组 ---
            breaks = zeros(total_breaks_len, 1);
            coefs  = zeros(total_coefs_len, 1);
            
            % --- 5. 第二遍遍历：填充数据 ---
            cpp_idx_ptr = 0; % 供 C++ 使用的 0-based 偏移量
            brk_ptr = 1;     % MATLAB 数组填充指针 (1-based)
            coe_ptr = 1;     % MATLAB 数组填充指针 (1-based)
            
            for tag = loop_keys
                if tag <= 0, continue; end
                
                mat = matLibData(tag);
                
                if strcmp(mat.Type, 'Linear')
                    linNu(tag) = mat.Nu_Linear;
                    isNon(tag) = 0;
                    % 线性材料不需要设置 starts/counts/maxB，保持为 0 即可
                else
                    linNu(tag) = mat.nu0; % 默认值或饱和值
                    isNon(tag) = 1;
                    
                    if isfield(mat, 'MaxBSq')
                        maxB(tag) = mat.MaxBSq; 
                    else
                        maxB(tag) = 1e10; 
                    end
                    
                    pp = mat.SplineNu;
                    b_vec = pp.breaks(:);
                    c_mat = pp.coefs; % [Pieces x 4]
                    
                    n_b = length(b_vec);
                    n_c = numel(c_mat);
                    n_pieces = n_b - 1;
                    
                    % 记录起始位置 (C++ 侧索引)
                    starts(tag) = cpp_idx_ptr;
                    counts(tag) = n_pieces;
                    
                    % 填充 breaks
                    breaks(brk_ptr : brk_ptr + n_b - 1) = b_vec;
                    
                    % 填充 coefs 
                    % 注意：pp.coefs 是 [Pieces x Order]，我们需要按行展平
                    % c_mat' 转置后变成 [Order x Pieces]，按列存储即为 [row1, row2...]
                    c_linear = c_mat'; 
                    coefs(coe_ptr : coe_ptr + n_c - 1) = c_linear(:);
                    
                    % 更新指针
                    brk_ptr = brk_ptr + n_b;
                    coe_ptr = coe_ptr + n_c;
                    cpp_idx_ptr = cpp_idx_ptr + n_b; % C++ breaks 指针前移
                end
            end
        end
    end
end