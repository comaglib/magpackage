classdef Assembler < handle
    % ASSEMBLER 有限元组装器 (v10.1 - Full MEX Support)
    %
    % 核心职责:
    %   负责协调网格(Mesh)、自由度(DofHandler)与底层积分内核(Kernels)，
    %   将物理方程的弱形式转化为全局稀疏矩阵或向量。
    %
    % 特性:
    %   1. 混合内核架构 (Hybrid Kernel Architecture): 
    %      优先调用 C++/MEX 高性能内核，若缺失则自动回退到 MATLAB 原生实现。
    %   2. 单元级材料支持: 支持传入非均匀材料属性 (Element-wise Material)。
    %   3. 向量化组装: 针对各类物理场 (Stiffness, Mass, Source, Coupling) 提供统一接口。
    
    properties
        Mesh        % 网格对象 (包含节点 P,拓扑 T,区域 RegionTags)
        DofHandler  % 自由度处理器 (管理全局/局部 DoF 映射)
        Config      % 配置对象 (积分阶数, 并行设置等)
    end
    
    methods
        function obj = Assembler(mesh, dofHandler, config)
            % 构造函数
            obj.Mesh = mesh;
            obj.DofHandler = dofHandler;
            if nargin < 3
                obj.Config = FemConfig.Default();
            else
                obj.Config = config;
            end
        end
        
        %  1. 刚度矩阵组装 (Stiffness Matrix) - CurlCurl
        function K = assembleStiffness(obj, space, materialMap)
            % ASSEMBLESTIFFNESS 组装旋度-旋度矩阵 S = (nu * curl N_i, curl N_j)
            %
            % 输入:
            %   space       - 有限元空间对象 (FunctionSpace)
            %   materialMap - 材料属性映射 (Map 或 向量)
            
            % 1. 准备组装数据
            packedData = obj.preparePackedData(space);
            
            % 2. 扩展材料属性 (支持 Map 映射或直接传入单元向量)
            packedData.Nu = obj.expandMaterialMap(materialMap, obj.Mesh.RegionTags);
            
            % 3. 安全性检查 (确保 DoF 映射完整)
            canUseMex = obj.checkMexSafety(packedData);
            
            if strcmpi(space.Type, 'Nedelec')
                % 优先尝试 C++/MEX 内核
                if canUseMex && exist('assemble_curl_curl_kernel_mex', 'file') == 3
                    % 获取积分参数
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order); 
                    [~, curl_ref] = nedelec_tet_p1(q_pts); % 参考单元旋度值
                    
                    % 调用 MEX
                    [I, J, V] = assemble_curl_curl_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), ...
                        double(packedData.Signs), packedData.Nu, q_w, curl_ref);
                else
                    % 回退到 MATLAB 内核
                    [I, J, V] = assemble_curl_curl_kernel(packedData, obj.Config);
                end
                
                % 构造稀疏矩阵
                valid = (I > 0) & (J > 0);
                K = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space type: %s', space.Type);
            end
        end
        
        %  2. 质量矩阵组装 (Mass Matrix)
        function M = assembleMass(obj, space)
            % ASSEMBLEMASS 组装标准质量矩阵 M = (N_i, N_j)
            
            packedData = obj.preparePackedData(space);
            canUseMex = obj.checkMexSafety(packedData);
            
            if strcmpi(space.Type, 'Nedelec')
                if canUseMex && exist('assemble_mass_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts); 
                    
                    [I, J, V] = assemble_mass_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), ...
                        double(packedData.Signs), [], q_w, val_ref);
                else
                    [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                end
                
                valid = (I > 0) & (J > 0);
                M = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space type: %s', space.Type);
            end
        end
        
        function M = assembleMassWeighted(obj, space, coeffMap)
            % ASSEMBLEMASSWEIGHTED 组装加权质量矩阵 M = (sigma * N_i, N_j)
            % 常用于涡流项或非均匀介质
            
            packedData = obj.preparePackedData(space);
            packedData.Coeff = obj.expandMaterialMap(coeffMap, obj.Mesh.RegionTags);
            canUseMex = obj.checkMexSafety(packedData);
            
            if strcmpi(space.Type, 'Nedelec')
                if canUseMex && exist('assemble_mass_kernel_mex', 'file') == 3
                    order = obj.Config.DefaultQuadratureOrder;
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts);
                    
                    % 传入系数向量 Coeff
                    [I, J, V] = assemble_mass_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), ...
                        double(packedData.Signs), packedData.Coeff, q_w, val_ref);
                else
                    [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                end
                
                valid = (I > 0) & (J > 0);
                M = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space type: %s', space.Type);
            end
        end
        
        %  3. 源向量组装 (Source Vector)
        function F = assembleSource(obj, space, sourceMap)
            % ASSEMBLESOURCE 组装源向量 F_i = (J_s, N_i)
            % 
            % 更新: 优先调用 assemble_source_kernel_mex
            
            % 1. 准备数据
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            packedData.Js = obj.expandSourceMap(sourceMap, obj.Mesh.RegionTags);
            
            % 2. 检查 MEX 可用性
            canUseMex = obj.checkMexSafety(packedData);
            
            if strcmpi(space.Type, 'Nedelec')
                if canUseMex && exist('assemble_source_kernel_mex', 'file') == 3
                    % --- A. MEX 路径 ---
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts); % Nedelec 基函数参考值
                    
                    % 调用 MEX 内核
                    % Input: P, T, CellDofs, Signs, Js, Qw, ValRef
                    [I, V] = assemble_source_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), ...
                        double(packedData.Signs), packedData.Js, q_w, val_ref);
                    
                    % 组装稀疏向量 (accumarray 逻辑)
                    valid = (I > 0);
                    F = sparse(I(valid), 1, V(valid), obj.DofHandler.NumGlobalDofs, 1);
                    
                else
                    % --- B. MATLAB 回退 ---
                    F_full = assemble_source_kernel(packedData, obj.Config);
                    F = F_full;
                    
                    % 调整大小 (MATLAB Kernel 可能返回压缩后的向量)
                    if size(F, 1) < obj.DofHandler.NumGlobalDofs
                        F(obj.DofHandler.NumGlobalDofs, 1) = 0;
                    end
                end
            else
                error('Unsupported space for Source Assembly: %s', space.Type);
            end
        end
        
        %  4. 场路耦合矩阵 (Coupling Matrix)
        function C = assembleCoupling(obj, spaceRow, spaceCol, coeffMap)
            % ASSEMBLECOUPLING 组装耦合矩阵 C_ij = (sigma * N_i, grad phi_j)
            % 通常用于 A-V 耦合系统
            % 
            % 更新: 优先调用 assemble_coupling_kernel_mex
            
            % 1. 准备数据
            packedRow = obj.preparePackedData(spaceRow); % Edge Space (A)
            packedCol = obj.preparePackedData(spaceCol); % Nodal Space (V)
            
            % 扩展系数 (Sigma)
            Sigma_vec = obj.expandMaterialMap(coeffMap, obj.Mesh.RegionTags);
            
            % 2. 检查 MEX 可用性 (需同时检查行列空间)
            canUseMex = obj.checkMexSafety(packedRow) && obj.checkMexSafety(packedCol);
            
            if strcmpi(spaceRow.Type, 'Nedelec') && strcmpi(spaceCol.Type, 'Lagrange')
                if canUseMex && exist('assemble_coupling_kernel_mex', 'file') == 3
                    % --- A. MEX 路径 ---
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [ned_ref, ~] = nedelec_tet_p1(q_pts);   % Nedelec 值
                    [~, grad_ref] = lagrange_tet_p1(q_pts); % Lagrange 梯度
                    
                    % 调用 MEX 内核
                    % Input: P, T, DofsE, Signs, DofsN, Coeffs, Qw, NedRef, GradRef
                    [I, J, V] = assemble_coupling_kernel_mex(...
                        packedRow.P, double(packedRow.T), double(packedRow.CellDofs), double(packedRow.Signs), ...
                        double(packedCol.CellDofs), Sigma_vec, q_w, ned_ref, grad_ref);
                    
                    valid = (I > 0) & (J > 0);
                    C = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                    
                else
                    % --- B. MATLAB 回退 ---
                    [I, J, V] = assemble_coupling_kernel(...
                        packedRow.P, packedRow.T, packedRow.CellDofs, packedRow.Signs, ...
                        packedCol.CellDofs, Sigma_vec, obj.Config);
                    
                    valid = (I > 0) & (J > 0);
                    C = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                end
            else
                error('Unsupported space combination for Coupling Assembly.');
            end
        end
        
        %  5. 标量场 Laplacian (Scalar Potential)
        function K = assembleScalarLaplacian(obj, space, coeffMap)
            % ASSEMBLESCALARLAPLACIAN 组装标量扩散矩阵 K = (sigma * grad phi_i, grad phi_j)
            
            packedData = obj.preparePackedData(space);
            packedData.Coeff = obj.expandMaterialMap(coeffMap, obj.Mesh.RegionTags);
            canUseMex = obj.checkMexSafety(packedData);
            
            if canUseMex && exist('assemble_scalar_laplacian_kernel_mex', 'file') == 3
                order = 1; % 标量场通常 P1 足够
                [q_pts, q_w] = get_quadrature_data('tet', order);
                [~, grad_ref] = lagrange_tet_p1(q_pts); 
                
                [I, J, V] = assemble_scalar_laplacian_kernel_mex(...
                    packedData.P, double(packedData.T), double(packedData.CellDofs), ...
                    packedData.Coeff, q_w, grad_ref);
            else
                [I, J, V] = assemble_scalar_laplacian_kernel(packedData, obj.Config);
            end
            
            valid = (I > 0) & (J > 0);
            K = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
        end
        
        %  6. 线圈源向量 (Winding)
        function C = assembleWinding(obj, space, windingObj)
            % ASSEMBLEWINDING 组装线圈电流投影向量 (绞线模型)
            
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            canUseMex = obj.checkMexSafety(packedData);
            
            if canUseMex && exist('assemble_winding_kernel_mex', 'file') == 3
                order = obj.Config.DefaultQuadratureOrder;
                [q_pts, q_w] = get_quadrature_data('tet', order);
                [val_ref, ~] = nedelec_tet_p1(q_pts);
                
                if isempty(windingObj.DirectionField), dirField = []; else, dirField = windingObj.DirectionField; end
                scale = windingObj.Turns / windingObj.CrossSectionArea;
                
                % 调用专门的 Winding MEX
                [I, V] = assemble_winding_kernel_mex(...
                    packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), ...
                    double(packedData.RegionTags), double(windingObj.RegionID), ...
                    dirField, windingObj.Direction, scale, q_w, val_ref);
                
                valid = (I > 0);
                C = sparse(I(valid), 1, V(valid), obj.DofHandler.NumGlobalDofs, 1);
            else
                C = assemble_winding_kernel(packedData, windingObj, obj.Config);
            end
        end
        
        %  7. Jacobian 矩阵 (Nonlinear Newton)
        function [J_mat, R_mat] = assembleJacobian(obj, space, solutionA, matLibData, calcJ)
            % ASSEMBLEJACOBIAN 组装非线性切线刚度矩阵与残差
            
            if nargin < 5, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            canUseMex = obj.checkMexSafety(packedData);
            
            if strcmpi(space.Type, 'Nedelec')
                if canUseMex && exist('assemble_jacobian_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [~, curl_ref] = nedelec_tet_p1(q_pts);
                    
                    % 打包材料数据用于 MEX
                    [m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt] = obj.packMaterialData(matLibData);
                    
                    [I, J, V, R_idx, R_val] = assemble_jacobian_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), ...
                        double(packedData.RegionTags), solutionA, q_w, curl_ref, ...
                        m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt, 0, double(calcJ));
                    
                    if calcJ
                        valid = (I > 0) & (J > 0); 
                        J_mat = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs); 
                    else
                        J_mat = []; 
                    end
                    
                    if isempty(R_idx)
                        R_mat = sparse(obj.DofHandler.NumGlobalDofs, 1); 
                    else
                        valid_r = (R_idx > 0); 
                        R_mat = sparse(R_idx(valid_r), 1, R_val(valid_r), obj.DofHandler.NumGlobalDofs, 1); 
                    end
                else
                    % MATLAB Fallback
                    [I, J, V, R_vec] = assemble_jacobian_kernel(packedData, solutionA, matLibData, obj.Config, calcJ);
                    if calcJ
                        valid = (I > 0) & (J > 0); 
                        J_mat = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs); 
                    else
                        J_mat = []; 
                    end
                    R_mat = R_vec; 
                end
            else
                error('Unsupported space type.');
            end
        end
        
        %  8. 谐波平衡矩阵 (HBFEM)
        function [J_triplets, Res_mat] = assembleHBFEM(obj, space, solHarmonics, aftObj, matLibData, calcJ)
            % ASSEMBLEHBFEM 组装谐波平衡 Jacobian 与残差
            
            if nargin < 6, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            canUseMex = obj.checkMexSafety(packedData);
            
            if strcmpi(space.Type, 'Nedelec')
                if canUseMex && exist('assemble_hbfem_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [~, curl_ref] = nedelec_tet_p1(q_pts);
                    
                    % 准备 HBFEM 专用数据
                    H_vec = aftObj.Harmonics; 
                    Scalings = ones(length(H_vec), 1); 
                    Scalings(H_vec > 0) = 2.0; 
                    [m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt] = obj.packMaterialData(matLibData);
                    
                    solH_full = full(solHarmonics); 
                    SolH_R = real(solH_full); SolH_I = imag(solH_full);
                    Dmat_R = real(aftObj.D_matrix); Dmat_I = imag(aftObj.D_matrix); 
                    Pmat_R = real(aftObj.P_matrix); Pmat_I = imag(aftObj.P_matrix);
                    
                    [I, J, V, R_rows, R_cols, R_val_R, R_val_I] = assemble_hbfem_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), ...
                        double(packedData.RegionTags), SolH_R, SolH_I, Dmat_R, Dmat_I, Pmat_R, Pmat_I, ...
                        q_w, curl_ref, m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt, 0, double(calcJ), Scalings);
                    
                    if calcJ
                        valid = (I > 0) & (J > 0); 
                        J_triplets.I = I(valid); J_triplets.J = J(valid); J_triplets.V = V(valid); 
                    else
                        J_triplets = []; 
                    end
                    
                    if isempty(R_rows)
                        Res_mat = sparse(obj.DofHandler.NumGlobalDofs, aftObj.NumHarmonics); 
                    else
                        valid_r = (R_rows > 0) & (R_cols > 0); 
                        R_val_Complex = complex(R_val_R, R_val_I); 
                        Res_mat = sparse(R_rows(valid_r), R_cols(valid_r), R_val_Complex(valid_r), obj.DofHandler.NumGlobalDofs, aftObj.NumHarmonics); 
                    end
                else
                    % MATLAB Fallback
                    [I, J, V, R_mat] = assemble_hbfem_kernel(packedData, solHarmonics, aftObj, matLibData, obj.Config, calcJ);
                    if calcJ
                        valid = (I > 0) & (J > 0); 
                        J_triplets.I = I(valid); J_triplets.J = J(valid); J_triplets.V = V(valid); 
                    else
                        J_triplets = []; 
                    end
                    Res_mat = R_mat;
                end
            else
                error('Unsupported space type.');
            end
        end
        
        %  Helpers (辅助方法)
        function packedData = preparePackedData(obj, space)
            % 提取几何与拓扑信息，打包为结构体
            packedData = obj.DofHandler.packForKernel(space);
            packedData.P = obj.Mesh.P;
            packedData.T = obj.Mesh.T;
        end
        
        function [J_mats, detJs, curl_ref] = precomputeGeometryData(obj)
            % 预计算几何数据 (仅用于调试或旧接口兼容)
            order = max(obj.Config.DefaultQuadratureOrder, 2);
            [q_pts, ~] = get_quadrature_data('tet', order);
            [~, curl_ref] = nedelec_tet_p1(q_pts);
            J_mats = []; detJs = [];
        end
        
        function vec = expandMaterialMap(~, mapOrArray, tags)
            % EXPANDMATERIALMAP 将材料映射或数组扩展为单元级向量
            % 1. Map: 根据 RegionTag 查找
            % 2. 标量: 广播到所有单元
            % 3. 向量 (长度==ElemNum): 直接使用 (Element-wise)
            
            numElems = length(tags);
            
            if isa(mapOrArray, 'containers.Map')
                vec = zeros(numElems, 1); 
                keys = mapOrArray.keys;
                for i = 1:length(keys)
                    k = keys{i}; 
                    if ischar(k), k_val = str2double(k); else, k_val = k; end
                    val = mapOrArray(k); 
                    if isscalar(val)
                        vec(tags == k_val) = val; 
                    end
                end
            elseif isnumeric(mapOrArray)
                if isscalar(mapOrArray)
                    vec = ones(numElems, 1) * mapOrArray;
                elseif isvector(mapOrArray)
                    if numel(mapOrArray) == numElems
                        vec = mapOrArray(:);
                    else
                        vec = mapOrArray(tags); % 尝试作为 Lookup Table
                    end
                else
                    error('Assembler:InvalidInput', 'Invalid material array format.');
                end
            else
                error('Assembler:InvalidInput', 'Material map must be a Map or numeric array.');
            end
        end
    end
    
    methods (Access = private)
        function safe = checkMexSafety(~, packedData)
            % 检查数据完整性，确保传入 MEX 不会崩溃
            % 主要是检查 CellDofs 是否含有 0 索引 (未分配自由度的单元)
            if any(packedData.CellDofs(:) == 0)
                safe = false;
            else
                safe = true;
            end
        end
        
        function mat = expandSourceMap(~, mapOrArray, tags)
            % 扩展源项映射为 [3 x NumElems] 矩阵
            numElems = length(tags); mat = zeros(3, numElems);
            if isa(mapOrArray, 'containers.Map')
                keys = mapOrArray.keys;
                for i = 1:length(keys)
                    k = keys{i}; if ischar(k), k_val = str2double(k); else, k_val = k; end
                    val = mapOrArray(k); mask = (tags == k_val);
                    if any(mask), mat(1, mask) = val(1); mat(2, mask) = val(2); mat(3, mask) = val(3); end
                end
            elseif isnumeric(mapOrArray)
                 if numel(mapOrArray) == 3, mat(1,:) = mapOrArray(1); mat(2,:) = mapOrArray(2); mat(3,:) = mapOrArray(3);
                 elseif size(mapOrArray, 1) == 3 && size(mapOrArray, 2) == numElems, mat = mapOrArray;
                 else, error('Source must be [3x1] or [3xN].'); 
                 end
            end
        end
        
        function [linNu, isNon, maxB, breaks, coefs, starts, counts] = packMaterialData(~, matLibData)
            % 将复杂的 MaterialLib 对象数组打包为 C++ 可读的扁平数组 (SoA Layout)
            if isa(matLibData, 'containers.Map'), loop_keys = cell2mat(matLibData.keys); maxTag = max(loop_keys);
            else, loop_keys = 1:length(matLibData); maxTag = length(matLibData); end
            
            linNu=zeros(maxTag+1,1); isNon=zeros(maxTag+1,1); maxB=zeros(maxTag+1,1); 
            starts=zeros(maxTag+1,1); counts=zeros(maxTag+1,1);
            
            total_breaks_len=0; total_coefs_len=0;
            for tag=loop_keys
                if tag<=0, continue; end; mat=matLibData(tag);
                if strcmp(mat.Type,'Nonlinear'), pp=mat.SplineNu; total_breaks_len=total_breaks_len+length(pp.breaks); total_coefs_len=total_coefs_len+numel(pp.coefs); end
            end
            breaks=zeros(total_breaks_len,1); coefs=zeros(total_coefs_len,1);
            
            cpp_idx_ptr=0; brk_ptr=1; coe_ptr=1;
            for tag=loop_keys
                if tag<=0, continue; end; mat=matLibData(tag);
                if strcmp(mat.Type,'Linear')
                    linNu(tag)=mat.Nu_Linear; isNon(tag)=0;
                else
                    linNu(tag)=mat.nu0; isNon(tag)=1; 
                    if isfield(mat,'MaxBSq'), maxB(tag)=mat.MaxBSq; else, maxB(tag)=1e10; end
                    pp=mat.SplineNu; b_vec=pp.breaks(:); c_mat=pp.coefs; n_b=length(b_vec); n_c=numel(c_mat);
                    starts(tag)=cpp_idx_ptr; counts(tag)=length(b_vec)-1;
                    breaks(brk_ptr:brk_ptr+n_b-1)=b_vec;
                    c_linear=c_mat'; coefs(coe_ptr:coe_ptr+n_c-1)=c_linear(:);
                    brk_ptr=brk_ptr+n_b; coe_ptr=coe_ptr+n_c; cpp_idx_ptr=cpp_idx_ptr+n_b;
                end
            end
        end
    end
end