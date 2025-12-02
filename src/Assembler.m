classdef Assembler < handle
    % ASSEMBLER 有限元组装器 (v10.5 - Data Integrity Check)
    
    properties
        Mesh        % 网格对象
        DofHandler  % 自由度处理器
        Config      % 配置对象
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
        
        % 1. assembleStiffness
        function K = assembleStiffness(obj, space, materialMap)
            packedData = obj.preparePackedData(space);
            packedData.Nu = obj.expandMaterialMap(materialMap, obj.Mesh.RegionTags);
            if exist('assemble_curl_curl_kernel_mex', 'file') == 3
                order = max(obj.Config.DefaultQuadratureOrder, 2);
                [q_pts, q_w] = get_quadrature_data('tet', order); 
                [~, curl_ref] = nedelec_tet_p1(q_pts);
                [I, J, V] = assemble_curl_curl_kernel_mex(packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), packedData.Nu, q_w, curl_ref);
            else
                [I, J, V] = assemble_curl_curl_kernel(packedData, obj.Config);
            end
            valid = (I > 0) & (J > 0);
            K = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
        end
        
        % 2. assembleMass
        function M = assembleMass(obj, space)
            packedData = obj.preparePackedData(space);
            if strcmpi(space.Type, 'Nedelec')
                if exist('assemble_mass_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts); 
                    [I, J, V] = assemble_mass_kernel_mex(packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), [], q_w, val_ref);
                else
                    [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                end
                valid = (I > 0) & (J > 0);
                M = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space type');
            end
        end
        
        % 3. assembleMassWeighted
        function M = assembleMassWeighted(obj, space, coeffMap)
            packedData = obj.preparePackedData(space);
            packedData.Coeff = obj.expandMaterialMap(coeffMap, obj.Mesh.RegionTags);
            if strcmpi(space.Type, 'Nedelec')
                if exist('assemble_mass_kernel_mex', 'file') == 3
                    order = obj.Config.DefaultQuadratureOrder;
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts);
                    [I, J, V] = assemble_mass_kernel_mex(packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), packedData.Coeff, q_w, val_ref);
                else
                    [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                end
                valid = (I > 0) & (J > 0);
                M = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space type');
            end
        end
        
        % 4. assembleSource
        function F = assembleSource(obj, space, sourceMap)
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            packedData.Js = obj.expandSourceMap(sourceMap, obj.Mesh.RegionTags);
            if strcmpi(space.Type, 'Nedelec')
                if exist('assemble_source_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [val_ref, ~] = nedelec_tet_p1(q_pts);
                    [I, V] = assemble_source_kernel_mex(packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), packedData.Js, q_w, val_ref);
                    valid = (I > 0);
                    F = sparse(I(valid), 1, V(valid), obj.DofHandler.NumGlobalDofs, 1);
                else
                    F = assemble_source_kernel(packedData, obj.Config);
                    if size(F, 1) < obj.DofHandler.NumGlobalDofs, F(obj.DofHandler.NumGlobalDofs, 1) = 0; end
                end
            else
                error('Unsupported space');
            end
        end
        
        % 5. assembleCoupling
        function C = assembleCoupling(obj, spaceRow, spaceCol, coeffMap)
            packedRow = obj.preparePackedData(spaceRow); 
            packedCol = obj.preparePackedData(spaceCol); 
            Sigma_vec = obj.expandMaterialMap(coeffMap, obj.Mesh.RegionTags);
            if strcmpi(spaceRow.Type, 'Nedelec') && strcmpi(spaceCol.Type, 'Lagrange')
                if exist('assemble_coupling_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [ned_ref, ~] = nedelec_tet_p1(q_pts); [~, grad_ref] = lagrange_tet_p1(q_pts); 
                    [I, J, V] = assemble_coupling_kernel_mex(packedRow.P, double(packedRow.T), double(packedRow.CellDofs), double(packedRow.Signs), double(packedCol.CellDofs), Sigma_vec, q_w, ned_ref, grad_ref);
                    valid = (I > 0) & (J > 0);
                    C = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                else
                    [I, J, V] = assemble_coupling_kernel(packedRow.P, packedRow.T, packedRow.CellDofs, packedRow.Signs, packedCol.CellDofs, Sigma_vec, obj.Config);
                    valid = (I > 0) & (J > 0);
                    C = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                end
            else
                error('Unsupported space combination');
            end
        end
        
        % 6. assembleScalarLaplacian
        function K = assembleScalarLaplacian(obj, space, coeffMap)
            packedData = obj.preparePackedData(space);
            packedData.Coeff = obj.expandMaterialMap(coeffMap, obj.Mesh.RegionTags);
            if exist('assemble_scalar_laplacian_kernel_mex', 'file') == 3
                order = 1; [q_pts, q_w] = get_quadrature_data('tet', order); [~, grad_ref] = lagrange_tet_p1(q_pts); 
                [I, J, V] = assemble_scalar_laplacian_kernel_mex(packedData.P, double(packedData.T), double(packedData.CellDofs), packedData.Coeff, q_w, grad_ref);
            else
                [I, J, V] = assemble_scalar_laplacian_kernel(packedData, obj.Config);
            end
            valid = (I > 0) & (J > 0);
            K = sparse(I(valid), J(valid), V(valid), obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
        end
        
        % 7. assembleWinding
        function C = assembleWinding(obj, space, windingObj)
            packedData = obj.preparePackedData(space); packedData.RegionTags = obj.Mesh.RegionTags;
            if exist('assemble_winding_kernel_mex', 'file') == 3
                order = obj.Config.DefaultQuadratureOrder; [q_pts, q_w] = get_quadrature_data('tet', order); [val_ref, ~] = nedelec_tet_p1(q_pts);
                if isempty(windingObj.DirectionField), dirField = []; else, dirField = windingObj.DirectionField; end
                scale = windingObj.Turns / windingObj.CrossSectionArea;
                [I, V] = assemble_winding_kernel_mex(packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), double(packedData.RegionTags), double(windingObj.RegionID), dirField, windingObj.Direction, scale, q_w, val_ref);
                valid = (I > 0);
                C = sparse(I(valid), 1, V(valid), obj.DofHandler.NumGlobalDofs, 1);
            else
                C = assemble_winding_kernel(packedData, windingObj, obj.Config);
            end
        end
        
        % 8. Jacobian
        function [J_mat, R_mat] = assembleJacobian(obj, space, solutionA, matLibData, calcJ)
            if nargin < 5, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            
            if strcmpi(space.Type, 'Nedelec')
                if exist('assemble_jacobian_kernel_mex', 'file') == 3
                    order = max(obj.Config.DefaultQuadratureOrder, 2);
                    [q_pts, q_w] = get_quadrature_data('tet', order);
                    [~, curl_ref] = nedelec_tet_p1(q_pts);
                    
                    [m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt] = obj.packMaterialData(matLibData);
                    
                    [I, J, V, R_idx, R_val] = assemble_jacobian_kernel_mex(...
                        packedData.P, double(packedData.T), double(packedData.CellDofs), double(packedData.Signs), ...
                        double(packedData.RegionTags), solutionA, q_w, curl_ref, ...
                        m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt, double(calcJ));
                    
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
        
        % 9. HBFEM
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
                    
                    % 调用 MEX (22 输入)
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
                        m_lin, m_isNon, m_maxB, m_brk, m_coe, m_start, m_cnt, ...
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
        
        % Helpers
        function packedData = preparePackedData(obj, space)
            packedData = obj.DofHandler.packForKernel(space);
            packedData.P = obj.Mesh.P;
            packedData.T = obj.Mesh.T;
        end
        
        function vec = expandMaterialMap(~, mapOrArray, tags)
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
                        vec = mapOrArray(tags); 
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
        function mat = expandSourceMap(~, mapOrArray, tags)
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
        
        function [linNu, isNon, maxB, breaks, coefs, starts, counts, coe_starts] = packMaterialData(~, matLibData)
            % [V7] Improved Packing + Warning
            if isa(matLibData, 'containers.Map'), loop_keys = cell2mat(matLibData.keys); maxTag = max(loop_keys);
            else, loop_keys = 1:length(matLibData); maxTag = length(matLibData); end
            
            linNu=zeros(maxTag+1,1); isNon=zeros(maxTag+1,1); maxB=zeros(maxTag+1,1); 
            starts=zeros(maxTag+1,1); counts=zeros(maxTag+1,1);
            coe_starts=zeros(maxTag+1,1); 
            
            total_breaks_len=0; total_coefs_len=0;
            for tag=loop_keys
                if tag<=0, continue; end; mat=matLibData(tag);
                if strcmp(mat.Type,'Nonlinear'), pp=mat.SplineNu; total_breaks_len=total_breaks_len+length(pp.breaks); total_coefs_len=total_coefs_len+numel(pp.coefs); end
            end
            breaks=zeros(total_breaks_len,1); coefs=zeros(total_coefs_len,1);
            
            cpp_idx_ptr=0; 
            coe_idx_ptr=0; 
            brk_ptr=1; coe_ptr=1;
            
            for tag=loop_keys
                if tag<=0, continue; end; mat=matLibData(tag);
                if strcmp(mat.Type,'Linear')
                    linNu(tag)=mat.Nu_Linear; isNon(tag)=0;
                else
                    linNu(tag)=mat.nu0; isNon(tag)=1; 
                    if isfield(mat,'MaxBSq'), maxB(tag)=mat.MaxBSq; else, maxB(tag)=1e10; end
                    pp=mat.SplineNu; b_vec=pp.breaks(:); c_mat=pp.coefs; n_b=length(b_vec); n_c=numel(c_mat);
                    
                    starts(tag)=cpp_idx_ptr; 
                    counts(tag)=length(b_vec)-1;
                    coe_starts(tag)=coe_idx_ptr; 
                    
                    breaks(brk_ptr:brk_ptr+n_b-1)=b_vec;
                    c_linear=c_mat'; 
                    coefs(coe_ptr:coe_ptr+n_c-1)=c_linear(:);
                    
                    % [DIAGNOSTIC] Check for zeros in coefficients (suspicious)
                    if all(c_linear(:) == 0)
                        warning('Assembler:ZeroCoefs', 'Material Tag %d has ALL ZERO coefficients! This will cause singularity.', tag);
                    end
                    
                    brk_ptr=brk_ptr+n_b; coe_ptr=coe_ptr+n_c; 
                    cpp_idx_ptr=cpp_idx_ptr+n_b;
                    coe_idx_ptr=coe_idx_ptr+n_c; 
                end
            end
        end
    end
end