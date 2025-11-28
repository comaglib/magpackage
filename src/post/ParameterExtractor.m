classdef ParameterExtractor < handle
    % PARAMETEREXTRACTOR 电路参数提取工具 (v2.0 - MRHS Optimized)
    % 
    % 优化:
    %   1. 使用多右端项 (MRHS) 技术，一次性求解所有绕组激励。
    %   2. 利用 MUMPS 的分解复用能力，极大提升多绕组系统的提取速度。
    
    properties
        Assembler
        LinearSolver
    end
    
    methods
        function obj = ParameterExtractor(assembler)
            obj.Assembler = assembler;
            obj.LinearSolver = LinearSolver('Auto');
            % 增加内存预估，因为 MRHS 需要更多内存存储因子
            obj.LinearSolver.MumpsICNTL.i14 = 40; 
        end
        
        function [L_matrix, R_matrix] = extractMatrices(obj, space, nuMap, windings)
            % EXTRACTMATRICES 提取 L 和 R 矩阵
            
            if iscell(windings)
                N_windings = length(windings);
                W_list = windings;
            else
                N_windings = length(windings);
                W_list = num2cell(windings);
            end
            
            fprintf('==============================================\n');
            fprintf('   Parameter Extraction (N=%d Windings)       \n', N_windings);
            fprintf('   Method: Vectorized MRHS Solve              \n');
            fprintf('==============================================\n');
            
            numDofs = obj.Assembler.DofHandler.NumGlobalDofs;
            
            % 1. 组装所有绕组向量 C
            fprintf('   [Extract] Assembling coupling vectors...\n');
            C_matrix = sparse(numDofs, N_windings);
            R_matrix = zeros(N_windings, N_windings);
            
            for i = 1:N_windings
                C_i = obj.Assembler.assembleWinding(space, W_list{i});
                C_matrix(:, i) = C_i;
                R_matrix(i, i) = W_list{i}.Resistance;
            end
            
            % 2. 组装刚度矩阵
            fprintf('   [Extract] Assembling Stiffness Matrix...\n');
            K = obj.Assembler.assembleStiffness(space, nuMap);
            
            % 3. 边界条件
            is_bnd = BoundaryCondition.findOuterBoundaryDofs(...
                obj.Assembler.Mesh, obj.Assembler.DofHandler, space);
            
            % 4. 构造多右端项 (MRHS)
            % F_all 每一列对应一个绕组的 1A 激励
            % F = C * 1.0
            F_all = C_matrix; 
            
            % 施加边界条件 (支持矩阵 F)
            [K_sys, F_sys] = BoundaryCondition.applyDirichlet(K, F_all, is_bnd);
            
            % 5. 一次性求解所有工况
            fprintf('   [Extract] Solving %d load cases simultaneously...\n', N_windings);
            tic;
            A_all = obj.LinearSolver.solve(K_sys, F_sys);
            t_solve = toc;
            fprintf('     -> Solved in %.4f s\n', t_solve);
            
            % 6. 计算电感矩阵
            % L_ij = Psi_i(j) = C_i' * A_j
            % Matrix form: L = C_matrix' * A_all
            fprintf('   [Extract] Computing Flux Linkages...\n');
            L_matrix = C_matrix' * A_all;
            
            % 确保 L 矩阵对称 (消除数值误差)
            L_matrix = 0.5 * (L_matrix + L_matrix');
            
            fprintf('   [Extract] Done.\n');
            fprintf('==============================================\n');
        end
    end
end