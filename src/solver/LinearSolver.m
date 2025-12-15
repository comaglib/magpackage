classdef LinearSolver < handle
    % LINEARSOLVER 线性方程组求解器封装 (v3.9 - Matrix Argument Fix)
    % 
    % 更新日志:
    %   v3.9: 修复了 MUMPS 调用参数不足的问题。
    %         solveFromFactors 现在强制要求传入矩阵 A，以适配 dmumps.m/zmumps.m 的接口要求。
    
    properties
        Method = 'Auto'     % 'Auto', 'MUMPS', 'Backslash', 'Iterative'
        Tolerance = 1e-6    
        MaxIter = 1000      
        MumpsSymmetry = 0;  % 0=Unsymmetric, 1=SPD, 2=General Symmetric
        MumpsICNTL = struct('i14', 100) 
    end
    
    methods
        function obj = LinearSolver(method)
            if nargin > 0, obj.Method = method; else, obj.Method = 'Auto'; end
            obj.MumpsICNTL.i7 = 5; % Metis ordering
        end
        
        function x = solve(obj, A, b)
            % SOLVE 标准一次性求解 Ax=b
            if any(isnan(A), 'all') || any(isinf(A), 'all'), error('Matrix A contains NaN/Inf'); end
            if any(isnan(b), 'all') || any(isinf(b), 'all'), error('RHS b contains NaN/Inf'); end
            
            use_method = obj.Method;
            if strcmpi(use_method, 'Auto')
                if exist('initmumps', 'file') == 2 || exist('initmumps', 'file') == 3
                    use_method = 'MUMPS';
                else
                    use_method = 'Backslash';
                end
            end
            
            switch upper(use_method)
                case 'MUMPS'
                    x = obj.solveMumpsOneShot(A, b);
                case {'BACKSLASH', 'DIRECT'}
                    x = A \ b;
                case 'ITERATIVE'
                    [x, ~] = gmres(A, b, [], obj.Tolerance, obj.MaxIter);
                otherwise
                    error('Unknown method: %s', obj.Method);
            end
        end
        
        % =========================================================
        %  新增接口: 分离式分解与求解 (针对 SDC 等算法优化)
        % =========================================================
        
        function mumpsID = factorize(obj, A)
            % FACTORIZE 执行 MUMPS 分析与分解 (JOB=4)
            % 返回 mumpsID 结构体，包含分解因子。
            
            % 1. 初始化
            id = initmumps;
            id.JOB = -1; 
            if obj.MumpsSymmetry == 1, id.SYM=1; elseif obj.MumpsSymmetry == 2, id.SYM=2; else, id.SYM=0; end
            
            is_complex = ~isreal(A);
            if is_complex, id = zmumps(id); else, id = dmumps(id); end
            
            % 2. 配置
            id.ICNTL(1:4) = 0; % Silence
            if ~isempty(obj.MumpsICNTL)
                fnames = fieldnames(obj.MumpsICNTL);
                for k=1:length(fnames)
                    idx = str2double(erase(fnames{k}, 'i'));
                    if ~isnan(idx), id.ICNTL(idx) = obj.MumpsICNTL.(fnames{k}); end
                end
            end
            
            % 3. 分析 + 分解 (JOB = 4)
            id.JOB = 4;
            try
                if is_complex
                    if isreal(A), A = complex(A); end
                    id = zmumps(id, A);
                else
                    id = dmumps(id, A);
                end
            catch ME
                obj.safeClear(id, is_complex);
                rethrow(ME);
            end
            
            if id.INFOG(1) < 0
                obj.safeClear(id, is_complex);
                error('MUMPS Factorization Failed: INFOG(1)=%d', id.INFOG(1));
            end
            
            mumpsID = id;
        end
        
        function x = solveFromFactors(~, mumpsID, b, A)
            % SOLVEFROMFACTORS 利用已分解的 ID 求解 (JOB=3)
            % [Update] 增加了参数 A，因为 dmumps.m 封装要求必须传入矩阵
            
            is_complex = (isfield(mumpsID, 'TYPE') && mumpsID.TYPE == 2);
            
            mumpsID.JOB = 3; % Solve
            if issparse(b), mumpsID.RHS = full(b); else, mumpsID.RHS = b; end
            
            % 调用求解
            % 注意: 即使是求解步，dmumps.m 也要求传入 A
            if is_complex
                if isreal(A), A = complex(A); end
                mumpsID = zmumps(mumpsID, A);
            else
                mumpsID = dmumps(mumpsID, A);
            end
            
            if mumpsID.INFOG(1) < 0
                error('MUMPS Solve Failed: INFOG(1)=%d', mumpsID.INFOG(1));
            end
            
            x = mumpsID.SOL;
        end
        
        function clearFactors(~, mumpsID)
            % CLEARFACTORS 释放 MUMPS 实例内存 (JOB=-2)
            if isempty(mumpsID), return; end
            
            is_complex = (isfield(mumpsID, 'TYPE') && mumpsID.TYPE == 2);
            mumpsID.JOB = -2;
            
            try
                % 清理时通常不需要传入矩阵，但如果 dmumps 报错，可以尝试传入 []
                if is_complex, mumpsID = zmumps(mumpsID); else, mumpsID = dmumps(mumpsID); end
            catch
                % Ignore cleanup errors
            end
        end
    end
    
    methods (Access = private)
        function x = solveMumpsOneShot(obj, A, b)
            % 一次性求解流程
            id = obj.factorize(A); 
            try
                % [Fix] 传入 A
                x = obj.solveFromFactors(id, b, A);
                obj.clearFactors(id);
            catch ME
                obj.clearFactors(id);
                rethrow(ME);
            end
        end
        
        function safeClear(~, id, is_complex)
            id.JOB = -2;
            try
                if is_complex, id = zmumps(id); else, id = dmumps(id); end
            catch; end
        end
    end
end