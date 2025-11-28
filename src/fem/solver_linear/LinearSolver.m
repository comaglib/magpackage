classdef LinearSolver < handle
    % LINEARSOLVER 线性方程组求解器封装 (v3.4 - Stability Guard)
    % 
    % 修复:
    %   1. 增加了 NaN/Inf 检查，防止传入坏数据导致 MUMPS 崩溃。
    %   2. 修正了 MUMPS 调用语法。
    
    properties
        Method = 'Auto'     % 'Auto', 'MUMPS', 'Backslash', 'Iterative'
        Tolerance = 1e-6    % 迭代求解器容差
        MaxIter = 1000      % 最大迭代步数
        Symmetric = true    % 矩阵是否对称
        
        % MUMPS 高级控制参数
        MumpsICNTL = struct() 
    end
    
    methods
        function obj = LinearSolver(method)
            if nargin > 0
                obj.Method = method;
            else
                obj.Method = 'Auto';
            end
            obj.MumpsICNTL.i7 = 5; 
        end
        
        function x = solve(obj, A, b)
            % [Safety Check] 防止 NaN/Inf 导致解释器崩溃
            if any(isnan(A), 'all') || any(isinf(A), 'all')
                error('LinearSolver:InputNaN', 'System matrix A contains NaN or Inf. Solver aborted to prevent crash.');
            end
            if any(isnan(b), 'all') || any(isinf(b), 'all')
                error('LinearSolver:InputNaN', 'RHS vector b contains NaN or Inf. Solver aborted to prevent crash.');
            end

            use_method = obj.Method;
            
            % 自动选择策略
            if strcmpi(use_method, 'Auto')
                if exist('initmumps', 'file') == 2 || exist('initmumps', 'file') == 3
                    use_method = 'MUMPS';
                else
                    use_method = 'Backslash';
                end
            end
            
            switch upper(use_method)
                case 'MUMPS'
                    x = obj.solveMumps(A, b);
                    
                case {'BACKSLASH', 'DIRECT'}
                    x = A \ b;
                    
                case 'ITERATIVE'
                    if obj.Symmetric
                        [x, flag] = minres(A, b, obj.Tolerance, obj.MaxIter);
                    else
                        [x, flag] = gmres(A, b, [], obj.Tolerance, obj.MaxIter);
                    end
                    if flag ~= 0
                        warning('Iterative solver did not converge (Flag: %d)', flag);
                    end
                    
                otherwise
                    error('Unknown solver method: %s', obj.Method);
            end
        end
    end
    
    methods (Access = private)
        function x = solveMumps(obj, A, b)
            % MUMPS 求解核心逻辑
            
            id = initmumps;
            id.JOB = -1; 
            id.SYM = 0;  
            if obj.Symmetric
                id.SYM = 2; 
            end
            
            is_complex = ~isreal(A) || ~isreal(b);
            if is_complex
                id = zmumps(id);
            else
                id = dmumps(id);
            end
            
            id.ICNTL(1:4) = 0; % 关闭输出
            
            if ~isempty(obj.MumpsICNTL)
                fnames = fieldnames(obj.MumpsICNTL);
                for k = 1:length(fnames)
                    name = fnames{k};
                    idx = str2double(erase(name, 'i')); 
                    if ~isnan(idx)
                        id.ICNTL(idx) = obj.MumpsICNTL.(name);
                    end
                end
            end
            
            id.JOB = 6; 
            id.RHS = b;
            
            try
                if is_complex
                    if isreal(A), A = complex(A); end 
                    id = zmumps(id, A);
                else
                    id = dmumps(id, A);
                end
            catch ME
                % 捕获可能的 MEX 错误
                fprintf('[MUMPS Error] MEX call failed: %s\n', ME.message);
                id.JOB = -2;
                if is_complex, zmumps(id); else, dmumps(id); end
                rethrow(ME);
            end
            
            % 错误检查
            if id.INFOG(1) < 0
                err1 = id.INFOG(1); err2 = id.INFOG(2);
                id.JOB = -2;
                if is_complex, zmumps(id); else, dmumps(id); end
                error('MUMPS Solver Failed: INFOG(1)=%d, INFOG(2)=%d', err1, err2);
            end
            
            x = id.SOL;
            
            % 清理
            id.JOB = -2; 
            if is_complex
                id = zmumps(id); 
            else
                id = dmumps(id); 
            end
        end
    end
end