classdef LinearSolver < handle
    % LINEARSOLVER 线性方程组求解器封装 (v3.6 - Robust Exception Handling)
    % 
    % 修复:
    %   1. [Critical] 修正 catch 块中的 MUMPS 清理调用语法 (id = dmumps(id))。
    %      这是导致报错 "输出参数数目不足" 的直接原因。
    %   2. 保持了 MRHS 支持和 NaN 检查。
    
    properties
        Method = 'Auto'     
        Tolerance = 1e-6    
        MaxIter = 1000      
        MumpsSymmetry= 0; % 0=Unsymmetric, 1=SPD, 2=General Symmetric
        MumpsICNTL = struct('i14', 100) % 默认内存增加百分比
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
            if any(isnan(A), 'all') || any(isinf(A), 'all')
                error('LinearSolver:InputNaN', 'System matrix A contains NaN or Inf.');
            end
            if any(isnan(b), 'all') || any(isinf(b), 'all')
                error('LinearSolver:InputNaN', 'RHS vector b contains NaN or Inf.');
            end

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
                    x = obj.solveMumps(A, b);
                    
                case {'BACKSLASH', 'DIRECT'}
                    x = A \ b;
                    
                case 'ITERATIVE'
                    if size(b, 2) > 1
                         warning('Iterative solver called with MRHS. Switching to Backslash.');
                         x = A \ b;
                    else
                        if obj.MumpsSymmetry
                            [x, flag] = minres(A, b, obj.Tolerance, obj.MaxIter);
                        else
                            [x, flag] = gmres(A, b, [], obj.Tolerance, obj.MaxIter);
                        end
                        if flag ~= 0
                            warning('Iterative solver did not converge (Flag: %d)', flag);
                        end
                    end
                    
                otherwise
                    error('Unknown solver method: %s', obj.Method);
            end
        end
    end
    
    methods (Access = private)
        function x = solveMumps(obj, A, b)
            id = initmumps;
            id.JOB = -1; 
            id.SYM = 0;  
            if obj.MumpsSymmetry, id.SYM = 2; end
            
            is_complex = ~isreal(A) || ~isreal(b);
            if is_complex, id = zmumps(id); else, id = dmumps(id); end
            
            id.ICNTL(1:4) = 0; 
            if ~isempty(obj.MumpsICNTL)
                fnames = fieldnames(obj.MumpsICNTL);
                for k = 1:length(fnames)
                    name = fnames{k};
                    idx = str2double(erase(name, 'i')); 
                    if ~isnan(idx), id.ICNTL(idx) = obj.MumpsICNTL.(name); end
                end
            end
            
            id.JOB = 6; 
            
            % MUMPS 要求 RHS 为满矩阵
            if issparse(b)
                id.RHS = full(b); 
            else
                id.RHS = b;
            end
            
            try
                if is_complex
                    if isreal(A), A = complex(A); end 
                    id = zmumps(id, A); % 必须接收 id
                else
                    id = dmumps(id, A); % 必须接收 id
                end
            catch ME
                fprintf('[MUMPS Error] MEX call failed: %s\n', ME.message);
                id.JOB = -2;
                % [Fix] 这里的调用之前缺少接收变量，导致掩盖了真实错误
                if is_complex
                    id = zmumps(id); 
                else
                    id = dmumps(id); 
                end
                rethrow(ME);
            end
            
            % 错误检查
            if id.INFOG(1) < 0
                err1 = id.INFOG(1); err2 = id.INFOG(2);
                id.JOB = -2;
                if is_complex, id = zmumps(id); else, id = dmumps(id); end
                error('MUMPS Solver Failed: INFOG(1)=%d, INFOG(2)=%d', err1, err2);
            end
            
            x = id.SOL;
            
            % 清理
            id.JOB = -2; 
            if is_complex
                id = zmumps(id); % [Fix] 必须接收 id
            else
                id = dmumps(id); % [Fix] 必须接收 id
            end
        end
    end
end