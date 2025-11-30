classdef LinearSolver < handle
    % LINEARSOLVER 线性方程组求解器封装 (v3.6 - Symmetry Support)
    %
    % 更新:
    %   1. [Symmetry] 增加 MumpsSymmetry 属性，支持设置矩阵对称性 (SYM)。
    %      0=非对称(默认), 1=正定对称, 2=一般对称。
    
    properties
        Method = 'Auto' % 'Auto', 'MUMPS', 'Backslash', 'GMRES'
        Tolerance = 1e-6
        MaxIterations = 1000
        MumpsICNTL = struct('i14', 40) 
        ReuseAnalysis = false
        MumpsSymmetry = 0; 
        MumpsID 
    end
    
    methods
        function obj = LinearSolver(method)
            if nargin > 0, obj.Method = method; end
        end
        
        function x = solve(obj, A, b)
            if isempty(A) || size(A, 1) == 0
                warning('LinearSolver:EmptySystem', 'System matrix is empty (0x0).');
                x = []; return;
            end
            
            use_method = obj.Method;
            if strcmpi(use_method, 'Auto')
                hasMumps = (exist('dmumps', 'file') >= 2) || (exist('zmumps', 'file') >= 2);
                if size(A, 1) > 100 && hasMumps
                    use_method = 'MUMPS';
                else
                    use_method = 'Backslash';
                end
            end
            
            switch upper(use_method)
                case 'BACKSLASH'
                    x = A \ b;
                    
                case 'MUMPS'
                    try
                        is_complex = ~isreal(A) || ~isreal(b);
                        
                        % 初始化 (如果需要)
                        if isempty(obj.MumpsID) || ~obj.ReuseAnalysis
                            % 1. 获取默认结构体
                            obj.MumpsID = initmumps;
                            
                            % 2. 初始化实例 (JOB = -1)
                            if is_complex
                                obj.MumpsID = zmumps(obj.MumpsID);
                            else
                                obj.MumpsID = dmumps(obj.MumpsID);
                            end
                            
                            % 3. [New] 设置对称性 (必须在 Analysis 之前设置)
                            % 注意: 这里将类的属性传递给 MUMPS 结构体
                            obj.MumpsID.SYM = obj.MumpsSymmetry;
                            
                            % 4. 设置控制参数 (ICNTL)
                            obj.MumpsID.JOB = 1; % 准备 Analysis
                            if isfield(obj.MumpsICNTL, 'i14')
                                if length(obj.MumpsID.ICNTL) >= 14
                                    obj.MumpsID.ICNTL(14) = obj.MumpsICNTL.i14;
                                end
                            end
                            
                            % 5. 执行分析 (Analysis)
                            if is_complex
                                obj.MumpsID = zmumps(obj.MumpsID, A);
                            else
                                obj.MumpsID = dmumps(obj.MumpsID, A);
                            end
                        end
                        
                        % 6. 因子分解与求解 (Factorization & Solve)
                        obj.MumpsID.JOB = 6; 
                        obj.MumpsID.RHS = b;
                        
                        if is_complex
                            obj.MumpsID = zmumps(obj.MumpsID, A);
                        else
                            obj.MumpsID = dmumps(obj.MumpsID, A);
                        end
                        
                        x = obj.MumpsID.SOL;
                        
                        % 7. 清理
                        if ~obj.ReuseAnalysis
                            obj.MumpsID.JOB = -2;
                            try 
                                if is_complex, zmumps(obj.MumpsID); else, dmumps(obj.MumpsID); end 
                            catch
                            end
                            obj.MumpsID = [];
                        end
                        
                    catch ME
                        fprintf('      [LinearSolver] Warning: MUMPS failed (%s). Fallback to Backslash.\n', ME.message);
                        obj.MumpsID = [];
                        x = A \ b;
                    end
                    
                case 'GMRES'
                    [x, ~] = gmres(A, b, [], obj.Tolerance, obj.MaxIterations);
                    
                otherwise
                    error('Unknown solver method: %s', use_method);
            end
        end
        
        function delete(obj)
            if ~isempty(obj.MumpsID)
                obj.MumpsID.JOB = -2;
                try dmumps(obj.MumpsID); catch; try zmumps(obj.MumpsID); catch; end; end
            end
        end
    end
end