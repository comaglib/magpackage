classdef LinearSolver < handle
    % LINEARSOLVER 线性方程组求解器封装 (v4.1 - MUMPS Fix)
    %
    % 修复:
    %   1. [Critical] 恢复了缺失的 JOB=-1 初始化步骤。
    %      MUMPS 必须先初始化才能执行后续的分析与求解。
    %   2. [Robustness] 优化了 clearMumps，防止对未初始化的实例调用销毁。
    
    properties
        Method = 'Auto'     
        Tolerance = 1e-6    
        MaxIter = 1000      
        Symmetric = true    
        MumpsICNTL = struct()
        
        % 如果为 true，假设矩阵稀疏结构不变，重用符号分析结果
        ReuseAnalysis = false 
    end
    
    properties (Access = private)
        MumpsID = []       % 缓存 MUMPS 实例结构体
        IsComplex = false  % 记录当前实例是否为复数
        Initialized = false
    end
    
    methods
        function obj = LinearSolver(method)
            if nargin > 0
                obj.Method = method;
            else
                obj.Method = 'Auto';
            end
            % 默认设置: 使用 Metis (5) 排序
            obj.MumpsICNTL.i7 = 5; 
            % 增加工作空间内存百分比
            obj.MumpsICNTL.i14 = 50; 
        end
        
        function delete(obj)
            obj.clearMumps();
        end
        
        function clearMumps(obj)
            if ~isempty(obj.MumpsID)
                % 仅当实例已在底层初始化 (INST ~= -9999) 时才调用销毁
                % 否则 dmumps.m 会打印 "Uninitialized instance"
                if isfield(obj.MumpsID, 'INST') && obj.MumpsID.INST ~= -9999
                    try
                        obj.MumpsID.JOB = -2; % Terminate
                        if obj.IsComplex
                            zmumps(obj.MumpsID);
                        else
                            dmumps(obj.MumpsID);
                        end
                    catch
                        % 忽略销毁时的错误
                    end
                end
            end
            obj.MumpsID = [];
            obj.Initialized = false;
        end
        
        function x = solve(obj, A, b)
            % 1. 基础检查
            if any(isnan(A), 'all') || any(isinf(A), 'all')
                error('LinearSolver:InputNaN', 'System matrix A contains NaN or Inf.');
            end
            if any(isnan(b), 'all') || any(isinf(b), 'all')
                error('LinearSolver:InputNaN', 'RHS vector b contains NaN or Inf.');
            end

            % 2. 方法选择
            use_method = obj.Method;
            if strcmpi(use_method, 'Auto')
                if exist('initmumps', 'file') == 2 || exist('initmumps', 'file') == 3
                    use_method = 'MUMPS';
                else
                    use_method = 'Backslash';
                end
            end
            
            % 3. 执行求解
            switch upper(use_method)
                case 'MUMPS'
                    x = obj.solveMumps(A, b);
                    
                case {'BACKSLASH', 'DIRECT'}
                    obj.clearMumps(); 
                    x = A \ b;
                    
                case 'ITERATIVE'
                    obj.clearMumps();
                    if size(b, 2) > 1
                         warning('Iterative solver called with MRHS. Switching to Backslash.');
                         x = A \ b;
                    else
                        if obj.Symmetric
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
            is_input_complex = ~isreal(A) || ~isreal(b);
            
            % 1. 初始化或重置
            if isempty(obj.MumpsID) || (obj.IsComplex ~= is_input_complex)
                obj.clearMumps(); 
                
                % Step 1.1: 创建结构体
                obj.MumpsID = initmumps;
                obj.IsComplex = is_input_complex;
                
                % Step 1.2: 设置初始化参数
                obj.MumpsID.JOB = -1; % Initialization
                obj.MumpsID.SYM = 0;
                if obj.Symmetric, obj.MumpsID.SYM = 2; end
                
                % Step 1.3: 执行初始化调用 (分配底层内存)
                try
                    if obj.IsComplex
                        obj.MumpsID = zmumps(obj.MumpsID);
                    else
                        obj.MumpsID = dmumps(obj.MumpsID);
                    end
                catch ME
                    error('MUMPS Init Failed: %s', ME.message);
                end
                
                obj.Initialized = true;
                
                % Step 1.4: 设置控制参数 (必须在 Init 之后)
                obj.MumpsID.ICNTL(1:4) = 0; % 关闭输出
                
                if ~isempty(obj.MumpsICNTL)
                    fnames = fieldnames(obj.MumpsICNTL);
                    for k = 1:length(fnames)
                        name = fnames{k};
                        idx = str2double(erase(name, 'i')); 
                        if ~isnan(idx), obj.MumpsID.ICNTL(idx) = obj.MumpsICNTL.(name); end
                    end
                end
                
                do_analysis = true; 
            else
                % 实例已存在
                if obj.ReuseAnalysis
                    do_analysis = false; 
                else
                    do_analysis = true;  
                end
            end
            
            % 2. 准备数据
            if issparse(b)
                obj.MumpsID.RHS = full(b); 
            else
                obj.MumpsID.RHS = b;
            end
            
            % 3. 设置 JOB
            if do_analysis
                obj.MumpsID.JOB = 6; % Analysis + Fact + Solve
            else
                obj.MumpsID.JOB = 5; % Fact + Solve
            end
            
            % 4. 调用 MUMPS MEX
            try
                if obj.IsComplex
                    if isreal(A), A = complex(A); end 
                    obj.MumpsID = zmumps(obj.MumpsID, A);
                else
                    obj.MumpsID = dmumps(obj.MumpsID, A);
                end
            catch ME
                fprintf('[MUMPS Error] MEX call failed: %s\n', ME.message);
                obj.clearMumps();
                rethrow(ME);
            end
            
            % 5. 错误检查
            if obj.MumpsID.INFOG(1) < 0
                err1 = obj.MumpsID.INFOG(1); 
                err2 = obj.MumpsID.INFOG(2);
                
                obj.clearMumps();
                
                if err1 == -6
                     error('MUMPS Failed: Matrix is singular in structure (INFOG(1)=-6).');
                elseif err1 == -9
                     error('MUMPS Failed: Not enough memory (INFOG(1)=-9). Increase ICNTL(14).');
                else
                     error('MUMPS Solver Failed: INFOG(1)=%d, INFOG(2)=%d', err1, err2);
                end
            end
            
            % 6. 获取解
            x = obj.MumpsID.SOL;
        end
    end
end