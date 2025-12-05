classdef LinearSolver < handle
    % LINEARSOLVER 线性方程组求解器封装 (v3.7 - Robust MUMPS & Comments)
    % 
    % 功能描述:
    %   该类为有限元求解提供统一的线性方程组求解接口 (Ax=b)。
    %   它自动管理 MATLAB 内置求解器 (\) 与高性能稀疏直接求解器 (MUMPS) 之间的切换。
    %
    % 核心特性:
    %   1. 自动回退: 若未安装 MUMPS，自动降级使用 MATLAB Backslash。
    %   2. 复数支持: 自动检测矩阵虚部，切换实数(dmumps)/复数(zmumps)内核。
    %   3. 健壮的错误处理: 修正了 MEX 调用异常时的内存清理逻辑，防止错误掩盖。
    %
    % 属性配置:
    %   MumpsSymmetry: 0=非对称(默认), 1=正定对称(SPD), 2=一般对称(Symmetric Indefinite)
    %   MumpsICNTL:    结构体，用于覆盖默认的 MUMPS 控制参数 (如内存分配率)。
    
    properties
        Method = 'Auto'     % 求解方法: 'Auto', 'MUMPS', 'Backslash', 'Iterative'
        Tolerance = 1e-6    % 迭代求解器的收敛容差
        MaxIter = 1000      % 迭代求解器的最大步数
        
        % --- MUMPS 配置 ---
        % SYM=0: Unsymmetric (处理 HBFEM、涡流场等非对称矩阵)
        % SYM=1: SPD (处理静磁场、拉普拉斯等正定对称矩阵，速度最快)
        % SYM=2: General Symmetric (处理 A-V 形式等不定对称矩阵)
        MumpsSymmetry = 0; 
        
        % ICNTL 控制参数 (默认 i14=100 表示内存预分配增加 100%)
        MumpsICNTL = struct('i14', 100) 
    end
    
    methods
        function obj = LinearSolver(method)
            % 构造函数
            if nargin > 0
                obj.Method = method;
            else
                obj.Method = 'Auto';
            end
            % 默认设置: ICNTL(7)=5 (Metis 排序)，通常对 3D 网格最优
            obj.MumpsICNTL.i7 = 5; 
        end
        
        function x = solve(obj, A, b)
            % SOLVE 求解 Ax = b
            
            % 1. 输入数据检查 (防御性编程)
            if any(isnan(A), 'all') || any(isinf(A), 'all')
                error('LinearSolver:InputNaN', 'System matrix A contains NaN or Inf.');
            end
            if any(isnan(b), 'all') || any(isinf(b), 'all')
                error('LinearSolver:InputNaN', 'RHS vector b contains NaN or Inf.');
            end
            
            % 2. 确定求解策略
            use_method = obj.Method;
            if strcmpi(use_method, 'Auto')
                % 检查是否有名为 initmumps 的 MEX 文件存在
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
                    % MATLAB 内置稀疏直接求解器 (UMFPACK/CHOLMOD)
                    % 极其稳定，但内存消耗通常高于 MUMPS
                    x = A \ b;
                    
                case 'ITERATIVE'
                    % 迭代求解器 (仅作备用，通常需要预条件子才能收敛)
                    if size(b, 2) > 1
                         warning('Iterative solver called with MRHS. Switching to Backslash.');
                         x = A \ b;
                    else
                        % 根据对称性选择算法
                        if obj.MumpsSymmetry == 1 || obj.MumpsSymmetry == 2
                            % 对称矩阵尝试 MINRES
                            [x, flag] = minres(A, b, obj.Tolerance, obj.MaxIter);
                        else
                            % 非对称矩阵使用 GMRES
                            [x, flag] = gmres(A, b, [], obj.Tolerance, obj.MaxIter);
                        end
                        
                        if flag ~= 0
                            warning('Iterative solver did not converge (Flag: %d). Result may be inaccurate.', flag);
                        end
                    end
                    
                otherwise
                    error('Unknown solver method: %s', obj.Method);
            end
        end
    end
    
    methods (Access = private)
        function x = solveMumps(obj, A, b)
            % SOLVEMUMPS 封装 MUMPS 调用流程
            % 流程: Init -> Analysis -> Factorization -> Solve -> Terminate
            
            % --- 1. 初始化 (JOB = -1) ---
            id = initmumps;
            id.JOB = -1; 
            
            % 设置对称性模式
            if obj.MumpsSymmetry == 1
                id.SYM = 1; % Symmetric Positive Definite
            elseif obj.MumpsSymmetry == 2
                id.SYM = 2; % General Symmetric
            else
                id.SYM = 0; % Unsymmetric (Default)
            end
            
            % 检测是否需要复数算术
            is_complex = ~isreal(A) || ~isreal(b);
            
            % 调用初始化 (实数或复数版本)
            if is_complex
                id = zmumps(id); 
            else
                id = dmumps(id); 
            end
            
            % --- 2. 配置控制参数 (ICNTL) ---
            % 先抑制标准输出 (除非调试需要)
            id.ICNTL(1:4) = 0; 
            
            % 应用用户自定义参数
            if ~isempty(obj.MumpsICNTL)
                fnames = fieldnames(obj.MumpsICNTL);
                for k = 1:length(fnames)
                    name = fnames{k};
                    % 解析 'i14' -> 14
                    idx = str2double(erase(name, 'i')); 
                    if ~isnan(idx) && idx > 0 && idx <= 40
                        id.ICNTL(idx) = obj.MumpsICNTL.(name); 
                    end
                end
            end
            
            % --- 3. 一步求解 (JOB = 6) ---
            % JOB=6 包含: Analysis(1) + Factorization(2) + Solve(3)
            id.JOB = 6; 
            
            % 准备右端项 RHS
            % MUMPS 要求 RHS 是满矩阵 (Full Matrix)
            if issparse(b)
                id.RHS = full(b); 
            else
                id.RHS = b;
            end
            
            try
                % 执行核心计算
                % 注意: zmumps 要求输入矩阵必须是 complex 类型，
                % 如果 A 是 sparse real 但 is_complex 为真 (例如 RHS 是复数)，
                % 必须显式转换 A，否则 MEX 可能崩溃。
                if is_complex
                    if isreal(A), A = complex(A); end 
                    id = zmumps(id, A); % [Critical] 必须接收返回的 id
                else
                    id = dmumps(id, A); % [Critical] 必须接收返回的 id
                end
                
            catch ME
                % --- 异常处理 ---
                fprintf('[MUMPS Error] MEX call failed: %s\n', ME.message);
                
                % 尝试释放内存 (JOB = -2)
                id.JOB = -2;
                try
                    if is_complex
                        id = zmumps(id); % [Critical Fix] 接收 id
                    else
                        id = dmumps(id); % [Critical Fix] 接收 id
                    end
                catch
                    % 如果清理也失败，通常意味着内存已损坏，无法做更多操作
                end
                
                % 抛出原始错误
                rethrow(ME);
            end
            
            % --- 4. 结果检查 ---
            % INFOG(1) < 0 表示出错
            if id.INFOG(1) < 0
                err1 = id.INFOG(1); 
                err2 = id.INFOG(2);
                
                % 尝试清理
                id.JOB = -2;
                if is_complex, id = zmumps(id); else, id = dmumps(id); end
                
                % 常见错误代码解释
                errMsg = sprintf('MUMPS Solver Failed: INFOG(1)=%d, INFOG(2)=%d', err1, err2);
                if err1 == -9
                    errMsg = [errMsg, ' (Error -9: Main memory allocation failed. Try increasing ICNTL(14))'];
                elseif err1 == -10
                    errMsg = [errMsg, ' (Error -10: Numerically singular matrix)'];
                end
                error('LinearSolver:MUMPSFailed', errMsg);
            end
            
            % 提取解向量
            x = id.SOL;
            
            % --- 5. 正常退出清理 (JOB = -2) ---
            id.JOB = -2; 
            if is_complex
                id = zmumps(id); % [Critical Fix] 接收 id
            else
                id = dmumps(id); % [Critical Fix] 接收 id
            end
        end
    end
end