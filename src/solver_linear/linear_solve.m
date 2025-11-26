function x = linear_solve(A, b, Model)
% LINEAR_SOLVE 线性方程组求解器接口 (集成 MUMPS)
% 升级: 显式开启 MUMPS 自动缩放，解决物理参数量级差异导致的病态问题

    method = 'Backslash';
    if isfield(Model, 'Solver') && isfield(Model.Solver, 'Linear')
        method = Model.Solver.Linear.Interface;
    end
    
    switch method
        case 'MUMPS'
            if exist('initmumps', 'file') ~= 2
                warning('未找到 initmumps。回退到 MATLAB 内置求解器。');
                x = A \ b; return;
            end
            
            is_complex_system = ~isreal(A) || ~isreal(b);
            id = initmumps;
            id.JOB = -1;
            id.SYM = 0; 
            if isfield(Model.Solver.Linear, 'Symmetric') && Model.Solver.Linear.Symmetric
                id.SYM = 2; 
            end
            
            if is_complex_system, id = zmumps(id); else, id = dmumps(id); end
            
            % ----------------------------------------------------
            % [核心升级] 开启自动缩放与匹配
            % ----------------------------------------------------
            id.ICNTL(1:4) = 0; 
            id.ICNTL(7) = 5;   % METIS
            id.ICNTL(16) = 0;  % OMP Threads
            
            % ICNTL(8): 缩放策略 (Scaling Strategy)
            % 77 = Automatic (通常结合了行列平衡和匹配算法)
            % 这对解决 K(1e6) 和 C(1) 混合的矩阵至关重要
            id.ICNTL(8) = 77; 
            
            % 允许用户覆盖
            if isfield(Model.Solver.Linear, 'MumpsICNTL')
                user_icntl = Model.Solver.Linear.MumpsICNTL;
                idx_list = fieldnames(user_icntl);
                for k = 1:length(idx_list)
                    idx_str = idx_list{k};
                    val = user_icntl.(idx_str);
                    idx = str2double(erase(idx_str, 'i')); 
                    if ~isnan(idx), id.ICNTL(idx) = val; end
                end
            end
            
            % 执行求解
            id.JOB = 6;
            id.RHS = b;
            
            if is_complex_system
                if isreal(A), A = complex(A); end
                id = zmumps(id, A);
            else
                id = dmumps(id, A);
            end
            
            if isempty(id.SOL)
                id.JOB = -2;
                if is_complex_system, zmumps(id); else, dmumps(id); end
                error('MUMPS 求解失败 (INFO(1)=%d)。', id.INFO(1));
            end
            x = id.SOL;
            
            id.JOB = -2;
            if is_complex_system, id = zmumps(id); else, id = dmumps(id); end
            
        case 'Backslash'
            x = A \ b;
            
        case 'GMRES'
            % GMRES 对缩放非常敏感，建议先做 Jacobi 预条件
            [L, U] = ilu(A, struct('type','ilutp','droptol',1e-3));
            [x, flag] = gmres(A, b, [], 1e-6, 1000, L, U);
            if flag ~= 0, warning('GMRES flag=%d', flag); end
            
        otherwise
            error('未知求解器: %s', method);
    end
end