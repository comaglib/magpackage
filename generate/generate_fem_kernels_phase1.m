function generate_fem_kernels_phase1()
% GENERATE_FEM_KERNELS_PHASE1
% 符号计算生成矢量化FEM内核：第一阶段（P1 + Nedelec）
%
% 功能：
%   1. 定义四面体几何与映射
%   2. 推导 P1 刚度和质量矩阵
%   3. 推导 Nedelec 旋度和质量矩阵（含符号修正）
%   4. 导出优化后的.m 文件

    %% 1. 环境初始化
    clear; clc;
    fprintf('启动第一阶段符号生成引擎...\n');
    
    % 定义实数符号变量（核心优化：避免复数运算分支）
    syms x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 real
    syms xi eta zeta real
    
    % 定义 Nédélec 棱边方向符号变量 (+1/-1)
    s = sym('s', [6 1], 'real'); 
    
    % 构造顶点矩阵
    V = [x1 y1 z1; x2 y2 z2; x3 y3 z3; x4 y4 z4];
    
    % 输入变量列表（用于 matlabFunction）
    % 将 12 个坐标分量独立列出，以支持 SoA 矢量输入
    vars_geo = {x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4};
    
    %% 2. 几何映射与雅可比计算
    fprintf('  - 计算几何映射与雅可比矩阵...\n');
    
    % 重心坐标 lambda (4x1)
    lambda = [1 - xi - eta - zeta; xi; eta; zeta];
          
    % 梯度（参考空间） d_lambda / d_xi
    Grad_ref = jacobian(lambda, [xi, eta, zeta]); 
    
    % 物理坐标映射 X = N * V
    X_phys = lambda.' * V;
    
    % 雅可比矩阵 J
    J = jacobian(X_phys, [xi, eta, zeta]);
    detJ = det(J);
    
    % 梯度变换矩阵：inv(J).'
    % 使用 adjoint/det 形式可能比直接 inv 更利于符号简化
    % 但 modern symbolic toolbox 对 inv 处理已足够好
    G = inv(J); 
    
    %% 3. P1 单元代码生成
    fprintf('  - 推导 P1 单元矩阵...\n');
    
    % 物理梯度：Grad_ref * G
    Grad_phys_P1 = Grad_ref * G;
    
    % 刚度矩阵 integrand
    % 标量积：sum((dNi/dx)*(dNj/dx))
    K_P1_sym = (Grad_phys_P1 * Grad_phys_P1.') * abs(detJ);
    
    % 积分：常数 integrand * Volume (1/6)
    % 注意：这里直接除以6，因为是在参考单元(Vol=1/6)积分常数
    K_P1 = K_P1_sym / 6;
    
    % 导出 P1 刚度矩阵
    matlabFunction(K_P1, 'File', 'P1_Stiffness_Tet',...
        'Vars', vars_geo, 'Optimize', true);
        
    %% 4. Nédélec 单元代码生成
    fprintf('  - 推导 Nédélec 单元矩阵...\n');
    
    % 定义棱边局部索引
    edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    n_edges = 6;
    
    W_ref = sym(zeros(n_edges, 3));
    
    % 构造参考基函数 (Whitney Forms)
    for k = 1:n_edges
        i = edges(k, 1);
        j = edges(k, 2);
        W_ref(k, :) = lambda(i)*Grad_ref(j,:) - lambda(j)*Grad_ref(i,:);
    end
    
    % 物理基函数映射 (协变 + 符号修正)
    % W_phys = s_k * (J^-T * W_ref)
    W_phys = sym(zeros(n_edges, 3));
    invJ_T = G.';
    
    W_phys_raw = (invJ_T * W_ref.').'; % 批量映射
    
    for k = 1:n_edges
        W_phys(k, :) = W_phys_raw(k, :) * s(k);
    end
    
    %% 4a. 旋度矩阵 (Curl-Curl)
    % 旋度变换：(1/detJ) * J * Curl_ref
    
    Curl_ref = sym(zeros(n_edges, 3));
    % 手动计算参考旋度（避免调用 curl 函数的坐标系歧义）
    % curl(u) = [du3/deta - du2/dzeta,...]
    for k = 1:n_edges
        u = W_ref(k, :);
        dudxi = jacobian(u, xi);
        dudeta = jacobian(u, eta);
        dudzeta = jacobian(u, zeta);
        
        c1 = dudeta(3) - dudzeta(2);
        c2 = dudzeta(1) - dudxi(3);
        c3 = dudxi(2) - dudeta(1);
        Curl_ref(k, :) = [c1 c2 c3];
    end
    
    % 映射到物理空间
    Curl_phys = (J * Curl_ref.').' / detJ;
    
    % 应用符号修正
    for k = 1:n_edges
        Curl_phys(k, :) = Curl_phys(k, :) * s(k);
    end
    
    % 积分：(Curl. Curl) * abs(detJ) / 6
    S_Ned = (Curl_phys * Curl_phys.') * abs(detJ) / 6;
    
    % 导出：输入需增加符号向量 s
    vars_ned = [vars_geo, {s}];
    matlabFunction(S_Ned, 'File', 'Nedelec_Stiffness_Tet',...
        'Vars', vars_ned, 'Optimize', true);
        
    %% 4b. 质量矩阵 (精确积分)
    fprintf('    > 正在计算 Nédélec 质量矩阵积分 (耗时操作)...\n');
    
    M_Ned = sym(zeros(n_edges));
    
    % 为了加速，先化简被积函数结构
    % M_ij = s_i s_j * integral( (J^-T W_ref_i). (J^-T W_ref_j) * abs(detJ) )
    
    % 这里我们演示串行积分，实际部署建议使用 parfor
    for i = 1:n_edges
        for j = i:n_edges
            % 构造被积函数 (关于 xi, eta, zeta 的二次型)
            integrand = dot(W_phys(i,:), W_phys(j,:)) * abs(detJ);
            
            % 三重积分：0->1, 0->1-xi, 0->1-xi-eta
            val = int(integrand, zeta, 0, 1 - xi - eta);
            val = int(val, eta, 0, 1 - xi);
            val = int(val, xi, 0, 1);
            
            M_Ned(i,j) = val;
            M_Ned(j,i) = val; % 对称性
        end
    end
    
    matlabFunction(M_Ned, 'File', 'Nedelec_Mass_Tet',...
        'Vars', vars_ned, 'Optimize', true);
        
    fprintf('第一阶段代码生成完成。\n');
end