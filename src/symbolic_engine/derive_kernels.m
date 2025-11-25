function derive_kernels()
% DERIVE_KERNELS 符号推导引擎 (修复版)
% 修复了 'sigma' 变量与 MATLAB 内置函数命名冲突的问题
% 使用 MATLAB Symbolic Toolbox 推导四面体单元矩阵并生成代码
% 输出目录: src/kernel_generated/

    % 0. 配置环境
    target_dir = fullfile('src', 'kernel_generated');
    if ~exist(target_dir, 'dir'), mkdir(target_dir); end
    
    fprintf('启动符号推导引擎...\n');
    
    % 1. 定义符号变量
    % 修改：使用 sigma_cond 代替 sigma，避免与内置函数冲突
    syms u v w real % 参考坐标 (Xi, Eta, Zeta)
    syms x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 real % 节点物理坐标
    syms nu sigma_cond real % 材料属性: 磁阻率, 电导率
    
    % 2. 定义几何映射 (参考四面体 -> 物理四面体)
    % 参考节点: N1(0,0,0), N2(1,0,0), N3(0,1,0), N4(0,0,1)
    
    % 节点形函数 (Lagrange)
    L1 = 1 - u - v - w;
    L2 = u;
    L3 = v;
    L4 = w;
    L = [L1; L2; L3; L4]; % 4x1 向量
    
    % 节点坐标矩阵 P (3x4)
    P_mat = [x1, x2, x3, x4;
             y1, y2, y3, y4;
             z1, z2, z3, z4];
         
    % 物理坐标 r (3x1) = P * L
    r = P_mat * L; 
    
    % 计算雅可比矩阵 J = d(x,y,z)/d(u,v,w)
    J = jacobian(r, [u, v, w]);
    detJ = det(J);      % 雅可比行列式
    invJ = inv(J);      % 雅可比逆
    J_invT = invJ.';    % 逆转置 J^{-T}
    
    % ---------------------------------------------------------------------
    % 3. 定义棱单元基函数 (Nédélec First Kind)
    % ---------------------------------------------------------------------
    % 定义6条棱的局部节点对: 
    % 1:(1,2), 2:(1,3), 3:(1,4), 4:(2,3), 5:(2,4), 6:(3,4)
    
    % 参考空间的梯度 (3x1)
    grad_L = [gradient(L1, [u,v,w]), gradient(L2, [u,v,w]), ...
              gradient(L3, [u,v,w]), gradient(L4, [u,v,w])];
    
    % 构造参考基函数 (Whitney 1-forms)
    edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
    N_edge_ref = sym(zeros(3, 6)); 
    
    for i = 1:6
        n_a = edge_defs(i, 1);
        n_b = edge_defs(i, 2);
        % N_e = L_a * grad(L_b) - L_b * grad(L_a)
        N_edge_ref(:, i) = L(n_a) * grad_L(:, n_b) - L(n_b) * grad_L(:, n_a);
    end
    
    % Piola 变换: 映射到物理空间 N_phy = J^{-T} * N_ref
    % 注意：H(curl) 变换不需要除以 detJ
    N_edge_phy = J_invT * N_edge_ref;
    
    % 计算物理空间的旋度: Curl(N)
    % Curl变换: Curl_phy(N) = (1/detJ) * J * Curl_ref(N_ref)
    Curl_N_ref = sym(zeros(3, 6));
    for i = 1:6
        Curl_N_ref(:, i) = curl(N_edge_ref(:, i), [u, v, w]);
    end
    Curl_N_phy = (1/detJ) * J * Curl_N_ref;
    
    % ---------------------------------------------------------------------
    % 4. 构造单元矩阵 (积分)
    % ---------------------------------------------------------------------
    fprintf('正在构造刚度矩阵 (Curl-Curl)...\n');
    
    % [K1] 磁场刚度矩阵: Integral [ (Curl N)' * nu * (Curl N) ] dOmega
    % 变量替换 dOmega = abs(detJ) du dv dw
    integrand_K = Curl_N_phy.' * (nu * eye(3)) * Curl_N_phy * abs(detJ);
    
    % 由于 Curl_N_phy 是常数(一阶棱单元特性)，不需要对u,v,w积分
    % 参考四面体体积 = 1/6
    Ke_curl_curl = simplify(integrand_K * (1/6)); 
    
    fprintf('正在构造质量矩阵 (Mass)...\n');
    
    % [M1] 磁场质量矩阵 (涡流项): Integral [ N' * sigma * N ] dOmega
    % 这里必须使用 sigma_cond
    integrand_M = N_edge_phy.' * (sigma_cond * eye(3)) * N_edge_phy * abs(detJ);
    % 需要在参考四面体上进行真正的积分
    Me_edge_edge = integrate_over_ref_tet(integrand_M, u, v, w);
    
    % [K2] 标量位拉普拉斯: Integral [ grad(L)' * sigma * grad(L) ] dOmega
    % 物理梯度: grad_phy = J^{-T} * grad_ref
    grad_L_phy = J_invT * grad_L; 
    integrand_K_grad = grad_L_phy.' * (sigma_cond * eye(3)) * grad_L_phy * abs(detJ);
    Ke_grad_grad = simplify(integrand_K_grad * (1/6));
    
    % [C] 耦合矩阵: Integral [ N' * sigma * grad(L) ] dOmega
    integrand_C = N_edge_phy.' * (sigma_cond * eye(3)) * grad_L_phy * abs(detJ);
    Ce_edge_grad = integrate_over_ref_tet(integrand_C, u, v, w);
    
    % ---------------------------------------------------------------------
    % 5. 代码生成 (Export)
    % ---------------------------------------------------------------------
    fprintf('正在生成优化代码...\n');
    
    % 注意：虽然符号变量叫 sigma_cond，但生成的函数参数名可以更简洁，但为了避免
    % 调用时混淆，我们在Vars里传入符号变量，生成的函数签名会自动匹配。
    
    % 5.1 生成 Ke_curl_curl
    matlabFunction(Ke_curl_curl, 'File', fullfile(target_dir, 'Ke_curl_curl'), ...
        'Vars', {P_mat, nu}, 'Optimize', true);
        
    % 5.2 生成 Me_edge_edge
    matlabFunction(Me_edge_edge, 'File', fullfile(target_dir, 'Me_edge_edge'), ...
        'Vars', {P_mat, sigma_cond}, 'Optimize', true);
    
    % 5.3 生成 Ke_grad_grad
    matlabFunction(Ke_grad_grad, 'File', fullfile(target_dir, 'Ke_grad_grad'), ...
        'Vars', {P_mat, sigma_cond}, 'Optimize', true);
        
    % 5.4 生成 Ce_edge_grad
    matlabFunction(Ce_edge_grad, 'File', fullfile(target_dir, 'Ce_edge_grad'), ...
        'Vars', {P_mat, sigma_cond}, 'Optimize', true);

    fprintf('代码生成完成! 请检查 %s\n', target_dir);
end

function I = integrate_over_ref_tet(expr, u, v, w)
    % 辅助函数：在参考四面体上积分
    % 积分顺序：先w (0 -> 1-u-v), 再v (0 -> 1-u), 最后u (0 -> 1)
    
    % 提示：对于大型矩阵，int可能会慢。使用 'IgnoreAnalyticConstraints' 有时能加速
    I_w = int(expr, w, 0, 1 - u - v);
    I_v = int(I_w, v, 0, 1 - u);
    I = int(I_v, u, 0, 1);
end