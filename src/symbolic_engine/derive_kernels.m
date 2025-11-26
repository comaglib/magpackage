function derive_kernels()
% DERIVE_KERNELS 符号推导引擎 (v2: 含各向异性内核)
% 输出目录: src/kernel_generated/

    target_dir = fullfile('src', 'kernel_generated');
    if ~exist(target_dir, 'dir'), mkdir(target_dir); end
    
    fprintf('启动符号推导引擎...\n');
    
    % 1. 定义符号变量
    syms u v w real 
    syms x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 real 
    syms nu sigma_cond real 
    
    % 新增: 各向异性张量 (对称 3x3)
    % nu_11, nu_12, ... 我们用向量输入构建矩阵
    syms nu_xx nu_xy nu_xz nu_yy nu_yz nu_zz real
    Nu_tensor = [nu_xx, nu_xy, nu_xz; 
                 nu_xy, nu_yy, nu_yz;
                 nu_xz, nu_yz, nu_zz];
    
    % 2. 几何映射
    L1 = 1 - u - v - w; L2 = u; L3 = v; L4 = w;
    L = [L1; L2; L3; L4];
    P_mat = [x1, x2, x3, x4; y1, y2, y3, y4; z1, z2, z3, z4];
    r = P_mat * L;
    J = jacobian(r, [u, v, w]);
    detJ = det(J);
    invJ_T = inv(J).';
    
    % 3. 棱单元基函数与旋度
    grad_L = [gradient(L1, [u,v,w]), gradient(L2, [u,v,w]), ...
              gradient(L3, [u,v,w]), gradient(L4, [u,v,w])];
          
    edge_defs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
    N_edge_ref = sym(zeros(3, 6)); 
    for i = 1:6
        na = edge_defs(i, 1); nb = edge_defs(i, 2);
        N_edge_ref(:, i) = L(na)*grad_L(:, nb) - L(nb)*grad_L(:, na);
    end
    
    % 物理空间基函数与旋度
    N_edge_phy = invJ_T * N_edge_ref;
    
    Curl_N_ref = sym(zeros(3, 6));
    for i = 1:6
        Curl_N_ref(:, i) = curl(N_edge_ref(:, i), [u, v, w]);
    end
    Curl_N_phy = (1/detJ) * J * Curl_N_ref;
    
    % 4. 构造单元矩阵
    
    % [K1] 标量各向同性刚度 (旧)
    integrand_K = Curl_N_phy.' * (nu * eye(3)) * Curl_N_phy * abs(detJ);
    Ke_curl_curl = simplify(integrand_K * (1/6)); 
    
    % [K_aniso] 各向异性/切线刚度 (新) - 用于 Newton-Raphson Jacobian
    % 积分: (Curl N)' * Nu_tensor * (Curl N) * Vol
    fprintf('正在构造各向异性刚度矩阵 (Jacobian Kernel)...\n');
    integrand_Kaniso = Curl_N_phy.' * Nu_tensor * Curl_N_phy * abs(detJ);
    Ke_aniso = simplify(integrand_Kaniso * (1/6));
    
    % [M1] 质量矩阵
    integrand_M = N_edge_phy.' * (sigma_cond * eye(3)) * N_edge_phy * abs(detJ);
    Me_edge_edge = integrate_over_ref_tet(integrand_M, u, v, w);
    
    % [K2] 节点拉普拉斯
    grad_L_phy = invJ_T * grad_L; 
    integrand_K_grad = grad_L_phy.' * (sigma_cond * eye(3)) * grad_L_phy * abs(detJ);
    Ke_grad_grad = simplify(integrand_K_grad * (1/6));
    
    % [C] 耦合矩阵
    integrand_C = N_edge_phy.' * (sigma_cond * eye(3)) * grad_L_phy * abs(detJ);
    Ce_edge_grad = integrate_over_ref_tet(integrand_C, u, v, w);
    
    % [Aux] 辅助函数: 计算物理空间旋度 B = Curl(A)
    % 给定单元内的 A_local (6x1), 计算 B (3x1, 常数)
    % B = Curl_N_phy * A_local
    % 为了速度，我们生成一个函数直接返回 B，而不是矩阵
    syms a1 a2 a3 a4 a5 a6 real
    A_local = [a1; a2; a3; a4; a5; a6];
    B_vec = Curl_N_phy * A_local;
    
    % 5. 代码生成
    fprintf('正在生成优化代码...\n');
    
    matlabFunction(Ke_curl_curl, 'File', fullfile(target_dir, 'Ke_curl_curl'), ...
        'Vars', {P_mat, nu}, 'Optimize', true);
        
    % 新增: Ke_aniso
    % 输入: P, [nu_xx, nu_xy, nu_xz, nu_yy, nu_yz, nu_zz]
    nu_vec = [nu_xx, nu_xy, nu_xz, nu_yy, nu_yz, nu_zz];
    matlabFunction(Ke_aniso, 'File', fullfile(target_dir, 'Ke_aniso'), ...
        'Vars', {P_mat, nu_vec}, 'Optimize', true);
        
    matlabFunction(Me_edge_edge, 'File', fullfile(target_dir, 'Me_edge_edge'), ...
        'Vars', {P_mat, sigma_cond}, 'Optimize', true);
    
    matlabFunction(Ke_grad_grad, 'File', fullfile(target_dir, 'Ke_grad_grad'), ...
        'Vars', {P_mat, sigma_cond}, 'Optimize', true);
        
    matlabFunction(Ce_edge_grad, 'File', fullfile(target_dir, 'Ce_edge_grad'), ...
        'Vars', {P_mat, sigma_cond}, 'Optimize', true);
        
    % 新增: 计算 B 的辅助函数
    matlabFunction(B_vec, 'File', fullfile(target_dir, 'calc_element_B'), ...
        'Vars', {P_mat, A_local}, 'Optimize', true);

    fprintf('代码生成完成!\n');
end

function I = integrate_over_ref_tet(expr, u, v, w)
    I_w = int(expr, w, 0, 1 - u - v);
    I_v = int(I_w, v, 0, 1 - u);
    I = int(I_v, u, 0, 1);
end