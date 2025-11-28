function run_team7()
% RUN_TEAM7 [Final Debug Run - Corrected Probe]
% 验证 TEAM 7 基准案例结果

    startup_fem;
    fprintf('==============================================\n');
    fprintf('   TEAM 7 Benchmark (Final Debug Run)         \n');
    fprintf('==============================================\n');

    % 1. 模型初始化
    [Model, Coil] = team7_setup_physics();
    
    % [DEBUG] 源场自检探针
    fprintf('\n  [Probe] 正在检查源场 (Biot-Savart) 强度...\n');
    % 选取线圈正下方一点 (0, 0, 50mm)
    ProbePt = [0; 0; 0.05]; 
    
    % [FIX] 参数顺序必须是 (Coil, Points)
    Bs_probe = compute_biot_savart_B_serial(Coil, ProbePt); 
    
    B_mag_T = norm(Bs_probe);
    fprintf('          探针点 (0,0,50mm) B_source = %.4f Gauss (预期值约 100-500)\n', B_mag_T * 1e4);
    
    if B_mag_T < 1e-6
        error('源场计算结果为 0！请检查 Coil.I 是否正确赋值。');
    end
    
    % 2. 求解 (50 Hz)
    FreqParams.Frequency = 50;
    fprintf('\n  [Solver] Solving Frequency Domain (50Hz)...\n');
    
    [Sol, Info] = solve_frequency_linear(Model, Coil, FreqParams);
    
    % 3. 结果对比
    compare_results_full(Model, Sol, Coil, Info);
end

function compare_results_full(Model, Sol, Coil, Info)
    % TEAM 7 Benchmark Data (mm)
    X_ref_mm = [0, 18, 36, 54, 72, 90, 108, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288]; 
    
    % Reference Values (Gauss)
    Bz_A1_Re = [-1.16, -18.46, 4.15, -21.59, -16.09, 0.23, 1.89, 75.53, 63.42, 14.15, 48.66, 12.40, 48.31, 51.26, 12.66, 46.11, 24.96];
    Bz_A1_Im = [3.63, 2.84, -23.62, -20.19, 3.07, 0.36, 44.35, 4.97, 71.55, 53.20, 13.04, 47.31, 12.05, 12.27, 53.61, 9.96, 2.36];
    Bz_A2_Re = [-1.83, -8.50, -13.60, -15.21, -14.48, -5.62, 28.77, 60.34, 61.84, 56.64, 53.40, 52.36, 53.93, 56.82, 59.48, 52.08, 26.56];
    Bz_A2_Im = [-1.63, -0.60, -0.43, 0.11, 1.26, 3.40, 6.53, 10.25, 11.83, 11.83, 11.01, 10.58, 10.80, 10.54, 10.62, 9.03, 1.79];
    
    X_J_mm = [0, 18, 126, 144, 162, 180, 198, 216, 234, 252, 270, 288];
    Jy_Top_Re = [-0.629, 0.794, -0.593, -0.229, -0.101, -0.011, 0.087, -0.000, 0.008, 0.043, 0.050, -0.687];
    Jy_Top_Im = [0.427, -0.873, -1.304, -0.035, 0.005, -0.001, 0.007, 0.182, 0.135, 0.033, 0.116, 0.855];
    Jy_Bot_Re = [0.461, 0.621, 1.573, 0.556, 0.237, 0.097, -0.034, -0.157, -0.305, -0.478, -0.660, -1.217];
    Jy_Bot_Im = [-0.662, 1.597, 4.163, -0.757, 0.672, -0.149, -0.050, 0.154, -0.749, -1.205, -1.575, -2.583];

    Line_A1 = [X_ref_mm/1000; repmat(0.072,1,17); repmat(0.034,1,17)];
    Line_A2 = [X_ref_mm/1000; repmat(0.144,1,17); repmat(0.034,1,17)];
    Line_J3 = [X_J_mm/1000; repmat(0.072,1,12); repmat(0.019,1,12)];
    Line_J4 = [X_J_mm/1000; repmat(0.072,1,12); zeros(1,12)];

    fprintf('  [Post] Interpolating Bz A1-B1...\n');
    [B1_Re, B1_Im] = get_B_field(Model, Sol, Coil, Line_A1, 3);
    
    fprintf('  [Post] Interpolating Bz A2-B2...\n');
    [B2_Re, B2_Im] = get_B_field(Model, Sol, Coil, Line_A2, 3);
    
    fprintf('  [Post] Interpolating Jy Top (A3-B3)...\n');
    [J3_Re, J3_Im] = get_J_field(Model, Sol, Coil, Line_J3, 2, Info);
    
    fprintf('  [Post] Interpolating Jy Bot (A4-B4)...\n');
    [J4_Re, J4_Im] = get_J_field(Model, Sol, Coil, Line_J4, 2, Info);

    figure('Name', 'TEAM 7 Comparison', 'Position', [50, 50, 1200, 800], 'Color', 'w');
    subplot(2,2,1); plot_compare(X_ref_mm, Bz_A1_Re, Bz_A1_Im, B1_Re*1e4, B1_Im*1e4, 'A1-B1 Bz', 'Gauss');
    subplot(2,2,2); plot_compare(X_ref_mm, Bz_A2_Re, Bz_A2_Im, B2_Re*1e4, B2_Im*1e4, 'A2-B2 Bz', 'Gauss');
    subplot(2,2,3); plot_compare(X_J_mm, Jy_Top_Re, Jy_Top_Im, J3_Re*1e-6, J3_Im*1e-6, 'Jy Top (A3-B3)', 'MA/m^2');
    subplot(2,2,4); plot_compare(X_J_mm, Jy_Bot_Re, Jy_Bot_Im, J4_Re*1e-6, J4_Im*1e-6, 'Jy Bot (A4-B4)', 'MA/m^2');
end

function plot_compare(X, RefRe, RefIm, FemRe, FemIm, Txt, YUnit)
    plot(X, RefRe, 'ks', 'MarkerFaceColor','k', 'MarkerSize', 6); hold on;
    plot(X, RefIm, 'kd', 'MarkerSize', 6);
    plot(X, FemRe, 'r-', 'LineWidth', 2);
    plot(X, FemIm, 'b--', 'LineWidth', 2);
    title(Txt); ylabel(YUnit); xlabel('x (mm)'); grid on;
    legend({'Ref Re','Ref Im','FEM Re','FEM Im'}, 'Location','best', 'FontSize', 9);
end

function [Val_0, Val_90] = get_B_field(Model, Sol, Coil, Pts, idx)
    Bs = compute_biot_savart_B_serial(Coil, Pts); 
    Mesh = Model.Mesh; Br = zeros(3, size(Pts, 2));
    Centers = (Mesh.P(:, Mesh.T(1,:)) + Mesh.P(:, Mesh.T(2,:)) + Mesh.P(:, Mesh.T(3,:)) + Mesh.P(:, Mesh.T(4,:))) / 4;
    for k = 1:size(Pts, 2)
        [~, e_id] = min(sum((Centers - Pts(:,k)).^2, 1));
        edges = Mesh.T2E(:, e_id); signs = double(Mesh.T2E_Sign(:, e_id));
        a_local = Sol.A(edges) .* signs;
        Br(:, k) = calc_element_B(Mesh.P(:, Mesh.T(:, e_id)), a_local); 
    end
    B_tot = Bs + Br; Val_0 = real(B_tot(idx, :)); Val_90 = -imag(B_tot(idx, :));
end

function [Val_0, Val_90] = get_J_field(Model, Sol, Coil, Pts, idx, Info)
    omega = 2 * pi * 50; sigma = 3.526e7;
    Mesh = Model.Mesh; 
    DofData = Info.DofData; 
    V_offset = DofData.NumEdges;
    Centers = (Mesh.P(:, Mesh.T(1,:)) + Mesh.P(:, Mesh.T(2,:)) + Mesh.P(:, Mesh.T(3,:)) + Mesh.P(:, Mesh.T(4,:))) / 4;
    As_all = compute_biot_savart_A(Coil, Pts); 
    vals = zeros(1, size(Pts, 2));
    for k = 1:size(Pts, 2)
        [~, e_id] = min(sum((Centers - Pts(:,k)).^2, 1));
        if Mesh.RegionTags(e_id) ~= 2, vals(k)=0; continue; end
        edges = Mesh.T2E(:, e_id); signs = double(Mesh.T2E_Sign(:, e_id));
        a_local = Sol.A(edges) .* signs;
        nodes = Mesh.T(:, e_id); p_elem = Mesh.P(:, nodes);
        invJ = inv([p_elem(:,2)-p_elem(:,1), p_elem(:,3)-p_elem(:,1), p_elem(:,4)-p_elem(:,1)])';
        G = invJ * [-1 1 0 0; -1 0 1 0; -1 0 0 1];
        L = [0.25;0.25;0.25;0.25];
        Ar_vec = zeros(3,1); pairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
        for i=1:6, na=pairs(i,1); nb=pairs(i,2); Ni = L(na)*G(:,nb) - L(nb)*G(:,na); Ar_vec = Ar_vec + Ni * a_local(i); end
        v_dofs = DofData.Node2DoF(nodes); V_elem = Sol.V_active(v_dofs - V_offset);
        J_vec = sigma * (-1i * omega * (Ar_vec + As_all(:,k)) - G * V_elem);
        vals(k) = J_vec(idx);
    end
    Val_0 = real(vals); Val_90 = -imag(vals);
end