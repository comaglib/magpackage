function test_transient_voltage()
% TEST_VOLTAGE_DRIVEN 测试电压驱动求解器
% 案例: 铁心电感器的电压阶跃响应 (RL 充电曲线)

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   Voltage Driven Test (RL Circuit)           \n');
    fprintf('==============================================\n');
    
    % 1. 网格 (复用 Iron-Air Mesh)
    msh_file = 'test_voltage.msh';
    create_iron_air_mesh_v2(msh_file); % 需确保此函数存在
    cleanup = onCleanup(@() delete(msh_file));
    
    Raw = read_msh(msh_file);
    Model.Mesh = build_topology(Raw);
    
    % 2. 材料 (铁心)
    Model.Materials.Lib(1).Name = 'Air';
    Model.Materials.Lib(1).Type = 'Linear';
    Model.Materials.Lib(1).Mu_r = 1.0; Model.Materials.Lib(1).Sigma = 0;
    
    Model.Materials.Lib(2).Name = 'Iron';
    Model.Materials.Lib(2).Type = 'Linear';
    Model.Materials.Lib(2).Mu_r = 1000.0; Model.Materials.Lib(2).Sigma = 0; % 无涡流
    
    Model.Materials.ActiveMap = Model.Mesh.RegionTags;
    
    % 3. 线圈
    CoilParams.Center = [0, 0, 0];
    CoilParams.Normal = [0, 0, 1];
    CoilParams.Length = 0.5; CoilParams.Radius = 0.2; 
    CoilParams.Current = 0; CoilParams.N_seg = 36;
    Coil = create_racetrack_coil(CoilParams);
    
    % 4. 电路参数
    CircuitProps.R = 1.0;      % 1 Ohm
    CircuitProps.L_leak = 1e-3; % 1 mH (漏感/空心电感)
    CircuitProps.V_func = @(t) 10.0 * (t>0); % 10V Step
    
    % 5. 时间设置
    % 预计 L_total = L_leak + L_iron. L_iron 应该比较大
    % 假设 L_tot ~ 10 mH. Tau = L/R ~ 0.01s.
    TimeParams.dt = 0.001; 
    TimeParams.N_steps = 50;
    
    % 6. 边界
    fixed_edges = identify_boundary_edges(Model.Mesh, Raw, 100);
    Model.Runtime.FixedEdges = fixed_edges;
    
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % 7. 求解
    [Sol, ~] = solve_transient_voltage(Model, Coil, CircuitProps, TimeParams);
    
    % 8. 验证
    I_final = Sol.I_history(end);
    I_theoretical = 10.0 / 1.0; % U/R
    
    fprintf('\nFinal Current: %.4f A (Target: %.4f A)\n', I_final, I_theoretical);
    
    % 绘制曲线
    t_vec = (1:TimeParams.N_steps) * TimeParams.dt;
    figure; plot(t_vec, Sol.I_history, 'o-'); grid on;
    title('Current Step Response'); xlabel('Time (s)'); ylabel('Current (A)');
    
    if abs(I_final - I_theoretical) < 0.5
        fprintf('[PASS] Current converges to U/R.\n');
    else
        fprintf('[WARN] Current did not reach steady state (Inductance too large?).\n');
    end
    
    % 验证电感效应
    % 如果只有 L_leak, tau = 1ms. 50步(50ms)应该早就稳态了。
    % 如果电流上升缓慢，说明 L_iron 起作用了。
    I_1tau = Sol.I_history(1); % t=1ms
    % I(1ms) for 1mH = 10*(1 - exp(-1)) = 6.32 A.
    fprintf('Current at t=1ms: %.4f A\n', I_1tau);
    
    if I_1tau < 6.0
        fprintf('[PASS] Slow rise detected. Iron core inductance is active.\n');
    end
end

function create_iron_air_mesh_v2(filename)
    % 简单的铁心网格 (ID X Y Z 格式)
    fid = fopen(filename, 'w');
    fprintf(fid, '$MeshFormat\n4.1 0 8\n$EndMeshFormat\n');
    fprintf(fid, '$Nodes\n1 5 1 5\n2 1 0 5\n');
    fprintf(fid, '1 0 0 0\n2 0.1 0 0\n3 0 0.1 0\n4 0 0 0.1\n5 0.2 0.2 0.2\n');
    fprintf(fid, '$EndNodes\n');
    fprintf(fid, '$Elements\n3 3 1 3\n');
    fprintf(fid, '2 100 2 1\n1 2 3 5\n');
    fprintf(fid, '3 2 4 1\n2 1 2 3 4\n'); % Iron (Tag 2)
    fprintf(fid, '3 1 4 1\n3 2 3 4 5\n'); % Air (Tag 1)
    fprintf(fid, '$EndElements\n');
    fclose(fid);
end