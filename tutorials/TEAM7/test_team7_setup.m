function test_team7_setup()
% TEST_TEAM7_SETUP 验证 TEAM 7 物理模型建立是否正确
%
% 修改记录:
% 1. 适配 plot_mesh_domain 无输出参数的情况。
% 2. 增加手动 hold on 管理，防止图像被覆盖。
% 3. 使用虚拟句柄 (Dummy Objects) 创建图例。
% 4. 验证铝板(Region 2)和线圈(Region 3-6)的网格位置与几何路径(Coil.P1)的重合度。
    
    fprintf('正在调用 setup_physics...\n');
    try
        [Model, Coil] = team7_setup_physics();
    catch ME
        fprintf('错误: 无法建立模型。\n%s\n', ME.message);
        return;
    end

    % 创建图形窗口
    figure('Name', 'TEAM 7 Geometry Check', 'Color', 'w', 'Position', [100, 100, 1000, 700]);
    
    % --- 1. 绘制铝板 (Region 2) ---
    fprintf('绘制铝板 (Region 2)...\n');
    % NewFigure=false 表示在当前窗口绘制
    % 注意: plot_mesh_domain 内部结束时会 hold off，所以下一行必须重新 hold on
    plot_mesh_domain(Model.Mesh, 2, ...
        struct('FaceColor','c', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'NewFigure', false, 'Title', 'TEAM 7 Geometry Check'));
    hold on; 
    
    % --- 2. 绘制网格中的线圈区域 (Region 3,4,5,6) ---
    fprintf('绘制线圈网格域 (Region 3-6)...\n');
    coil_regions = [3, 4, 5, 6];
    existing_regions = unique(Model.Mesh.RegionTags);
    
    if all(ismember(coil_regions, existing_regions))
        plot_mesh_domain(Model.Mesh, coil_regions, ...
            struct('FaceColor','y', 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'NewFigure', false));
        hold on; % 再次保持，准备画线
    else
        warning('未在网格中找到线圈区域 ID (3-6)。跳过绘制网格线圈。');
    end

    % --- 3. 绘制几何线圈路径 (激励源) ---
    fprintf('绘制几何线圈路径...\n');
    % 这是 Biot-Savart 积分使用的实际电流路径
    hGeomCoil = plot3(Coil.P1(1,:), Coil.P1(2,:), Coil.P1(3,:), ...
        'r-', 'LineWidth', 2);
    
    % --- 4. 绘制参考边框 ---
    rectangle('Position',[0 0 0.294 0.294], 'EdgeColor','b', 'LineStyle','--', 'LineWidth', 1);

    % --- 5. 构建图例 (由于 plot_mesh_domain 无返回，使用虚拟对象) ---
    % 创建不可见的 patch 对象用于图例显示
    hLegend_Plate = fill3(nan, nan, nan, 'c', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hLegend_MeshCoil = fill3(nan, nan, nan, 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    legend([hLegend_Plate, hLegend_MeshCoil, hGeomCoil], ...
           {'Al Plate (Mesh)', 'Coil (Mesh Region)', 'Coil Source (Biot-Savart)'}, ...
           'Location', 'best');

    % --- 6. 视图设置 ---
    axis equal; grid on;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    view(3); % 3D 视图
    
    fprintf('完成。请检查:\n');
    fprintf('1. 红色线条 (几何源) 应位于 黄色块 (网格线圈) 的中心。\n');
    fprintf('2. 黄色块应悬浮于 青色板 上方。\n');
end