% test_team7_coil_analysis.m
% 读取 Team 7 网格，分析线圈几何特征，计算并绘制电流矢量
% [Refactored] 使用 FileUtils 管理路径，代码更整洁
% ---------------------------------------------------------
clear; clc; close all;

% 1. 环境初始化
% 确保 src 和 tutorials 在路径中
base_paths = {'../src', 'magpackage/src', '../tutorials/TEAM7', 'magpackage/tutorials/TEAM7'};
for i = 1:length(base_paths)
    if exist(base_paths{i}, 'dir'), addpath(genpath(base_paths{i})); end
end

fprintf('==============================================\n');
fprintf('   TEAM 7 Coil Analysis & Visualization       \n');
fprintf('==============================================\n');

% 2. 加载网格
% 定义额外的搜索目录 (相对于脚本可能的运行位置)
searchDirs = {'../tutorials/TEAM7', 'magpackage/tutorials/TEAM7', '../../tutorials/TEAM7'};

try
    meshPath = FileUtils.locateFile('Team7.mphtxt', searchDirs);
    fprintf('[Step 1] Loading Mesh: %s...\n', meshPath);
    mesh = Mesh.load(meshPath, 'mm'); % 自动处理单位 (mm -> m)
    mesh.generateEdges();
catch ME
    error('Setup failed: %s', ME.message);
end

% 3. 区域合并 (Region Merging)
% Team 7 线圈由区域 3, 4, 5, 6 组成
coil_tags = [3, 4, 5, 6];
coil_id_merged = 999;

mask_coil = ismember(mesh.RegionTags, coil_tags);
if ~any(mask_coil)
    error('No elements found for coil regions [3, 4, 5, 6].');
end
mesh.RegionTags(mask_coil) = coil_id_merged;

fprintf('[Step 2] Merged regions %s into ID %d for analysis.\n', ...
    mat2str(coil_tags), coil_id_merged);

% 4. 几何特征分析
fprintf('[Step 3] Analyzing Geometry...\n');
[center, R_in, R_out, Lx, Ly, ~, ax_idx] = CoilGeometryUtils.autoDetectRoundedRect(mesh, coil_id_merged);

fprintf('   -> Analysis Result (Units: m):\n');
fprintf('      Center = [%.4f, %.4f, %.4f]\n', center);
fprintf('      R_in   = %.4f m\n', R_in);
fprintf('      R_out  = %.4f m\n', R_out);
fprintf('      Lx (Straight Dir 1) = %.4f m\n', Lx);
fprintf('      Ly (Straight Dir 2) = %.4f m\n', Ly);

% 5. 计算电流矢量
fprintf('[Step 4] Computing Current Vectors...\n');
dir_map = CoilGeometryUtils.computeRoundedRectDirection(mesh, coil_id_merged, center, Lx, Ly, ax_idx);

% 6. 可视化
fprintf('[Step 5] Visualizing...\n');

dof = DofHandler(mesh);
asm = Assembler(mesh, dof);
post = PostProcessor(asm);
viz = Visualizer(post);

figure('Name', 'TEAM 7 Coil Analysis', 'Position', [100, 100, 1000, 800]);

% 绘制线圈外表面 (半透明铜色)
viz.plotMesh('RegionTags', coil_id_merged, ...
    'FaceColor', [0.85 0.5 0.2], 'Alpha', 0.4, 'EdgeColor', 'none', 'Style', 'surface'); 
hold on;

% 绘制网格线框 (辅助观察结构)
viz.plotMesh('RegionTags', coil_id_merged, ...
    'Style', 'wireframe', 'EdgeColor', 'k', 'Alpha', 0.1);

% 绘制电流矢量
viz.plotElementVectors(dir_map, ...
    'MaxArrows', 1500, 'Color', 'r', 'Scale', 0.7);

title(sprintf('TEAM 7 (m): Lx=%.3f, Ly=%.3f, R_{in}=%.3f, R_{out}=%.3f', ...
    Lx, Ly, R_in, R_out));
view(3); axis equal; grid on; axis vis3d;
lighting gouraud; camlight;

fprintf('[Done] Analysis completed.\n');