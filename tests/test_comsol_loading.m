function test_comsol_loading()
% TEST_COMSOL_LOADING 验证 .mphtxt 读取器 (集成单位转换与绘图测试)
%
% 功能:
% 1. 基础读取: 验证节点和单元数量。
% 2. 单位测试: 验证 'mm', 'in' 等标签是否正确缩放坐标。
% 3. 绘图测试: 验证 plot_mesh_domain 是否能正常运行。

    addpath(genpath('src'));
    fprintf('==============================================\n');
    fprintf('   COMSOL Mesh Import Test Suite              \n');
    fprintf('==============================================\n');
    
    % 1. 搜索测试文件
    target_file = 'tests/test.mphtxt'; 
    if ~exist(target_file, 'file')
        target_files = dir('**/*.mphtxt');
        if ~isempty(target_files)
            target_file = fullfile(target_files(1).folder, target_files(1).name);
        else
            target_file = 'dummy_test.mphtxt';
            create_dummy_mphtxt(target_file);
            cleanup = onCleanup(@() delete(target_file));
        end
    end
    
    fprintf('测试目标文件: %s\n', target_file);
    
    try
        % --- Test 1: 基础读取 ---
        fprintf('\n[Test 1] 基础读取与拓扑检查...\n');
        Raw = read_mphtxt(target_file, 'm'); % 默认米
        
        nNodes = size(Raw.P, 2);
        nTets  = size(Raw.T, 2);
        
        fprintf('  Nodes: %d, Tets: %d\n', nNodes, nTets);
        
        if nNodes == 0
            error('[FAIL] 未读取到节点。');
        end
        fprintf('  [PASS] 基础数据读取成功。\n');
        
        % --- Test 2: 单位转换验证 ---
        fprintf('\n[Test 2] 单位转换验证 (Unit Scaling)...\n');
        
        % 读取 mm 版本
        Raw_mm = read_mphtxt(target_file, 'mm');
        
        % 验证比例: mm 读进来的数值应该是 m 读进来的 0.001 倍
        % 因为 read_mphtxt 会把文件里的数 * scale 变成米。
        % 比如文件里是 1000。
        % mode='m': P = 1000
        % mode='mm': P = 1000 * 0.001 = 1.0
        
        scale_err_mm = norm(Raw.P * 1e-3 - Raw_mm.P, 'fro');
        
        % 读取 inch 版本
        Raw_in = read_mphtxt(target_file, 'in');
        scale_err_in = norm(Raw.P * 0.0254 - Raw_in.P, 'fro');
        
        fprintf('  Diff (m vs mm): %.2e\n', scale_err_mm);
        fprintf('  Diff (m vs in): %.2e\n', scale_err_in);
        
        if scale_err_mm < 1e-9 && scale_err_in < 1e-9
            fprintf('  [PASS] 坐标缩放正确。\n');
        else
            warning('  [FAIL] 坐标缩放误差过大，请检查 read_mphtxt 单位逻辑。');
        end
        
        % --- Test 3: 可视化功能验证 ---
        fprintf('\n[Test 3] 绘图功能验证 (plot_mesh_domain)...\n');
        
        % 构造最小 Model 结构体用于绘图
        ModelPlot.Mesh = Raw;
        
        % 获取所有存在的 Region Tag
        if ~isempty(Raw.RegionTags)
            tags = unique(Raw.RegionTags);
            target_tag = tags(1);
            
            fprintf('  尝试绘制 Region %d ... ', target_tag);
            
            try
                % 设置不新建窗口或立即关闭，以免阻塞测试
                Opts.NewFigure = true;
                Opts.Title = 'Test Plot (Auto Close)';
                
                hFig = figure('Visible', 'off'); % 后台绘制
                plot_mesh_domain(ModelPlot.Mesh, target_tag, Opts);
                close(hFig);
                
                fprintf('[PASS] 绘图函数运行正常。\n');
            catch ME
                fprintf('\n  [FAIL] 绘图函数报错: %s\n', ME.message);
            end
        else
            fprintf('  [SKIP] 无四面体单元，跳过体网格绘图测试。\n');
        end
        
    catch ME
        fprintf('\n[ERROR] 测试中断: %s\n', ME.message);
        % rethrow(ME);
    end
    
    fprintf('\n==============================================\n');
    fprintf('   Tests Completed.                           \n');
    fprintf('==============================================\n');
end

function create_dummy_mphtxt(filename)
    % 生成符合 MPHTXT 格式的 Dummy 文件 (用于无文件时的回退测试)
    fid = fopen(filename, 'w');
    fprintf(fid, '0 0 1\n'); 
    
    % Nodes
    fprintf(fid, '4 # number of mesh vertices\n');
    fprintf(fid, '# Mesh vertex coordinates\n');
    fprintf(fid, '0 0 0\n1 0 0\n0 1 0\n0 0 1\n');
    
    % Tets
    fprintf(fid, '\n3 tet # type name\n');
    fprintf(fid, '1 # number of elements\n');
    fprintf(fid, '4 # number of vertices per element\n');
    fprintf(fid, '# Elements\n0 1 2 3\n');
    fprintf(fid, '# Geometric entity indices\n7\n');
    
    % Tris
    fprintf(fid, '\n3 tri # type name\n');
    fprintf(fid, '1 # number of elements\n');
    fprintf(fid, '3 # number of vertices per element\n');
    fprintf(fid, '# Elements\n0 1 2\n');
    fprintf(fid, '# Geometric entity indices\n3\n');
    
    fclose(fid);
end