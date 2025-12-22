classdef PostProcessor < handle
    % POSTPROCESSOR 有限元后处理模块 (v3.5 - Matrix Source Support)
    % 
    % 描述:
    %   该类封装了有限元分析(FEM)求解后的所有数据处理逻辑。
    %   它负责将求解得到的自由度系数(DOF coefficients)转换为具有物理意义的场量
    %   (如 B 场、J 场、能量等)，并提供点探测(Probe)和空间插值功能。
    %
    % 主要功能:
    %   1. 场量计算: 
    %      - computeElementB: 计算单元中心的磁通密度 B = curl(A)。
    %      - computeElementJ: 计算单元中心的电流密度 (包含源电流和涡流)。
    %      - computeMagneticEnergy: 计算全场磁能 (支持非线性BH曲线)。
    %   2. 场量探测 (Probe): 
    %      - 支持任意空间点的 B 场、A 场及总电流密度 J 的插值查询。
    %      - 支持沿线段采样 (probeLine)。
    %   3. 数据平滑: 
    %      - mapElementsToNodes: 将单元常数场平滑映射到节点，用于云图绘制。
    %        特别处理了材料分界面处的场不连续性。
    %
    % 更新日志:
    %   v3.5: probeTotalCurrent 支持传入全场源电流矩阵 (3 x Ne)，
    %         可直接处理空间变化的电流分布（如自动生成的线圈电流场）。
    
    properties
        Assembler       % 组装器对象，包含刚度矩阵组装逻辑
        Mesh            % 网格对象，包含 P(节点), T(拓扑), RegionTags(区域标记)
        DofHandler      % 自由度处理器，管理局部到全局的自由度映射
        Triangulation   % MATLAB内置 triangulation 对象，用于加速空间点定位 (pointLocation)
    end
    
    methods
        function obj = PostProcessor(assembler)
            % 构造函数
            % 输入:
            %   assembler - 包含 Mesh 和 DofHandler 的组装器实例
            
            obj.Assembler = assembler;
            obj.Mesh = assembler.Mesh;
            obj.DofHandler = assembler.DofHandler;
            
            % 初始化三角剖分对象 (Triangulation)
            % 这是一个非常重要的步骤，后续所有的空间点探测 (Probe) 
            % 都依赖于 pointLocation 函数来快速查找点所在的四面体单元。
            if ~isempty(obj.Mesh.T)
                T_double = double(obj.Mesh.T'); % 转置为 [Ne x 4]
                P_double = double(obj.Mesh.P'); % 转置为 [Np x 3]
                obj.Triangulation = triangulation(T_double, P_double);
            end
        end
        
        function centers = computeCentroids(obj)
            % COMPUTECENTROIDS 计算所有四面体单元的几何中心
            % 输出:
            %   centers - [3 x Ne] 矩阵，每列为一个单元的重心坐标
            
            P = obj.Mesh.P; 
            T = obj.Mesh.T;
            % 重心坐标公式: (P1 + P2 + P3 + P4) / 4
            centers = (P(:, T(1,:)) + P(:, T(2,:)) + P(:, T(3,:)) + P(:, T(4,:))) / 4.0;
        end
        
        % POST
        function B_elems = computeElementB(obj, A_sol, spaceName)
            % COMPUTEELEMENTB 计算每个单元中心的磁通密度 B = curl(A)
            % 
            % 输入: 
            %   A_sol - [Ndofs x 1] 或 [Ndofs x Ncols] 磁矢位解向量
            %   spaceName - (可选) 有限元空间名称，默认取第一个 (通常是 'Nedelec_P1')
            % 输出:
            %   B_elems - [3 x Ne x Ncols] 每个单元中心的 B 矢量
            
            if nargin < 3
                keys = obj.DofHandler.DofMaps.keys; 
                spaceName = keys{1}; 
            end
            
            dofMap = obj.DofHandler.DofMaps(spaceName);
            
            numElems = fix(full(double(obj.Mesh.NumElements)));
            numCols = size(A_sol, 2); % 支持多列解 (如多谐波分析或多时间步结果)
            B_elems = zeros(3, numElems, numCols);
            
            % --- 并行计算准备 ---
            % 使用 Parallel Pool Constant 优化大数据传输。
            % 如果直接在 parfor 中使用大数组 (如 obj.Mesh.P)，
            % MATLAB 会在每次迭代或每个 worker 中复制数据，导致巨大开销。
            C_P = parallel.pool.Constant(obj.Mesh.P);
            C_T = parallel.pool.Constant(obj.Mesh.T);
            C_Signs = parallel.pool.Constant(double(obj.Mesh.T2E_Sign)); % 棱边方向修正符号
            C_Sol = parallel.pool.Constant(A_sol);
            C_DofMap = parallel.pool.Constant(dofMap);
            
            % 预计算参考单元的 Curl 值 (常数, 任意点即可，因为 Nedelec P1 的 curl 是常数)
            q_center = [0.25; 0.25; 0.25];
            [~, curl_ref] = nedelec_tet_p1(q_center); 
            C_ref = reshape(curl_ref, 6, 3);
            
            % --- 并行循环计算 ---
            parfor e = 1:numElems
                loc_P = C_P.Value; loc_T = C_T.Value; loc_Signs = C_Signs.Value;
                loc_Sol = C_Sol.Value; loc_DofMap = C_DofMap.Value;
                
                % 1. 提取单元节点坐标
                nodes = loc_T(:, e); 
                p_elem = loc_P(:, nodes); 
                
                % 2. 计算雅可比矩阵
                % 这是从参考单元(Reference)映射到物理单元(Physical)的关键
                [J_mat, detJ] = compute_jacobian_tet(p_elem); 
                invDetJ = 1.0 / detJ;
                
                % 3. Piola 变换 (H-curl 空间的变换法则)
                % curl_phy = (1/detJ) * J * curl_ref
                % 这里计算物理空间下的基函数旋度矩阵
                C_phy = (C_ref * J_mat') * invDetJ;
                
                % 4. 提取局部自由度系数
                dofs = loc_DofMap(:, e); 
                s = loc_Signs(:, e); % 棱边方向符号修正
                A_local_all = loc_Sol(dofs, :) .* s;
                
                % 5. 叠加计算 B = Sum( Ni_curl * ai )
                B_elems(:, e, :) = C_phy' * A_local_all;
            end
        end
        
        function J_elems = computeElementJ(obj, A_sol, V_sol, omega, sigmaMap, sourceMap)
            % COMPUTEELEMENTJ 计算每个单元的电流密度 J = Js + J_eddy
            % 
            % 物理公式:
            %   J_total = J_source + J_eddy
            %   J_eddy  = sigma * E = sigma * (-j*w*A - grad(V))
            %
            % 输入:
            %   A_sol: 磁矢位解向量
            %   V_sol: 电标量位解向量 (可选，若无则为 [])
            %   omega: 角频率 (2*pi*f)
            %   sigmaMap: 导电率映射 (Map: Tag->Value) 或 向量
            %   sourceMap: 源电流映射 (Map) 或 全场电流矩阵 [3 x Ne]
            % 输出:
            %   J_elems: [3 x Ne] 每个单元的电流密度矢量
            
            numElems = fix(full(double(obj.Mesh.NumElements)));
            J_elems = zeros(3, numElems);
            
            % --- 1. 准备并行常量 (减少通信开销) ---
            C_P = parallel.pool.Constant(obj.Mesh.P);
            C_T = parallel.pool.Constant(obj.Mesh.T);
            C_Signs = parallel.pool.Constant(double(obj.Mesh.T2E_Sign));
            C_Tags = parallel.pool.Constant(obj.Mesh.RegionTags);
            
            % Nedelec 自由度映射 (用于 A 场)
            dofMapA = obj.DofHandler.DofMaps('Nedelec_P1');
            C_DofA = parallel.pool.Constant(dofMapA);
            C_SolA = parallel.pool.Constant(A_sol);
            
            % Lagrange 自由度映射 (用于 V 场，可选)
            hasV = ~isempty(V_sol);
            if hasV
                dofMapV = obj.DofHandler.DofMaps('Lagrange_P1');
                C_DofV = parallel.pool.Constant(dofMapV);
                C_SolV = parallel.pool.Constant(V_sol);
                
                % 获取 Lagrange 空间在全局解向量中的偏移量
                if obj.DofHandler.SpaceOffsets.isKey('Lagrange_P1')
                    offset_V = obj.DofHandler.SpaceOffsets('Lagrange_P1');
                else
                    offset_V = 0;
                end
            else
                offset_V = 0;
                C_DofV = []; C_SolV = [];
            end
            
            % 准备材料与源 (处理 Map/Matrix/Vector 多种输入格式)
            % 处理导电率 Sigma
            if isa(sigmaMap, 'containers.Map')
                C_SigmaMap = parallel.pool.Constant(sigmaMap);
                useSigmaMap = true;
                C_SigmaVec = [];
            else
                C_SigmaVec = parallel.pool.Constant(sigmaMap);
                useSigmaMap = false;
                C_SigmaMap = [];
            end
            
            % 处理源电流 Source
            useSourceMatrix = (isnumeric(sourceMap) && size(sourceMap, 2) == numElems);
            if useSourceMatrix
                C_SourceMat = parallel.pool.Constant(sourceMap); % 矩阵形式 [3xNe]
                C_SourceMap = [];
            else
                C_SourceMap = parallel.pool.Constant(sourceMap); % Map形式
                C_SourceMat = [];
            end
            
            % --- 2. 预计算参考单元基函数 (中心点) ---
            q_center = [0.25; 0.25; 0.25];
            [val_ref_A, ~] = nedelec_tet_p1(q_center); % Nedelec 基函数值 [6 x 3]
            [~, grad_ref_V] = lagrange_tet_p1(q_center); % Lagrange 梯度值 [4 x 3]
            
            C_RefA = parallel.pool.Constant(val_ref_A);
            C_RefGradV = parallel.pool.Constant(grad_ref_V);
            
            % --- 3. 并行循环计算 ---
            parfor e = 1:numElems
                % 获取几何信息
                loc_P = C_P.Value; loc_T = C_T.Value; 
                nodes = loc_T(:, e); p_elem = loc_P(:, nodes);
                
                tag = C_Tags.Value(e);
                
                % A. 获取当前单元的 Sigma
                sigma = 0;
                if useSigmaMap
                    map = C_SigmaMap.Value;
                    if map.isKey(tag), sigma = map(tag); end
                else
                    vec = C_SigmaVec.Value;
                    if length(vec) >= e, sigma = vec(e); 
                    elseif length(vec) >= tag, sigma = vec(tag); end
                end
                
                % B. 获取当前单元的 Js (源电流)
                Js = [0;0;0];
                if useSourceMatrix
                    src = C_SourceMat.Value;
                    Js = src(:, e);
                else
                    map = C_SourceMap.Value;
                    if map.isKey(tag), val = map(tag); Js = val(:); end
                end
                
                % C. 计算涡流 J_eddy
                J_eddy = [0;0;0];
                if abs(sigma) > 1e-12 % 仅在导电区域计算
                    [J_mat, ~] = compute_jacobian_tet(p_elem);
                    invJt = inv(J_mat)';
                    
                    % C1. 计算 A 场值 (Nedelec插值)
                    loc_SolA = C_SolA.Value; loc_DofA = C_DofA.Value; loc_Signs = C_Signs.Value;
                    dofs_A = loc_DofA(:, e);
                    s_A = loc_Signs(:, e);
                    A_coeffs = loc_SolA(dofs_A) .* double(s_A);
                    
                    % H(curl) 变换: A = J^-T * A_ref (注意: A 是协变矢量)
                    ref_A_val = C_RefA.Value' * A_coeffs; % [3x6]*[6x1] = [3x1]
                    A_val = invJt * ref_A_val;
                    
                    % C2. 计算 V 场梯度 (Lagrange梯度)
                    GradV_val = [0;0;0];
                    if hasV
                        loc_SolV = C_SolV.Value; loc_DofV = C_DofV.Value;
                        dofs_V = loc_DofV(:, e);
                        local_dofs_V = dofs_V - offset_V;
                        
                        % 确保索引有效 (空气中可能无V自由度，取决于编号策略)
                        if all(local_dofs_V > 0) && all(local_dofs_V <= length(loc_SolV))
                            V_coeffs = loc_SolV(local_dofs_V);
                            
                            % H(grad) 变换: gradV = J^-T * gradV_ref
                            ref_grad_val = C_RefGradV.Value' * V_coeffs; % [3x4]*[4x1] = [3x1]
                            GradV_val = invJt * ref_grad_val;
                        end
                    end
                    
                    % 欧姆定律: J_eddy = sigma * (-jwA - gradV)
                    J_eddy = sigma * (-1j * omega * A_val - GradV_val);
                end
                
                % D. 总电流叠加
                J_elems(:, e) = Js + J_eddy;
            end
        end
        
        function RegionDataMap = mapElementsToNodes(obj, Val_elems)
            % MAPELEMENTSTONODES_SEPARATED 按材料区域独立平滑单元值到节点
            %
            % 原理: 
            %   有限元计算得到的 B 场和 J 场是"基于单元"的常数(P0)或分片函数。
            %   在不同材料界面处(如空气/铁，空气/导体)，物理场通常是不连续的(切向或法向分量)。
            %   如果直接对所有围绕节点的单元取平均，会导致边界处的数值被错误地"抹平"。
            %   本函数遍历所有 RegionTag，仅在同一材料区域内部进行节点平均平滑。
            %
            % 输入:
            %   Val_elems: [dims x NumElements] 单元中心值矩阵
            %
            % 输出:
            %   RegionDataMap: containers.Map
            %     - Key: Region Tag (double)
            %     - Value: [dims x NumNodes] 矩阵。
            %              其中属于该 Region 的节点存有平滑后的值，
            %              不属于该 Region 的节点值为 0。
            %              (保留全尺寸矩阵是为了方便直接使用全局节点索引进行绘图)
            
            mesh = obj.Mesh;
            numNodes = mesh.NumNodes;
            [dims, numElems] = size(Val_elems);
            
            if numElems ~= mesh.NumElements
                error('PostProcessor:DimensionMismatch', ...
                    'Val_elems columns (%d) must match mesh elements (%d).', numElems, mesh.NumElements);
            end
            
            % 初始化输出 Map
            RegionDataMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            uniqueTags = unique(mesh.RegionTags);
            
            % 提取网格拓扑
            T = mesh.T;
            RegionTags = mesh.RegionTags;
            
            for k = 1:length(uniqueTags)
                tag = uniqueTags(k);
                
                % 1. 筛选属于当前 Region 的单元
                mask_elem = (RegionTags == tag);
                T_sub = T(:, mask_elem);       % [4 x N_sub]
                Val_sub = Val_elems(:, mask_elem); % [dims x N_sub]
                
                if isempty(T_sub), continue; end
                
                % 2. 准备累加器 (全尺寸)
                sum_vals = zeros(dims, numNodes);
                
                % 3. 累加单元值到节点 (Vectorized - 向量化操作加速)
                % 将单元值广播到其 4 个节点
                nodes_flat = T_sub(:); % [4*N_sub x 1] 展开所有节点索引
                
                for d = 1:dims
                    % 取出第 d 维分量，并复制 4 次以对应 4 个节点
                    v_d = Val_sub(d, :); 
                    vals_flat = repmat(v_d, 4, 1);
                    vals_flat = vals_flat(:);
                    
                    % 使用 accumarray 进行快速分组求和
                    % 累加到全局节点索引位置
                    sum_vals(d, :) = accumarray(nodes_flat, vals_flat, [numNodes, 1])';
                end
                
                % 4. 计算每个节点被该区域内多少个单元共享 (作为分母)
                ones_flat = ones(numel(T_sub), 1);
                count_nodes = accumarray(nodes_flat, ones_flat, [numNodes, 1])';
                
                % 5. 求平均 (Sum / Count)
                valid_mask = count_nodes > 0;
                region_node_vals = zeros(dims, numNodes);
                
                for d = 1:dims
                    % 仅计算有效节点，避免除以零
                    region_node_vals(d, valid_mask) = sum_vals(d, valid_mask) ./ count_nodes(valid_mask);
                end
                
                % 6. 存入 Map
                RegionDataMap(tag) = region_node_vals;
            end
        end
        
        function CombinedData = combineRegionData(obj, RegionDataMap, targetRegions)
            % COMBINEREGIONDATA 将多个区域的数据合并为一个全局数组
            %
            % 用途: 用于 Visualizer.plotFieldOnSurface 等需要全场节点数据的场合。
            % 输入:
            %   RegionDataMap: mapElementsToNodes 返回的 Map
            %   targetRegions: 需要合并的区域 Tag 数组
            % 输出:
            %   CombinedData: [dims x NumNodes] 合并后的数据矩阵
            
            mesh = obj.Mesh;
            numNodes = mesh.NumNodes;
            
            % 自动推断数据维度 (dims)
            dims = 1;
            for k = 1:length(targetRegions)
                tag = targetRegions(k);
                if RegionDataMap.isKey(tag)
                    dims = size(RegionDataMap(tag), 1);
                    break;
                end
            end
            
            CombinedData = zeros(dims, numNodes);
            
            % 遍历目标区域并合并数据
            for k = 1:length(targetRegions)
                tag = targetRegions(k);
                if RegionDataMap.isKey(tag)
                    vals = RegionDataMap(tag);
                    % 找出该区域有数据的节点 (非零)
                    if dims == 1
                        mask = (vals ~= 0);
                        CombinedData(mask) = vals(mask);
                    else
                        mask = (sum(abs(vals), 1) > 0);
                        CombinedData(:, mask) = vals(:, mask);
                    end
                end
            end
        end
        
        % PROBE
        function [B_val, elem_idx] = probeB(obj, A_sol, point)
            % PROBEB 计算空间任意点处的磁通密度 B = curl(A)
            % 输入:
            %   A_sol: 解向量
            %   point: 探测点坐标 [x, y, z]
            % 输出:
            %   B_val: B 场矢量
            %   elem_idx: 点所在的单元索引
            
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            % 1. 查找点所在的四面体单元
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                B_val = zeros(3, size(A_sol, 2)); % 点在网格外部
                return;
            end
            
            % 2. 复用内部逻辑进行插值
            dofMap = obj.DofHandler.DofMaps('Nedelec_P1');
            nodes = obj.Mesh.T(:, elem_idx); 
            p_elem = obj.Mesh.P(:, nodes);
            
            [J_mat, detJ] = compute_jacobian_tet(p_elem); 
            invDetJ = 1.0 / detJ;
            
            q_center = [0.25; 0.25; 0.25]; % P1 单元内 curl 为常数，坐标无关
            [~, curl_ref] = nedelec_tet_p1(q_center);
            C_ref = reshape(curl_ref, 6, 3);
            
            % Piola Transform
            C_phy = (C_ref * J_mat') * invDetJ;
            
            dofs = dofMap(:, elem_idx); 
            s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            A_local = A_sol(dofs, :) .* s;
            
            B_val = C_phy' * A_local;
        end
        
        function [A_val, elem_idx] = probeSolution(obj, A_sol, point)
            % PROBESOLUTION 计算空间任意点处的磁矢位 A (Nedelec 插值)
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                A_val = zeros(3, size(A_sol, 2));
                return;
            end
            
            [A_val, ~] = obj.internalProbeNedelec(A_sol, elem_idx, point, 'Nedelec_P1');
        end
        
        function [J_total, elem_idx] = probeTotalCurrent(obj, point, A_sol, V_sol, omega, sigmaMap, sourceMap)
            % PROBETOTALCURRENT 计算空间任意点处的总电流密度
            % 公式: J_total = J_source + J_eddy
            %              = J_source + sigma * (-j*omega*A - grad(V))
            
            % 1. 查找点所在的单元
            if isempty(obj.Triangulation), error('Triangulation not initialized.'); end
            elem_idx = pointLocation(obj.Triangulation, point(:)');
            
            if isnan(elem_idx)
                J_total = [0;0;0];
                return;
            end
            
            tag = obj.Mesh.RegionTags(elem_idx);
            
            % 2. 获取源电流 (J_source)
            % [Update] 支持全场矩阵 (3 x Ne) 或 传统的 Map
            if isnumeric(sourceMap) && size(sourceMap, 2) == obj.Mesh.NumElements
                % 如果是全场矩阵，直接根据单元索引取值 (适用于复杂线圈)
                J_s = sourceMap(:, elem_idx);
            else
                % 否则按区域 Tag 查找 (适用于简单块状电流)
                J_s = obj.getPropertyValue(sourceMap, tag, [0,0,0]);
            end
            J_s = J_s(:);
            
            % 3. 获取导电率
            sigma = obj.getPropertyValue(sigmaMap, tag, 0);
            
            % 4. 计算涡流 J_eddy
            if abs(sigma) > 1e-12
                % 仅在导电区域计算电场项，避免在空气中查找 Lagrange 自由度导致越界
                [A_val, ~] = obj.internalProbeNedelec(A_sol, elem_idx, point, 'Nedelec_P1');
                
                if ~isempty(V_sol)
                    GradV_val = obj.internalProbeLagrangeGrad(V_sol, elem_idx, 'Lagrange_P1');
                else
                    GradV_val = [0;0;0];
                end
                
                E_eddy = -1j * omega * A_val - GradV_val;
                J_total = J_s + sigma * E_eddy;
            else
                % 非导电区域 (如空气、线圈源区域)，没有涡流
                J_total = J_s;
            end
        end
        
        function result = probeLine(obj, A_sol, startPt, endPt, numPts, type, varargin)
            % PROBELINE 沿线段进行采样
            % 输入:
            %   type: 'B' (磁密), 'A' (磁矢位), 'J' (电流密度)
            %   varargin: 传递给 J 计算所需的额外参数 (V_sol, omega, sigmaMap, sourceMap)
            
            pts = [linspace(startPt(1), endPt(1), numPts)', ...
                   linspace(startPt(2), endPt(2), numPts)', ...
                   linspace(startPt(3), endPt(3), numPts)'];
            
            result = zeros(numPts, 3); 
            
            for k = 1:numPts
                pt = pts(k, :);
                if strcmpi(type, 'B')
                    val = obj.probeB(A_sol, pt);
                elseif strcmpi(type, 'A')
                    val = obj.probeSolution(A_sol, pt);
                elseif strcmpi(type, 'J')
                    % 对于 J，varargin 需包含 V_sol, omega, sigmaMap, sourceMap
                    val = obj.probeTotalCurrent(pt, A_sol, varargin{:});
                else
                    error('Unknown probe type');
                end
                result(k, :) = val.'; 
            end
        end
        
        % ENERGY
        function mag = computeMagnitude(~, VectorField)
            % COMPUTEMAGNITUDE 计算矢量场的模值
            % 输入: VectorField [3 x N x K]
            % 输出: mag [N x K]
            mag_sq = sum(abs(VectorField).^2, 1);
            mag = sqrt(mag_sq);
            [~, Ne, K] = size(VectorField);
            mag = reshape(mag, Ne, K);
        end
        
        function Energy = computeMagneticEnergy(obj, A_sol, MatLibData)
            % COMPUTEMAGNETICENERGY 计算全场磁能 (支持 Map 和 Array 输入)
            % 
            % 原理:
            %   如果全是线性材料，能量 E = 0.5 * Integral( B * H ) dV
            %   离散后等价于 E = 0.5 * A' * K * A (其中 K 为刚度矩阵)
            %   这比逐单元积分快得多。
            %
            %   如果包含非线性材料 (B-H 曲线)，刚度矩阵 K 不再适用能量计算，
            %   必须使用逐单元数值积分 W = Integral( w_density ) dV。
            
            % 1. 检查材料库类型并提取 Keys
            is_map = isa(MatLibData, 'containers.Map');
            if is_map
                matKeys = MatLibData.keys;
                numMats = length(matKeys);
            else
                numMats = length(MatLibData);
            end
            
            % 2. 检查是否全线性材料
            all_linear = true;
            for i = 1:numMats
                if is_map
                    mat = MatLibData(matKeys{i});
                else
                    mat = MatLibData(i);
                end
                
                if ~strcmp(mat.Type, 'Linear')
                    all_linear = false; break;
                end
            end
            
            space = FunctionSpace('Nedelec', 1); 
            
            if all_linear
                % --- 线性优化分支 (基于刚度矩阵) ---
                % 组装 NuMap (支持 Map 或 Vector 格式以匹配 Assembler)
                if is_map
                    NuMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
                    for i = 1:numMats
                        tag = matKeys{i};
                        NuMap(tag) = MatLibData(tag).Nu_Linear;
                    end
                else
                    NuMap = zeros(numMats, 1);
                    for i = 1:numMats
                        NuMap(i) = MatLibData(i).Nu_Linear;
                    end
                end
                
                % 组装线性刚度矩阵
                K = obj.Assembler.assembleStiffness(space, NuMap);
                % 矩阵乘法快速计算能量
                Energy = 0.5 * real(A_sol' * K * A_sol);
            else
                % --- 非线性通用分支 (基于单元积分) ---
                % 准备数据并调用积分核函数
                packedData = obj.Assembler.preparePackedData(space);
                packedData.RegionTags = obj.Mesh.RegionTags;
                % integrate_energy_kernel 内部支持 Map 类型的 MatLibData
                Energy = obj.integrate_energy_kernel(packedData, A_sol, MatLibData);
            end
        end
    end
    
    methods (Access = private)
        function [val, J_mat] = internalProbeNedelec(obj, solVec, elem_idx, point, spaceName)
            % internalProbeNedelec 封装 Nedelec 单元内的插值逻辑
            dofMap = obj.DofHandler.DofMaps(spaceName);
            nodes = obj.Mesh.T(:, elem_idx); 
            p_elem = obj.Mesh.P(:, nodes);
            
            [J_mat, ~] = compute_jacobian_tet(p_elem);
            
            % 物理坐标 -> 参考坐标
            % x_ref = J^-1 * (x_phy - P0)
            x_ref = J_mat \ (point(:) - p_elem(:,1));
            
            [val_ref, ~] = nedelec_tet_p1(x_ref); 
            
            % H(curl) Piola 变换: u = J^{-T} * u_ref
            invJt = inv(J_mat)'; 
            N_phy = val_ref * invJt'; 
            
            dofs = dofMap(:, elem_idx); 
            s = double(obj.Mesh.T2E_Sign(:, elem_idx));
            
            local_coeffs = solVec(dofs) .* s; 
            val = N_phy' * local_coeffs;
        end
        
        function grad_val = internalProbeLagrangeGrad(obj, solVec, elem_idx, spaceName)
            % internalProbeLagrangeGrad 封装 Lagrange P1 梯度计算逻辑
            dofMap = obj.DofHandler.DofMaps(spaceName);
            
            % 处理全局索引到局部向量索引的转换
            global_dofs = dofMap(:, elem_idx);
            if obj.DofHandler.SpaceOffsets.isKey(spaceName)
                offset = obj.DofHandler.SpaceOffsets(spaceName);
            else
                offset = 0;
            end
            local_dofs = global_dofs - offset;
            
            nodes = obj.Mesh.T(:, elem_idx); 
            p_elem = obj.Mesh.P(:, nodes);
            [J_mat, ~] = compute_jacobian_tet(p_elem);
            
            % P1 梯度参考值 (常数)
            q_ref = [0.25; 0.25; 0.25];
            [~, grad_ref] = lagrange_tet_p1(q_ref); % [4x3]
            
            % H(grad) Piola 变换: grad = J^{-T} * grad_ref
            invJt = inv(J_mat)';
            Grad_phy = grad_ref * invJt'; % [4x3]
            
            local_coeffs = solVec(local_dofs);
            grad_val = Grad_phy' * local_coeffs; % [3x1]
        end
        
        function val = getPropertyValue(~, mapOrArray, tag, defaultVal)
            % 统一处理 Map 或 Array 形式的属性查询
            % 允许用户传入容器 Map 或直接的数组作为材料/源属性
            if isa(mapOrArray, 'containers.Map')
                if mapOrArray.isKey(tag)
                    val = mapOrArray(tag);
                else
                    val = defaultVal;
                end
            elseif isnumeric(mapOrArray)
                if isscalar(mapOrArray)
                    val = mapOrArray; % 统一标量
                elseif numel(mapOrArray) == 3 && isvector(mapOrArray)
                     % 可能是 [Jx, Jy, Jz] 向量作为标量传入? (虽然有点歧义，但兼容旧逻辑)
                     val = mapOrArray;
                elseif tag <= length(mapOrArray)
                    val = mapOrArray(tag); % 按索引取值
                else
                    val = defaultVal;
                end
            else
                val = defaultVal;
            end
        end
        
        function total_W = integrate_energy_kernel(obj, PackedData, A_sol, MatLibData)
            % INTEGRATE_ENERGY_KERNEL 逐单元积分计算磁场能量 (并行加速)
            
            C_P = parallel.pool.Constant(PackedData.P);
            C_T = parallel.pool.Constant(PackedData.T);
            C_Dofs = parallel.pool.Constant(PackedData.CellDofs);
            C_Signs = parallel.pool.Constant(PackedData.Signs);
            C_Tags = parallel.pool.Constant(PackedData.RegionTags);
            C_Sol = parallel.pool.Constant(A_sol);
            C_MatLib = parallel.pool.Constant(MatLibData);
            
            % 获取积分点数据 (2阶精度足够)
            [q_pts, q_w] = get_quadrature_data('tet', 2);
            [~, curl_ref] = nedelec_tet_p1(q_pts);
            C_Qw = parallel.pool.Constant(q_w);
            C_CurlRef = parallel.pool.Constant(curl_ref);
            
            numElems = fix(full(double(obj.Mesh.NumElements)));
            n_q = length(q_w);
            
            W_local = zeros(numElems, 1);
            
            parfor e = 1:numElems
                loc_P = C_P.Value; loc_T = C_T.Value; loc_Sol = C_Sol.Value;
                loc_MatLib = C_MatLib.Value; loc_Tags = C_Tags.Value;
                loc_Signs = C_Signs.Value; loc_Dofs = C_Dofs.Value;
                loc_qw = C_Qw.Value; loc_curl_ref = C_CurlRef.Value;
                
                nodes = loc_T(:, e); p_elem = loc_P(:, nodes);
                dofs = loc_Dofs(:, e); s = loc_Signs(:, e);
                tag = loc_Tags(e); matData = loc_MatLib(tag);
                A_local = loc_Sol(dofs) .* double(s);
                
                [J_mat, detJ] = compute_jacobian_tet(p_elem); 
                invDetJ = 1.0 / detJ;
                
                w_sum = 0;
                for q = 1:n_q
                    weight = loc_qw(q) * abs(detJ); % 物理积分权重 = 权重 * 雅可比行列式
                    
                    % 计算积分点处的 B
                    C_ref = reshape(loc_curl_ref(:, :, q), 6, 3);
                    C_phy = (C_ref * J_mat') * invDetJ;
                    B_vec = C_phy' * A_local;
                    B_mag = norm(B_vec);
                    B_sq = B_mag^2;
                    
                    w_density = 0;
                    if strcmp(matData.Type, 'Linear')
                        % 线性能量密度: 0.5 * nu * B^2
                        w_density = 0.5 * matData.Nu_Linear * B_sq;
                    else
                        % 非线性能量密度 (共能): W = \int H dB = \int nu(B)*B dB
                        % 使用辛普森/梯形公式进行简单的数值积分近似
                        b_samples = [0, 0.5*B_mag, B_mag];
                        b_sq_samples = b_samples.^2;
                        nu_vals = zeros(1,3);
                        for k=1:3
                             [nu_k, ~] = MaterialLib.evaluate(b_sq_samples(k), matData);
                             nu_vals(k) = nu_k;
                        end
                        h_vals = nu_vals .* b_samples;
                        % Simpson / Trapezoidal rule
                        w_density = (B_mag / 6.0) * (0 + 4*h_vals(2) + h_vals(3));
                    end
                    w_sum = w_sum + weight * w_density;
                end
                W_local(e) = w_sum;
            end
            total_W = sum(W_local);
        end
    end
end