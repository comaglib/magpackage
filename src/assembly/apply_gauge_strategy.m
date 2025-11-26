function [SysK, Sysb, Info] = apply_gauge_strategy(Model, K, b, fixed_edges, fixed_vals, method)
% APPLY_GAUGE_STRATEGY 应用库伦规范并构建最终线性系统
%
% 自动化处理:
% 1. Lagrange: 组装耦合矩阵 C，构建鞍点系统，自动查找边界节点并固定 Lambda=0
% 2. Penalty:  组装质量矩阵 M，自动计算 Adaptive Alpha，构建正则化系统
%
% 输入:
%   Model       - 模型结构体 (用于组装 M/C 和查找拓扑)
%   K           - 原始刚度矩阵 (N_edge x N_edge)
%   b           - 原始右端项
%   fixed_edges - 边界棱索引 (Dirichlet BC)
%   fixed_vals  - 边界棱数值 (通常为标量或向量)
%   method      - 'Lagrange' (默认) 或 'Penalty'
%
% 输出:
%   SysK, Sysb  - 可以直接送入 linear_solve 的最终系统
%   Info        - 包含 Alpha 值、额外 DoF 信息等诊断数据

    if nargin < 6, method = 'Lagrange'; end
    
    numEdges = size(K, 1);
    
    % 统一 fixed_vals 格式
    if isscalar(fixed_vals)
        fixed_vals = repmat(fixed_vals, length(fixed_edges), 1);
    end

    switch lower(method)
        case 'lagrange'
            fprintf('[Gauge] Applying Lagrange Multiplier Strategy...\n');
            
            % 1. 组装耦合矩阵
            C = assemble_coupling_matrix(Model);
            numNodes = size(C, 2);
            
            % 2. 构建鞍点矩阵 [K C; C' 0]
            % 使用 sparse 拼接以优化内存
            SysK = [K, C; C', sparse(numNodes, numNodes)];
            
            % 3. 扩展右端项
            Sysb = [b; zeros(numNodes, 1)];
            
            % 4. 自动推导 Lambda 的边界条件
            % 物理逻辑: 在 A 被固定的边界上，虚构电势 Lambda 必须接地 (0)
            % 从边界棱中提取关联的节点
            edge_defs = Model.Mesh.Edges(:, fixed_edges); % 2 x N_fixed
            boundary_nodes = unique(edge_defs(:));
            
            fprintf('        Identified %d boundary nodes for Lambda pinning.\n', length(boundary_nodes));
            
            % 5. 施加所有边界条件
            % 列表A: 原始棱 BC
            dofs_A = fixed_edges;
            vals_A = fixed_vals;
            
            % 列表B: Lambda BC (偏移 numEdges)
            dofs_L = numEdges + boundary_nodes;
            vals_L = zeros(size(boundary_nodes));
            
            % 合并施加
            all_dofs = [dofs_A; dofs_L];
            all_vals = [vals_A; vals_L];
            
            [SysK, Sysb] = apply_dirichlet_bc(SysK, Sysb, all_dofs, all_vals);
            
            % 输出信息
            Info.Method = 'Lagrange';
            Info.NumEdges = numEdges;
            Info.NumNodes = numNodes;
            
        case 'penalty'
            fprintf('[Gauge] Applying Adaptive Penalty Strategy...\n');
            
            % 1. 组装质量矩阵
            M = assemble_mass_matrix(Model, 'All');
            
            % 2. 计算自适应 Alpha
            % 基准: 刚度矩阵平均强度
            ref_scale = mean(abs(nonzeros(K)));
            epsilon = 1e-6; % 经验系数: 既能压制零空间，又不影响物理精度
            alpha = ref_scale * epsilon;
            
            fprintf('        Adaptive Alpha = %.2e (Ref=%.2e)\n', alpha, ref_scale);
            
            % 3. 构建正则化矩阵
            SysK = K + alpha * M;
            Sysb = b; % RHS 不变
            
            % 4. 施加边界条件
            [SysK, Sysb] = apply_dirichlet_bc(SysK, Sysb, fixed_edges, fixed_vals);
            
            % 输出信息
            Info.Method = 'Penalty';
            Info.Alpha = alpha;
            
        otherwise
            error('Unknown gauge method: %s', method);
    end
end