function [K, b] = apply_dirichlet_bc(K, b, dof_indices, values)
% APPLY_DIRICHLET_BC 对稀疏矩阵施加 Dirichlet 边界条件
% 采用保持对称性的置 1 置 0 法。
%
% 输入:
%   K           - 全局刚度矩阵 (N x N)
%   b           - 右端项向量 (N x 1)
%   dof_indices - 需要固定的自由度索引 (向量)
%   values      - 固定的数值 (标量 0，或与 dof_indices 等长的向量)
%
% 输出:
%   K, b        - 修改后的矩阵和向量

    if isempty(dof_indices)
        return;
    end
    
    fprintf('  - 正在施加 Dirichlet 边界条件 (Fixed DoFs: %d)...\n', length(dof_indices));
    
    N = size(K, 1);
    
    % 统一 values 格式
    if isscalar(values)
        values = repmat(values, length(dof_indices), 1);
    end
    
    % --------------------------------------------------------
    % 步骤 1: 处理非零边界值对 RHS 的影响 (Move to RHS)
    % b = b - A(:, boundary) * values
    % 仅当 values 不全为 0 时需要计算
    % --------------------------------------------------------
    if any(values ~= 0)
        % 提取受约束的列
        K_constrained = K(:, dof_indices);
        % 更新 RHS
        b = b - K_constrained * values;
        
        % 修正被约束行自身的 RHS (将在下面强制覆盖，这里先不管)
    end
    
    % --------------------------------------------------------
    % 步骤 2: 矩阵置零 (Zero rows and columns)
    % --------------------------------------------------------
    % 策略: 利用稀疏矩阵索引操作
    % K(dof, :) = 0; K(:, dof) = 0;
    % 但直接操作稀疏矩阵索引很慢。
    % 高效方法: 构建一个保留掩码 (Keep Mask) 对角阵
    
    % 找出所有不受约束的自由度
    free_dofs = true(N, 1);
    free_dofs(dof_indices) = false;
    
    % 构建对角掩码矩阵 I_mask (1 for free, 0 for fixed)
    % I_mask = spdiags(double(free_dofs), 0, N, N);
    
    % K_mod = I_mask * K * I_mask; 
    % 这会将行和列都置零。
    % 然后再把对角线设为 1。
    
    % 使用逻辑索引直接置零（MATLAB 近版本优化过，速度尚可）
    % 为了保持对称性，行和列都要置零
    % 注意：这会改变稀疏结构，但在 MUMPS 中通常没问题
    
    % 更快的方法：只修改对角线为 1，非对角线为 0
    % 我们可以遍历 K 的三元组来重建，或者直接置零
    
    % 为了简便且保持对称性：
    K(:, dof_indices) = 0; % 列置零
    K(dof_indices, :) = 0; % 行置零
    
    % --------------------------------------------------------
    % 步骤 3: 对角线置 1
    % --------------------------------------------------------
    new_diagonals = sparse(dof_indices, dof_indices, 1, N, N);
    K = K + new_diagonals;
    
    % --------------------------------------------------------
    % 步骤 4: 强制设定 RHS
    % --------------------------------------------------------
    b(dof_indices) = values;
    
end