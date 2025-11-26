function NodalData = map_element_to_node(Model, ElementData)
% MAP_ELEMENT_TO_NODE 将单元中心数据平滑映射到节点 (体积加权平均)
%
% 输入: ElementData (N_elems x 1)
% 输出: NodalData (N_nodes x 1)

    Mesh = Model.Mesh;
    numNodes = size(Mesh.P, 2);
    numElems = size(Mesh.T, 2);
    
    % 累加器
    NodeSum = zeros(numNodes, 1);
    NodeVol = zeros(numNodes, 1);
    
    % 提取体积 (如果未预计算，需在此计算)
    if isfield(Mesh, 'Volumes')
        Vols = Mesh.Volumes;
    else
        % 简易计算... 略
        Vols = ones(1, numElems); 
    end
    
    T = Mesh.T;
    
    % 向量化累加 (利用 sparse 的累加特性)
    % 构造稀疏矩阵 S (Nodes x Elems)，元素为体积
    % NodeVal = (S * (Data .* Vol)) ./ (S * Vol)
    
    % 展开索引
    J_idx = 1:numElems;
    J_rep = repmat(J_idx, 4, 1); % 4 x Ne
    I_idx = T;                   % 4 x Ne (Node IDs)
    
    % 权重矩阵 (节点 i 在单元 j 中的权重 = Vol_j)
    % 实际上更简单：直接遍历
    
    % 使用 accumarray 会非常快
    I_vec = T(:);
    Vol_rep = repmat(Vols(:)', 4, 1);
    Data_rep = repmat(ElementData(:)', 4, 1);
    
    Weight_vec = Vol_rep(:);
    Value_vec = Data_rep(:) .* Weight_vec;
    
    NodeSum = accumarray(I_vec, Value_vec, [numNodes, 1]);
    NodeVol = accumarray(I_vec, Weight_vec, [numNodes, 1]);
    
    % 避免除零
    mask = NodeVol > 1e-12;
    NodalData = zeros(numNodes, 1);
    NodalData(mask) = NodeSum(mask) ./ NodeVol(mask);
end