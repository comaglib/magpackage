function ElemA = prepare_element_A(Model, x)
% PREPARE_ELEMENT_A 将全局解向量映射为单元局部向量 (含符号修正)
%
% 输入:
%   Model - 模型结构体 (包含 T2E 和 T2E_Sign)
%   x     - 全局解向量 (N_edges x 1)
%
% 输出:
%   ElemA - 单元局部系数矩阵 (6 x N_elems)
%           每一列对应一个单元的 [a1; a2; ...; a6]
%           已乘入符号修正，可直接用于 calc_element_B

    % fprintf('  - 正在映射局部解向量 (Gathering A)...\n');
    
    T2E = Model.Mesh.T2E;       % 6 x Ne
    Signs = double(Model.Mesh.T2E_Sign); % 6 x Ne
    
    % 1. 向量化提取
    % x(T2E) 会生成与 T2E 同维度的矩阵 (6 x Ne)
    % 这里的 x 必须是列向量
    if size(x, 2) > 1, x = x(:); end
    
    RawA = x(T2E);
    
    % 2. 符号修正
    % Local_A = Global_A * Sign
    ElemA = RawA .* Signs;
    
end