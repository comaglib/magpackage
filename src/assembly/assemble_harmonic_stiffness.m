function K_global = assemble_harmonic_stiffness(Model, Nu_harmonics, Harmonics, Info_Nu)
% ASSEMBLE_HARMONIC_STIFFNESS 组装 HBFEM 的 Block-Toeplitz 刚度矩阵 (修正版)
%
% 输入:
%   Model        - 模型结构体
%   Nu_harmonics - (N_elems x N_nu_harm) 磁阻率的频域系数矩阵
%   Harmonics    - 求解的谐波次数列表 (A场)
%   Info_Nu      - [新增] 磁阻率谐波的索引信息，用于映射物理阶次到矩阵列

    numEdges = size(Model.Mesh.Edges, 2);
    numHarm = length(Harmonics);
    
    % 1. 建立查找表: 物理谐波阶数 (Physical Order) -> Nu矩阵列索引 (Column Index)
    % Info_Nu.Harmonics 包含如 [0, 2, 4, 6]
    Nu_Order_Map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:length(Info_Nu.Harmonics)
        Nu_Order_Map(Info_Nu.Harmonics(k)) = k;
    end
    
    K_blocks = cell(numHarm, numHarm);
    
    % 2. 填充 Block-Toeplitz 矩阵
    % T_{mn} = Nu_{m-n}
    for r = 1:numHarm
        m = Harmonics(r);
        for c = 1:numHarm
            n = Harmonics(c);
            
            diff_order = m - n;       % 差频
            abs_order = abs(diff_order);
            
            % 检查该阶次是否存在于计算出的 Nu 谐波中
            if isKey(Nu_Order_Map, abs_order)
                col_idx = Nu_Order_Map(abs_order);
                Nu_vec = Nu_harmonics(:, col_idx);
                
                % 共轭逻辑: 
                % Nu(t) 是实数，频域满足 Nu(-k) = conj(Nu(k))
                % 我们存储的是正频率部分。如果 diff_order < 0，取共轭。
                if diff_order < 0
                    Nu_vec = conj(Nu_vec);
                end
                
                % 组装子块
                K_sub = assemble_magnetic_stiffness(Model, Nu_vec);
                K_blocks{r, c} = K_sub;
            else
                % 如果差频阶次太高未被覆盖 (例如 m=5, n=-5 -> 10次谐波)，
                % 且 Info_Nu 没算到 10 次，则忽略 (视为0)
                % 这种情况通常发生在 Info_Nu 截断过早时。
                K_blocks{r, c} = sparse(numEdges, numEdges);
            end
        end
    end
    
    K_global = cell2mat(K_blocks);
end