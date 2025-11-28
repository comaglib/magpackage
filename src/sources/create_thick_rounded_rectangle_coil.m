function Coil = create_thick_rounded_rectangle_coil(CoilParams, N_sub)
% CREATE_THICK_ROUNDED_RECTANGLE_COIL 生成有厚度的多丝线圈模型
%
% 原理: 在矩形截面上生成 NxN 的点阵，将单根中心路径复制并平移到这些点上。
%       这能有效解决近场计算中的奇异性问题。

    if nargin < 2, N_sub = 5; end % 默认 5x5 = 25 根丝
    
    % 基础参数
    I_total = CoilParams.Current;
    I_strand = I_total / (N_sub * N_sub); % 电流均分
    
    W = CoilParams.Width;
    H = CoilParams.Height;
    
    % 生成局部截面网格 (均匀分布的中心点)
    du = W / N_sub; 
    dv = H / N_sub;
    u_vec = linspace(-W/2 + du/2, W/2 - du/2, N_sub);
    v_vec = linspace(-H/2 + dv/2, H/2 - dv/2, N_sub);
    
    [U_grid, V_grid] = ndgrid(u_vec, v_vec);
    offsets_u = U_grid(:);
    offsets_v = V_grid(:);
    n_strands = length(offsets_u);
    
    % 容器
    All_P1 = [];
    All_P2 = [];
    All_I  = [];
    
    % 循环生成每根丝
    for k = 1:n_strands
        d_radial = offsets_u(k);
        d_z      = offsets_v(k);
        
        % 构建该丝的参数
        StrandParams = CoilParams;
        StrandParams.Current = I_strand;
        
        % 修改几何参数 (径向偏移影响 Lx, Ly 和 R)
        StrandParams.Lx_mean = CoilParams.Lx_mean + 2 * d_radial;
        StrandParams.Ly_mean = CoilParams.Ly_mean + 2 * d_radial;
        StrandParams.R_mean  = CoilParams.R_mean  + d_radial;
        
        % Z向偏移直接加到 Center 上
        StrandParams.Center(3) = CoilParams.Center(3) + d_z;
        
        % 保护 R (防止内缩过头)
        if StrandParams.R_mean < 0, StrandParams.R_mean = 0; end
        
        % 调用单丝生成器 (依赖 create_rounded_rectangle_coil.m)
        Strand = create_rounded_rectangle_coil(StrandParams);
        
        % 收集
        All_P1 = [All_P1, Strand.P1]; %#ok<AGROW>
        All_P2 = [All_P2, Strand.P2]; %#ok<AGROW>
        All_I  = [All_I,  Strand.I];  %#ok<AGROW>
    end
    
    % 合并结果
    Coil.P1 = All_P1;
    Coil.P2 = All_P2;
    Coil.I  = All_I;
    Coil.Turns = 1; 
end