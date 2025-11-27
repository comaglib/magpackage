function [Solution, Info] = solve_hbfem_voltage(Model, Coil, Circuit, HBFEMParams)
% SOLVE_HBFEM_VOLTAGE 电压驱动型谐波平衡有限元求解器 (HBFEM)
%
% =========================================================================
% 理论基础 (Theoretical Basis):
%   求解非线性时谐涡流方程 (A-Formulation):
%   curl( nu(B) * curl(A) ) + sigma * dA/dt = Js
%
%   采用分解式谐波平衡法 (Decomposed HBFEM / Fixed-Point):
%   将非线性磁阻率 nu(B) 分解为 直流分量 nu_dc 和 扰动分量 nu_pert(t)。
%   方程被线性化为一系列解耦的谐波方程:
%   
%   [K(nu_dc) + j*w_k*M_sigma] * A_k = J_source_k - curl( M_resid_k )
%
%   其中 M_resid (残差磁化强度) 等效于非线性修正源:
%   M_resid(t) = ( nu(t) - nu_dc ) * B(t)
%
% 核心特性 (Key Features):
%   1. **连续源步进 (Continuous Source Stepping)**: 
%      融合了"负载推进"与"非线性迭代"。在迭代初期电压从 0->100% 线性增加，
%      引导解始终沿着物理路径演化，极大提高了强饱和问题的收敛性。
%
%   2. **物理修正磁化源 (Corrected Magnetization)**: 
%      严格推导的源项公式 M_src = (nu0 - nu_dc)*Bs，保证了源场在介质中的
%      正确分布，解决了开路磁芯 B 场过低的问题。
%
%   3. **快速边界处理 (Sparse Optimization)**:
%      使用稀疏矩阵掩码运算替代传统的行列索引操作，大网格下速度提升显著。
%
% =========================================================================

    fprintf('==============================================\n');
    fprintf('   HBFEM Solver (Field-Circuit Coupled)       \n');
    fprintf('   Method: Fixed-Point + Continuous Ramping   \n');
    fprintf('==============================================\n');

    %% 1. 初始化与频域参数设置
    base_freq = HBFEMParams.Frequency;
    harmonics = HBFEMParams.Harmonics; 
    numHarm = length(harmonics);
    
    % 电压向量维度对齐
    if length(Circuit.Voltage) < numHarm
        Circuit.Voltage = [Circuit.Voltage(:); zeros(numHarm-length(Circuit.Voltage), 1)];
    end
    
    % FFT 时间采样点数 (满足 Shannon 采样定理: N > 2*MaxHarm)
    max_h = max(harmonics);
    min_req = 2 * max_h + 1;
    n_time = max(32, 2^nextpow2(min_req + 2));
    if isfield(HBFEMParams, 'TimeSteps') && ~isempty(HBFEMParams.TimeSteps) && HBFEMParams.TimeSteps >= min_req
        n_time = HBFEMParams.TimeSteps;
    end
    fprintf('  [Setup] FFT Time Samples: %d\n', n_time);
    
    % 自由度映射 (Edge Elements)
    DofData = create_dof_map(Model);
    numEdges = DofData.NumEdges;
    
    % 启动并行池 (用于材料属性更新)
    if isempty(gcp('nocreate')), try parpool('local'); catch; end; end
    
    % AFT 模块 (Alternating Frequency/Time)
    AFT = aft_module();
    Info_A = AFT.prepare_indices(harmonics, n_time);
    
    %% 2. 组装时不变矩阵 (Time-Invariant Matrices)
    % 涡流质量矩阵: M_ij = int( sigma * Ni * Nj ) dV
    M_sigma = assemble_mass_matrix(Model, 'ConductorOnly');
    % 正则化质量矩阵 (用于稳定空气域零频解)
    M_reg   = assemble_mass_matrix(Model, 'All');
    
    if isfield(Coil, 'Turns'), N_turns = Coil.Turns; else, N_turns = 1; end
    fprintf('  [Init] Coil Turns: %d\n', N_turns);

    % 绕组耦合向量 W (几何因子)
    % 物理意义: 磁链 Psi = W' * A. 
    % 注意: 必须乘以 N_turns，因为 FEM 几何仅建模了单匝分布。
    fprintf('  [Init] Assembling Winding Vector...\n');
    W_vec_geom = assemble_winding_vector(Model, Coil);
    W_vec = W_vec_geom * N_turns; 
    
    % 预计算源磁场 (Biot-Savart Law)
    % 用于计算开路磁芯的源激励项 Bs
    fprintf('  [Init] Pre-computing Source Field (Biot-Savart)...\n');
    Mesh = Model.Mesh;
    numElems = size(Mesh.T, 2);
    Centers = zeros(3, numElems);
    P = Mesh.P; T = Mesh.T;
    for i = 1:numElems, Centers(:,i) = sum(P(:, T(:,i)), 2) / 4; end
    
    CoilUnit = Coil; 
    CoilUnit.I(:) = 1.0 * N_turns; % 单位电流对应 N 安匝
    
    Bs_unit_raw = compute_biot_savart_B_serial(CoilUnit, Centers);
    Bs_unit_global = Bs_unit_raw;
    As_unit_vec = project_source_A_on_edges(Model, CoilUnit);
    
    % 源电感 (漏感+空心电感)
    L_source = full(W_vec' * As_unit_vec);
    
    % 边界条件索引
    if isfield(Model.Runtime,'FixedEdges'), fe=Model.Runtime.FixedEdges; else, fe=[]; end
    
    %% 3. 求解状态初始化
    X_curr = zeros(numEdges, numHarm); % 磁矢位谐波系数 [A_1, A_3, ...]
    
    % 初始磁阻率 Nu_DC (线性启动)
    % 使用材料定义的初始 Mu_r，避免从真空启动导致的数值震荡
    Nu_DC = zeros(numElems, 1);
    mu0 = 4*pi*1e-7;
    nu0 = 1/mu0; 
    MatMap = Model.Materials.ActiveMap;
    MatLib = Model.Materials.Lib;
    
    for i = 1:numElems
        mat_id = MatMap(i);
        if isfield(MatLib(mat_id), 'Mu_r') && ~isempty(MatLib(mat_id).Mu_r)
            Nu_DC(i) = 1.0 / (MatLib(mat_id).Mu_r * mu0);
        else
            Nu_DC(i) = nu0;
        end
    end
    
    % 初始电流猜测 (基于线性阻抗)
    I_curr = zeros(numHarm, 1);
    for k = 1:numHarm
        wk = harmonics(k) * 2 * pi * base_freq;
        Zk = Circuit.R + 1i * wk * (Circuit.L_leak + L_source * 100); 
        if abs(Zk)<1e-9, Zk=1e-3; end
        I_curr(k) = Circuit.Voltage(k) / Zk;
    end
    
    % 求解器配置
    Model.Solver.Linear.Interface = 'MUMPS';
    Model.Solver.Linear.Symmetric = true;
    use_mumps = exist('linear_solve', 'file');
    
    fprintf('  [HBFEM] Starting Iterations (Ramping 0->100%%)...\n');
    
    %% 4. 非线性迭代主循环 (Continuous Ramping)
    % 策略: 
    %   Iter 1 ~ RampSteps: 电压线性增加，同时更新 Nu。
    %   Iter > RampSteps: 电压维持 100%，直至收敛。
    
    MaxIter = 100;      
    RampSteps = 25;     % 负载爬坡步数
    Tol = 1e-4;         % 收敛容差
    
    for iter = 1:MaxIter
        % 4.1 确定当前负载比例
        if iter <= RampSteps
            load_scale = iter / RampSteps;
        else
            load_scale = 1.0;
        end
        V_target = Circuit.Voltage * load_scale;
        
        % 4.2 组装基准刚度矩阵 (LHS)
        % K(nu_dc) = int( curlN * nu_dc * curlN ) dV
        % 这是定点迭代的左端算子，保持正定性。
        K_base = assemble_magnetic_stiffness(Model, Nu_DC);
        
        % 计算正则化因子
        scale_K = mean(abs(nonzeros(K_base)));
        eps_val = 1e-6 * scale_K;
        
        % 4.3 更新材料属性 & 计算非线性残差 (RHS)
        % 核心步骤: A -> B(t) -> H(t) -> nu(t) -> FFT -> M_resid
        [Nu_DC_sug, M_resid_harm, M_src_unit_harm, max_B, mean_B_iron] = ...
            update_properties_voltage_mode(Model, X_curr, AFT, Info_A, I_curr, Nu_DC, Bs_unit_global, nu0);
        
        % 4.4 准备线性源项驱动力
        % M_lin = (nu0 - nu_dc) * Bs. 
        % 物理意义: 源电流产生的真空场 Bs 在介质中产生的等效磁化强度。
        M_lin_vec = Bs_unit_global .* repmat((nu0 - Nu_DC)', 3, 1);
        RHS_src_linear = assemble_rhs_from_magnetization(Model, M_lin_vec);
        
        X_new = zeros(size(X_curr));
        I_new = zeros(size(I_curr));
        
        % 4.5 逐谐波求解 (Parallel/Loop)
        for k = 1:numHarm
            h_order = harmonics(k);
            wk = h_order * 2 * pi * base_freq;
            
            % 构建线性系统矩阵 [K + jwM]
            if wk == 0
                K_sys = K_base + eps_val * M_reg;
            else
                K_sys = K_base + 1i * wk * M_sigma + eps_val * M_reg;
            end
            
            % --- 组装对偶右端项 (Dual-RHS) ---
            
            % RHS 1: 非线性残差驱动项 (Nonlinear Residual Force)
            % M_res = (nu(t) - nu_dc) * B_tot
            % 这一项补偿了 LHS 使用常数 nu_dc 带来的误差。
            M_vec_res = reshape(squeeze(M_resid_harm(:, k, :)), 3, numElems);
            RHS_resid = - assemble_rhs_from_magnetization(Model, M_vec_res);
            
            % RHS 2: 单位电流驱动项 (Unit Current Response)
            % 包含: 线性源场驱动 + 涡流反作用 + 源场非线性修正
            RHS_eddy = -1i * wk * M_sigma * As_unit_vec;
            
            M_vec_src_harm = reshape(squeeze(M_src_unit_harm(:, k, :)), 3, numElems);
            RHS_src_harm = assemble_rhs_from_magnetization(Model, M_vec_src_harm);
            
            if k == 1
                 % 基波包含线性源项
                 RHS_unit = RHS_src_linear + RHS_src_harm + RHS_eddy;
            else
                 RHS_unit = RHS_src_harm + RHS_eddy;
            end
            
            % 快速施加边界条件
            RHS_dual = [RHS_resid, RHS_unit];
            [K_bc, R_dual_bc] = apply_dirichlet_bc_multi(K_sys, RHS_dual, fe, zeros(length(fe), 1));
            
            % 线性方程组求解
            if use_mumps
                try
                    X_dual = linear_solve(K_bc, R_dual_bc, Model);
                catch
                    X_dual = K_bc \ R_dual_bc;
                end
            else
                X_dual = K_bc \ R_dual_bc;
            end
            
            % 分离解向量: A = A_res + I * A_unit
            A_resid = X_dual(:, 1);
            A_unit  = X_dual(:, 2);
            
            % 计算磁链耦合
            Psi_resid = W_vec' * A_resid;
            Psi_unit  = W_vec' * A_unit;
            
            % --- 4.6 求解电路耦合方程 ---
            % 电路方程: V = (R + jw*L_leak)*I + jw * Psi_total
            % Psi_total = Psi_resid + I * (Psi_unit + L_source)
            
            Z_leak = Circuit.R + 1i * wk * Circuit.L_leak;
            Z_total = Z_leak + 1i * wk * (Psi_unit + L_source); 
            
            if abs(Z_total) < 1e-6, Z_total = 1e-3; end
            
            % 计算电流
            I_calc = (V_target(k) - 1i * wk * Psi_resid) / Z_total;
            I_new(k) = I_calc;
            
            % 合成最终磁矢位
            X_new(:, k) = A_resid + I_new(k) * A_unit;
        end
        
        % 4.7 收敛性检查
        err_A = norm(X_new - X_curr, 'fro') / (norm(X_new, 'fro') + 1e-10);
        err_I = norm(I_new - I_curr) / (norm(I_new) + 1e-10);
        err = max(err_A, err_I);
        
        fprintf('    Step %d: Load=%.0f%%, Max|B|=%.2fT, I_fund=%.1fA, Err=%.2e\n', ...
            iter, load_scale*100, max_B, abs(I_new(1)), err);
        
        % 仅在满负载且误差满足要求时退出
        if load_scale >= 1.0 && err < Tol
            fprintf('  [Converged] Solution stabilized.\n');
            X_curr = X_new; I_curr = I_new; Nu_DC = Nu_DC_sug;
            break;
        end
        
        % 4.8 状态更新与松弛 (Relaxation)
        alpha = 0.5;
        if iter > 5
            if err < 0.1, alpha = 0.8; end
            if err < 0.01, alpha = 1.0; end
        end
        
        X_curr = (1-alpha)*X_curr + alpha*X_new;
        I_curr = (1-alpha)*I_curr + alpha*I_new;
        Nu_DC  = (1-alpha)*Nu_DC  + alpha*Nu_DC_sug;
    end
    
    Solution.Harmonics = X_curr; 
    Solution.Current = I_curr;
    Info.Iterations = iter;
    Info.TimeSteps = n_time;
end

% -------------------------------------------------------------------------
% 辅助函数
% -------------------------------------------------------------------------

function [Nu_DC_sug, M_resid_harm, M_src_unit_harm, max_B, mean_B_iron] = update_properties_voltage_mode(Model, X_curr, AFT, Info_A, I_vec, Nu_DC_In, Bs_unit_global, nu0)
% UPDATE_PROPERTIES_VOLTAGE_MODE 更新非线性材料属性
%
% 输入:
%   X_curr, I_vec: 当前迭代的场与电流谐波
%   Nu_DC_In: 上一步的直流磁阻率
%
% 输出:
%   Nu_DC_sug: 建议的新直流磁阻率 (基于能量等效或时间平均)
%   M_resid_harm: 非线性残差磁化强度的谐波系数
%   M_src_unit_harm: 源场非线性修正项

    Mesh = Model.Mesh; numElems = size(Mesh.T, 2); numHarmA = Info_A.NumHarmonics; Nt = Info_A.Nt; 
    
    % IFFT: 重构时域电流波形
    I_t = zeros(1, Nt);
    for k = 1:numHarmA, c=zeros(1,numHarmA); c(k)=I_vec(k); I_t=I_t+AFT.freq2time(c, Info_A); end
    
    % 并行计算分块
    chunkSize = 10000; numChunks = ceil(numElems/chunkSize); ElementBlocks = cell(numChunks,1);
    T_raw=Mesh.T; T2E_raw=Mesh.T2E; Signs_raw=double(Mesh.T2E_Sign); MatMap_raw=Model.Materials.ActiveMap;
    
    for k=1:numChunks
        idx=(k-1)*chunkSize+1 : min(k*chunkSize,numElems);
        ElementBlocks{k}.T=T_raw(:,idx); ElementBlocks{k}.T2E=T2E_raw(:,idx);
        ElementBlocks{k}.Signs=Signs_raw(:,idx); ElementBlocks{k}.MatIDs=MatMap_raw(idx);
        ElementBlocks{k}.Count=length(idx); ElementBlocks{k}.Nu_DC_In=Nu_DC_In(idx); ElementBlocks{k}.Bs_Unit=Bs_unit_global(:,idx);
    end
    
    C_P=parallel.pool.Constant(Mesh.P); C_A=parallel.pool.Constant(X_curr);
    C_MatLib=parallel.pool.Constant(Model.Materials.Lib); C_AFT=parallel.pool.Constant(AFT); C_It=parallel.pool.Constant(I_t);
    
    Nu_DC_Out=cell(numChunks,1); M_resid_Out=cell(numChunks,1); M_src_Out=cell(numChunks,1); 
    B_max_Out=cell(numChunks,1); B_iron_Out=cell(numChunks,1); Count_iron_Out=cell(numChunks,1);
    
    parfor k=1:numChunks
        Block=ElementBlocks{k}; cnt=Block.Count;
        loc_P=C_P.Value; loc_A=C_A.Value; loc_MatLib=C_MatLib.Value; loc_It=C_It.Value; loc_AFT=C_AFT.Value;
        loc_Nu_In=Block.Nu_DC_In; loc_Bs_Unit=Block.Bs_Unit;
        
        nu_loc=zeros(cnt,1); m_res=zeros(3,numHarmA,cnt); m_src_unit=zeros(3,numHarmA,cnt); b_max=0;
        b_sum_iron=0; cnt_iron=0;
        G_ref=[-1 1 0 0; -1 0 1 0; -1 0 0 1];
        
        for i=1:cnt
            % 1. 几何雅可比计算 (Gradient Matrix)
            nodes=Block.T(:,i); p_vec=loc_P(:,nodes);
            v21=p_vec(:,2)-p_vec(:,1); v31=p_vec(:,3)-p_vec(:,1); v41=p_vec(:,4)-p_vec(:,1);
            J=[v21,v31,v41]; detJ=det(J); c34=cross(v31,v41); c42=cross(v41,v21); c23=cross(v21,v31);
            invJ_T=[c34,c42,c23]/detJ; G_phy=invJ_T*G_ref;
            
            % 2. 反应场 Br(t) 重构 (IFFT)
            edges=Block.T2E(:,i); s=Block.Signs(:,i); A_elem=loc_A(edges,:).*s; 
            A_t=loc_AFT.freq2time(A_elem,Info_A);
            
            Br_t=zeros(3,Nt); pairs=[1 2;1 3;1 4;2 3;2 4;3 4];
            for e=1:6, c=2*cross(G_phy(:,pairs(e,1)),G_phy(:,pairs(e,2))); Br_t=Br_t+c*A_t(e,:); end
            
            % 3. 总场合成 B_tot = Br + Bs
            Bs_unit_vec=loc_Bs_Unit(:,i); B_tot=Br_t+Bs_unit_vec*loc_It;
            B_mag_sq=sum(B_tot.^2,1); 
            if max(sqrt(B_mag_sq))>b_max, b_max=max(sqrt(B_mag_sq)); end
            
            mat_id=Block.MatIDs(i); 
            if strcmp(loc_MatLib(mat_id).Type, 'Nonlinear')
                b_sum_iron = b_sum_iron + mean(sqrt(B_mag_sq)); cnt_iron = cnt_iron + 1;
            end
            
            % 4. 查表更新磁阻率
            B_clamp_sq=min(B_mag_sq, 5.0^2); 
            [nu_t,~]=eval_material_nu(B_clamp_sq,loc_MatLib(mat_id));
            nu_loc(i)=mean(nu_t);
            
            % 5. 计算等效磁化源项 (Equivalent Magnetization)
            % M_res: 定点迭代残差 = (nu(t) - nu_dc) * B_tot
            M_resid_t_pure=(nu_t-loc_Nu_In(i)).*B_tot;
            
            % M_src: 源场修正 = (nu0 - nu(t)) * Bs - (nu0 - nu_dc) * Bs = (nu_dc - nu(t)) * Bs
            M_src_t = (loc_Nu_In(i) - nu_t) .* (Bs_unit_vec*loc_It);
            
            for d=1:3
                m_res(d,:,i)=loc_AFT.time2freq(M_resid_t_pure(d,:),Info_A);
                m_src_unit(d,:,i)=loc_AFT.time2freq(M_src_t(d,:),Info_A);
            end
        end
        Nu_DC_Out{k}=nu_loc; M_resid_Out{k}=m_res; M_src_Out{k}=m_src_unit; 
        B_max_Out{k}=b_max; B_iron_Out{k}=b_sum_iron; Count_iron_Out{k}=cnt_iron;
    end
    Nu_DC_sug=vertcat(Nu_DC_Out{:}); M_resid_harm=cat(3,M_resid_Out{:}); M_src_unit_harm=cat(3,M_src_Out{:}); 
    max_B=max([B_max_Out{:}]);
    
    total_iron_B = sum([B_iron_Out{:}]);
    total_iron_cnt = sum([Count_iron_Out{:}]);
    if total_iron_cnt > 0, mean_B_iron = total_iron_B / total_iron_cnt; else, mean_B_iron = 0; end
end

function [K, B] = apply_dirichlet_bc_multi(K, B, dofs, vals)
% APPLY_DIRICHLET_BC_MULTI 快速边界条件应用 (Sparse Masking)
% 
% 优化说明:
%   传统的 K(:,dofs)=0 索引操作在 MATLAB 中对稀疏矩阵非常慢。
%   此处采用稀疏对角矩阵乘法: K_new = D * K * D + I_fixed
%   其中 D 为掩码矩阵 (Free=1, Fixed=0)。这种全矩阵运算经过高度优化。

    if isempty(dofs), return; end
    
    n = size(K, 1);
    
    % 1. 构造掩码对角阵 D
    mask = true(n, 1);
    mask(dofs) = false;
    D = spdiags(double(mask), 0, n, n);
    
    % 2. 矩阵行列置零 (Sparse Multiplication)
    K = D * K * D;
    
    % 3. 对角线置 1
    I_fixed = spdiags(double(~mask), 0, n, n);
    K = K + I_fixed;
    
    % 4. 修正 RHS
    if size(vals,1)~=length(dofs), vals=repmat(vals(1),length(dofs),1); end
    if all(vals == 0)
        B(dofs, :) = 0;
    else
        B(dofs, :) = repmat(vals, 1, size(B, 2));
    end
end

function RHS_vec = assemble_rhs_from_magnetization(Model, M_vec)
% ASSEMBLE_RHS_FROM_MAGNETIZATION 组装磁化电流产生的源向量
% 公式: b_i = int( curl(Ni) . M ) dV
    Mesh=Model.Mesh; numEdges=size(Mesh.Edges,2); numElems=size(Mesh.T,2);
    I_idx=zeros(6*numElems,1); V_val=zeros(6*numElems,1); G_ref=[-1 1 0 0;-1 0 1 0;-1 0 0 1];
    P=Mesh.P; T=Mesh.T; T2E=Mesh.T2E; Signs=double(Mesh.T2E_Sign);
    for e=1:numElems
        nodes=T(:,e); pts=P(:,nodes); v21=pts(:,2)-pts(:,1); v31=pts(:,3)-pts(:,1); v41=pts(:,4)-pts(:,1);
        J=[v21,v31,v41]; detJ=det(J); invJ_T=inv(J)'; G_phy=invJ_T*G_ref; Vol=abs(detJ)/6.0; M_e=M_vec(:,e);
        for i=1:6
            na=[1 1 1 2 2 3]; nb=[2 3 4 3 4 4]; curl_ni=2*cross(G_phy(:,na(i)),G_phy(:,nb(i)));
            val=Vol*dot(curl_ni,M_e)*Signs(i,e); idx_store=(e-1)*6+i; I_idx(idx_store)=T2E(i,e); V_val(idx_store)=val;
        end
    end
    RHS_vec=sparse(I_idx,1,V_val,numEdges,1);
end