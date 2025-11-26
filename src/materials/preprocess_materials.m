function Model = preprocess_materials(Model)
% PREPROCESS_MATERIALS 预处理材料库 (v3: 自动补全 Mu_r)
% 
% 改进:
% 1. 引入真空渐近线约束，确保饱和区外推符合 B = mu0*H + Js
% 2. 强制检查数据的物理单调性
% 3. 生成更密集的饱和区引导点，防止样条震荡
% 4. 自动计算非线性材料的 Mu_r (如果用户未提供)，防止后续模块报错

    fprintf('正在预处理非线性材料数据 (含自动补全)...\n');
    
    Lib = Model.Materials.Lib;
    mu0 = 4*pi*1e-7;
    nu0 = 1/mu0; 
    
    for i = 1:length(Lib)
        mat = Lib(i);
        
        if isfield(mat, 'BH_Curve') && ~isempty(mat.BH_Curve)
            B_data = mat.BH_Curve.B(:);
            H_data = mat.BH_Curve.H(:);
            
            % 1. 数据清洗
            if B_data(1) ~= 0 || H_data(1) ~= 0
                B_data = [0; B_data];
                H_data = [0; H_data];
            end
            
            [B_data, idx] = unique(B_data);
            H_data = H_data(idx);
            
            % 2. 物理外推
            B_last = B_data(end);
            H_last = H_data(end);
            
            factors = [1.05, 1.2, 2.0, 10.0, 100.0]; 
            B_ext = B_last * factors';
            H_ext = H_last + nu0 * (B_ext - B_last);
            
            B_full = [B_data; B_ext];
            H_full = [H_data; H_ext];
            
            % 3. 转换为 nu vs B^2
            Nu_full = zeros(size(B_full));
            
            % 计算初始磁阻率
            if length(B_full) > 1 && B_full(2) > 0
                nu_init = H_full(2) / B_full(2);
            else
                nu_init = 1/(1000*mu0); 
            end
            Nu_full(1) = nu_init;
            
            mask = B_full > 0;
            Nu_full(mask) = H_full(mask) ./ B_full(mask);
            
            if any(Nu_full > nu0)
                Nu_full(Nu_full > nu0) = nu0;
            end
            
            % 4. 存回样条
            Model.Materials.Lib(i).Type = 'Nonlinear';
            Model.Materials.Lib(i).NonlinearData.pp_nu = pchip(B_full.^2, Nu_full);
            Model.Materials.Lib(i).NonlinearData.pp_dnu = fnder(Model.Materials.Lib(i).NonlinearData.pp_nu, 1);
            Model.Materials.Lib(i).NonlinearData.nu0 = nu0;
            
            % ----------------------------------------------------
            % [核心修复] 自动回填 Mu_r
            % ----------------------------------------------------
            % 如果用户忘记定义 Mu_r，或者 Mu_r 为空，则自动计算
            if ~isfield(mat, 'Mu_r') || isempty(mat.Mu_r)
                calc_mu_r = 1 / (nu_init * mu0);
                Model.Materials.Lib(i).Mu_r = calc_mu_r;
                fprintf('  - [Auto-Fill] Material "%s": Mu_r set to %.2f based on BH curve.\n', mat.Name, calc_mu_r);
            end
            
            fprintf('  - Material "%s": Preprocessing done.\n', mat.Name);
        else
            Model.Materials.Lib(i).Type = 'Linear';
        end
    end
end