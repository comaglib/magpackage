function [nu, dnu] = eval_material_nu(B_sq, MatInfo)
% EVAL_MATERIAL_NU 计算磁阻率及其导数 (v2: 含安全截断)
% 
% 增强:
% 1. 强制 nu <= nu0 (真空磁阻率)
% 2. 当 nu 被截断时，强制 dnu = 0 (真空是线性的)

    if strcmp(MatInfo.Type, 'Linear')
        mu0 = 4*pi*1e-7;
        val = 1.0 / (MatInfo.Mu_r * mu0);
        
        if isscalar(B_sq)
            nu = val; dnu = 0;
        else
            nu = repmat(val, size(B_sq));
            dnu = zeros(size(B_sq));
        end
    else
        % 非线性样条
        pp_nu = MatInfo.NonlinearData.pp_nu;
        pp_dnu = MatInfo.NonlinearData.pp_dnu;
        
        nu = ppval(pp_nu, B_sq);
        dnu = ppval(pp_dnu, B_sq);
        
        % [安全防御] 硬约束: 任何材料的磁阻率不应超过真空
        % (除非是抗磁性材料，但工程电磁场通常忽略抗磁性)
        nu0 = MatInfo.NonlinearData.nu0;
        
        % 找出越界的点
        mask_overflow = nu > nu0;
        
        if any(mask_overflow)
            nu(mask_overflow) = nu0;
            dnu(mask_overflow) = 0; % 饱和区的导数为 0 (nu是常数)
        end
        
        % [安全防御] 硬约束: nu 不能为负 (虽然 pchip 应该保证了)
        mask_neg = nu < 0;
        if any(mask_neg)
            nu(mask_neg) = nu0; % 回退到空气，防止崩溃
            dnu(mask_neg) = 0;
        end
    end
end