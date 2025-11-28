classdef MaterialLib
    % MATERIALLIB 无状态材料库 (v3.3 - Clamped Safe)
    % 
    % 修复:
    % 1. 增加 spline 范围检查，防止 B 场过大时的样条外插发散。
    % 2. 保证饱和区导数行为的物理一致性。
    
    methods (Static)
        function [nu, dnu_db2] = evaluate(B_sq, matData)
            % 提取物理常数
            if isfield(matData, 'nu0')
                nu0 = matData.nu0;
            else
                nu0 = 1.0 / (4*pi*1e-7); 
            end

            if strcmp(matData.Type, 'Linear')
                val = matData.Nu_Linear;
                if isscalar(B_sq)
                    nu = val; dnu_db2 = 0;
                else
                    nu = ones(size(B_sq)) * val; dnu_db2 = zeros(size(B_sq));
                end
            else
                % [关键修复] 输入截断
                % 防止 B_sq 超出样条定义的范围导致外插数值爆炸
                % 对于超出范围的 B，我们认为材料性质保持为最后一点的状态 (完全饱和)
                B_sq_clamped = B_sq;
                
                if isfield(matData, 'MaxBSq')
                    mask_high = B_sq > matData.MaxBSq;
                    if any(mask_high)
                        B_sq_clamped(mask_high) = matData.MaxBSq;
                    end
                end
                
                % 计算 Spline
                nu = ppval(matData.SplineNu, B_sq_clamped);
                dnu_db2 = ppval(matData.SplineDNu, B_sq_clamped);
                
                % [额外防御] 对于超出范围的点，导数设为 0 (假设进入线性饱和区)
                % 或者保留最后一点的导数？保留最后一点导数更符合 Newton 法的光滑要求
                % 但物理上，深度饱和后 nu 趋于 nu0 (常数)，所以导数趋于 0 是对的
                % 这里我们保持 Spline 的最后导数，防止雅可比矩阵突变
                
                % 安全防御: nu 不超过 nu0
                mask_overflow = nu > nu0;
                if any(mask_overflow)
                    nu(mask_overflow) = nu0;
                    dnu_db2(mask_overflow) = 0; 
                end
                
                % 安全防御: nu 不为负
                mask_neg = nu < 0;
                if any(mask_neg)
                    nu(mask_neg) = nu0; 
                    dnu_db2(mask_neg) = 0;
                end
            end
        end
        
        function matData = createLinear(mu_r)
            mu0 = 4*pi*1e-7;
            matData.Type = 'Linear';
            matData.Nu_Linear = 1.0 / (mu_r * mu0);
            matData.nu0 = 1.0 / mu0;
        end
        
        function matData = createNonlinear(B_curve, H_curve)
            mu0 = 4*pi*1e-7;
            nu0 = 1.0 / mu0;
            
            B_curve = B_curve(:);
            H_curve = H_curve(:);
            
            % 记录定义域上限
            matData.MaxBSq = max(B_curve.^2);
            
            B_sq = B_curve.^2;
            Nu = zeros(size(B_curve));
            
            mask_nonzero = B_curve > 1e-12;
            Nu(mask_nonzero) = H_curve(mask_nonzero) ./ B_curve(mask_nonzero);
            
            if ~mask_nonzero(1)
                if length(B_curve) > 1
                    Nu(1) = H_curve(2) / B_curve(2);
                else
                    Nu(1) = nu0;
                end
            end
            
            Nu(Nu > nu0) = nu0;
            
            matData.Type = 'Nonlinear';
            matData.SplineNu = pchip(B_sq, Nu);
            matData.SplineDNu = fnder(matData.SplineNu, 1);
            matData.nu0 = nu0;
        end
    end
end