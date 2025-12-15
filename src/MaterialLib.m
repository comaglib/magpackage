classdef MaterialLib
    % MATERIALLIB 材料库 (v4.0 - PCHIP Smoothness)
    % 
    % 核心修复:
    %   1. 强制使用 PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) 
    %      替代线性插值。保证 Nu(B^2) 及其导数 C1 连续。
    %   2. 修复了线性插值导致的 Jacobian 跳变和 SDC 发散问题。
    %   3. 保证外推区与插值区在连接点处切线连续。
    
    methods (Static)
        function [nu, dnu_db2] = evaluate(B_sq, matData)
            % EVALUATE 计算磁阻率 nu 及其对 B^2 的导数
            
            if isfield(matData, 'nu0'), nu0 = matData.nu0; else, nu0 = 1/(4*pi*1e-7); end
            
            if strcmp(matData.Type, 'Linear')
                val = matData.Nu_Linear;
                if isscalar(B_sq), nu = val; dnu_db2 = 0;
                else, nu = ones(size(B_sq)) * val; dnu_db2 = zeros(size(B_sq)); end
            else
                % 1. 区分区域
                if isfield(matData, 'MaxBSq'), mask_ext = B_sq > matData.MaxBSq;
                else, mask_ext = false(size(B_sq)); end
                mask_in = ~mask_ext;
                
                nu = zeros(size(B_sq)); dnu_db2 = zeros(size(B_sq));
                
                % 2. 插值区 (使用 PCHIP 样条)
                if any(mask_in)
                    nu(mask_in) = ppval(matData.SplineNu, B_sq(mask_in));
                    dnu_db2(mask_in) = ppval(matData.SplineDNu, B_sq(mask_in));
                end
                
                % 3. 外推区 (线性外推，但斜率锁定为样条末端切线)
                if any(mask_ext)
                    x_end = matData.MaxBSq;
                    val_end = ppval(matData.SplineNu, x_end);
                    slope_end = ppval(matData.SplineDNu, x_end); 
                    
                    delta_x = B_sq(mask_ext) - x_end;
                    nu(mask_ext) = val_end + slope_end * delta_x;
                    dnu_db2(mask_ext) = slope_end;
                end
                
                % 4. 物理限幅 (防止数值下溢或非物理上溢)
                % 注意: PCHIP 通常不会产生过冲，但外推可能
                mask_overflow = nu > nu0;
                if any(mask_overflow)
                    nu(mask_overflow) = nu0; dnu_db2(mask_overflow) = 0; 
                end
                mask_neg = nu < 1e-6; 
                if any(mask_neg)
                    nu(mask_neg) = nu0; dnu_db2(mask_neg) = 0;
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
            
            % 1. 数据预处理
            B_curve = abs(B_curve(:)); H_curve = abs(H_curve(:));
            [B_curve, sortIdx] = sort(B_curve); H_curve = H_curve(sortIdx);
            [B_unique, uniqueIdx] = unique(B_curve, 'last');
            H_unique = H_curve(uniqueIdx);
            
            % 2. 转换到 Nu - B^2 空间
            B_sq = B_unique.^2;
            Nu = zeros(size(B_unique));
            mask_nz = B_unique > 1e-12;
            Nu(mask_nz) = H_unique(mask_nz) ./ B_unique(mask_nz);
            
            % 处理原点奇异性 (L'Hopital: dH/dB)
            if ~mask_nz(1)
                if length(B_unique) > 1
                    % 用第一个非零点的斜率近似初始磁导率
                    Nu(1) = Nu(2); 
                else
                    Nu(1) = nu0;
                end
            end
            
            % 3. 自动真空修正 (Physics-Correction)
            % 确保数据末端指向真空，防止 PCHIP 在末端 "摆尾"
            if length(B_sq) > 2
                B_end = B_unique(end);
                Nu_end = Nu(end);
                slope_vacuum = (nu0 - Nu_end) / (2 * B_end^2);
                
                % 追加一个非常接近真空的点，引导 PCHIP 的导数方向
                delta_B_sq = 0.1; 
                B_sq_new = B_sq(end) + delta_B_sq;
                Nu_new = Nu(end) + slope_vacuum * delta_B_sq;
                if Nu_new > nu0, Nu_new = nu0; end
                
                B_sq = [B_sq; B_sq_new];
                Nu   = [Nu;   Nu_new];
            end
            
            matData.MaxBSq = max(B_sq);
            matData.Type = 'Nonlinear';
            matData.nu0 = nu0;
            
            % 4. 构建 PCHIP 样条 (关键步骤!)
            % pchip 保证单调性，且生成的分段多项式是 C1 连续的
            if length(B_sq) > 1
                matData.SplineNu = pchip(B_sq, Nu);
            else
                % 单点退化为常数
                matData.SplineNu = mkpp(B_sq, [0 0 0 Nu(1)]);
            end
            
            % 计算导数样条 (用于 Jacobian)
            matData.SplineDNu = fnder(matData.SplineNu, 1);
        end
    end
end