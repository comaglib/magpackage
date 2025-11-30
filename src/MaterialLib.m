classdef MaterialLib
    % MATERIALLIB 无状态材料库 (v3.4 - Robust Interpolation)
    % 
    % 修复日志:
    % 1. [Crash Fix] createNonlinear 中增加了 unique() 处理。
    %    原因: 当 B-H 曲线进入深度饱和区 (如 tanh) 时，B 值可能在机器精度下重复，
    %    导致 pchip 插值报错 "The first input must contain unique values"。
    % 2. 强制输入数据为正值并排序，确保单调性。
    
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
                % [输入截断] 防止 B_sq 超出样条范围
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
                
                % [安全防御] 物理约束
                mask_overflow = nu > nu0;
                if any(mask_overflow)
                    nu(mask_overflow) = nu0;
                    dnu_db2(mask_overflow) = 0; 
                end
                
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
            
            % 1. 预处理: 确保列向量且为非负
            B_curve = abs(B_curve(:));
            H_curve = abs(H_curve(:));
            
            % 2. 排序: 按 B 升序排列
            [B_curve, sortIdx] = sort(B_curve);
            H_curve = H_curve(sortIdx);
            
            % 3. [关键修复] 去重
            % pchip 要求自变量严格单调。如果 B 饱和变平，会有重复值。
            % 'first' 策略: 保留进入饱和区的第一个点，丢弃后续重复点。
            % (物理上 B-H 曲线即使饱和也不应完全平坦，应有 mu0 斜率，
            % 但为了支持合成数据，去重是必要的)
            [B_unique, uniqueIdx] = unique(B_curve, 'first');
            H_unique = H_curve(uniqueIdx);
            
            % 记录定义域上限
            matData.MaxBSq = max(B_unique.^2);
            
            B_sq = B_unique.^2;
            Nu = zeros(size(B_unique));
            
            % 计算 Nu = H / B
            mask_nonzero = B_unique > 1e-12;
            Nu(mask_nonzero) = H_unique(mask_nonzero) ./ B_unique(mask_nonzero);
            
            % 处理原点 (0/0) 或极小值
            if ~mask_nonzero(1)
                if length(B_unique) > 1
                    % 简单处理: 假设原点斜率与第二点连续
                    Nu(1) = Nu(2); 
                else
                    Nu(1) = nu0;
                end
            end
            
            % 物理截断: 磁阻率不能超过真空 (相对磁导率不能小于1)
            Nu(Nu > nu0) = nu0;
            
            % 构建样条
            matData.Type = 'Nonlinear';
            matData.SplineNu = pchip(B_sq, Nu);
            matData.SplineDNu = fnder(matData.SplineNu, 1);
            matData.nu0 = nu0;
        end
    end
end