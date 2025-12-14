classdef MaterialLib
    % MATERIALLIB 无状态材料库 (v3.8 - Physics-Based Extrapolation Fix)
    % 
    % 变更日志:
    % 1. [Critical Fix] 引入"真空切线修正" (Vacuum Tangent Correction)。
    %    原因: 简单的 B-H 数据往往缺少深饱和段。如果在 B-H 末端直接线性外推 nu(B^2)，
    %    斜率通常过小(过软)，导致仿真电流峰值远低于实际(如 COMSOL)。
    %    方案: 代码自动计算末端点趋向真空的理论斜率 d(nu)/d(B^2) = (nu0 - nu)/(2*B^2)。
    %    如果数据本身的斜率小于该理论值，代码会自动追加一个"物理修正点"，
    %    强制外推斜率符合真空物理特性。
    % 2. [Compatibility] 保持 Linear Interpolation + Cubic Stride 结构不变。
    
    methods (Static)
        function [nu, dnu_db2] = evaluate(B_sq, matData)
            % EVALUATE 计算磁阻率 nu 及其对 B^2 的导数 (保持 v3.7 逻辑不变)
            
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
                
                % 2. 插值区
                if any(mask_in)
                    nu(mask_in) = ppval(matData.SplineNu, B_sq(mask_in));
                    dnu_db2(mask_in) = ppval(matData.SplineDNu, B_sq(mask_in));
                end
                
                % 3. 外推区 (线性)
                if any(mask_ext)
                    x_end = matData.MaxBSq;
                    val_end = ppval(matData.SplineNu, x_end);
                    slope_end = ppval(matData.SplineDNu, x_end); % 这里将获取到修正后的陡峭斜率
                    
                    delta_x = B_sq(mask_ext) - x_end;
                    nu(mask_ext) = val_end + slope_end * delta_x;
                    dnu_db2(mask_ext) = slope_end;
                end
                
                % 4. 物理限幅
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
            
            % 1. 基础预处理
            B_curve = abs(B_curve(:)); H_curve = abs(H_curve(:));
            [B_curve, sortIdx] = sort(B_curve); H_curve = H_curve(sortIdx);
            [B_unique, uniqueIdx] = unique(B_curve, 'last'); % 使用 'last' 保留更饱和的点
            H_unique = H_curve(uniqueIdx);
            
            B_sq = B_unique.^2;
            Nu = zeros(size(B_unique));
            mask_nz = B_unique > 1e-12;
            Nu(mask_nz) = H_unique(mask_nz) ./ B_unique(mask_nz);
            if ~mask_nz(1), if length(B_unique)>1, Nu(1)=Nu(2); else, Nu(1)=nu0; end; end
            Nu(Nu > nu0) = nu0;
            
            % =============================================================
            % [核心修正] 自动追加真空切线点 (Physics-Correction Point)
            % =============================================================
            if length(B_sq) > 2
                % A. 计算数据本身的末端斜率 (割线)
                dB2 = B_sq(end) - B_sq(end-1);
                dNu = Nu(end) - Nu(end-1);
                slope_data = dNu / dB2;
                
                % B. 计算理论真空切线斜率
                % nu(B) = H/B -> d(nu)/d(B^2) = (nu0 - nu) / (2*B^2)
                B_end = B_unique(end);
                Nu_end = Nu(end);
                slope_vacuum = (nu0 - Nu_end) / (2 * B_end^2);
                
                % C. 判定与修正
                % 如果数据斜率显著小于理论斜率 (例如小于 50%)，说明不够"硬"
                if slope_data < 0.8 * slope_vacuum && slope_vacuum > 1e-10
                    % 追加一个极近的点，强制拉高斜率
                    delta_B_sq_add = 0.01; % B^2 增加 0.01 (很小)
                    B_sq_new = B_sq(end) + delta_B_sq_add;
                    
                    % 使用理论斜率计算新 Nu 值
                    Nu_new = Nu(end) + slope_vacuum * delta_B_sq_add;
                    
                    % 防止超限
                    if Nu_new > nu0, Nu_new = nu0; end
                    
                    % 拼接到数组
                    B_sq = [B_sq; B_sq_new];
                    Nu   = [Nu;   Nu_new];
                    
                    % 更新 MaxBSq 到这个新点
                    % 这样 C++ 代码就会使用这个新线段作为外推基准
                    % fprintf('   [MaterialLib] Auto-corrected saturation slope: %.2e -> %.2e\n', slope_data, slope_vacuum);
                end
            end
            % =============================================================
            
            matData.MaxBSq = max(B_sq);
            matData.Type = 'Nonlinear';
            
            % 构建兼容 C++ 的伪三次样条
            if length(B_sq) > 1
                pp_linear = interp1(B_sq, Nu, 'linear', 'pp');
                coefs_lin = pp_linear.coefs;
                [num_segments, ~] = size(coefs_lin);
                coefs_cubic = [zeros(num_segments, 2), coefs_lin];
                matData.SplineNu = mkpp(pp_linear.breaks, coefs_cubic);
            else
                matData.SplineNu = mkpp(B_sq, [0 0 0 Nu(1)]); 
            end
            
            matData.SplineDNu = fnder(matData.SplineNu, 1);
            matData.nu0 = nu0;
        end
    end
end