classdef LossCalculator < handle
    % LOSSCALCULATOR 损耗计算模块 (v2.0 - HBFEM Support)
    % 
    % 功能:
    %   1. 计算实体导体的涡流损耗 (Ohmic Loss)。
    %   2. 计算铁磁材料的铁损 (Iron Loss - Steinmetz)。
    %   3. 支持频域 (Frequency) 和 谐波平衡 (HBFEM) 结果。
    
    properties
        Assembler
        PostProcessor
    end
    
    methods
        function obj = LossCalculator(assembler)
            obj.Assembler = assembler;
            obj.PostProcessor = PostProcessor(assembler);
        end
        
        function P_ohmic = computeOhmicLoss_Frequency(obj, A_phasor, sigmaMap, freq, space)
            % COMPUTEOHMICLOSS_FREQUENCY 计算频域下的涡流损耗
            % Formula: P = 1/2 * w^2 * (A' * M_sigma * A)
            
            omega = 2 * pi * freq;
            M_sigma = obj.Assembler.assembleMassWeighted(space, sigmaMap);
            term = A_phasor' * M_sigma * A_phasor;
            P_ohmic = 0.5 * (omega^2) * real(term);
        end
        
        function P_ohmic = computeOhmicLoss_HBFEM(obj, X_harmonics, sigmaMap, aft, space)
            % COMPUTEOHMICLOSS_HBFEM 计算 HBFEM 下的总涡流损耗
            % Formula: P = Sum( 1/2 * (k*w)^2 * A_k' * M_sigma * A_k )
            
            base_omega = 2 * pi * aft.BaseFreq;
            harmonics = aft.Harmonics;
            
            M_sigma = obj.Assembler.assembleMassWeighted(space, sigmaMap);
            P_ohmic = 0;
            
            for k = 1:length(harmonics)
                h = harmonics(k);
                if h == 0, continue; end 
                
                w = h * base_omega;
                A_k = X_harmonics(:, k);
                
                term = A_k' * M_sigma * A_k;
                P_k = 0.5 * (w^2) * term;
                P_ohmic = P_ohmic + P_k;
            end
        end
        
        function P_iron = computeIronLoss_Steinmetz(obj, A_sol, freq, k_h, alpha, k_e, space)
            % COMPUTEIRONLOSS_STEINMETZ (单频版)
            spaceName = space.toString();
            B_elems = obj.PostProcessor.computeElementB(A_sol, spaceName);
            B_mag = obj.PostProcessor.computeMagnitude(B_elems); 
            Vols = obj.computeElementVolumes(); 
            
            if size(B_mag, 2) == 1
                B_peak = B_mag;
                p_dens = k_h * freq * (B_peak .^ alpha) + k_e * (freq^2) * (B_peak .^ 2);
                P_iron = sum(p_dens .* Vols);
            else 
                warning('Multi-harmonic input detected in single-freq method. Using Column 1.');
                B_peak = B_mag(:, 1);
                p_dens = k_h * freq * (B_peak .^ alpha) + k_e * (freq^2) * (B_peak .^ 2);
                P_iron = sum(p_dens .* Vols);
            end
        end
        
        function P_total = computeIronLoss_HBFEM(obj, A_sol, aft, k_h, alpha, k_e, space)
            % COMPUTEIRONLOSS_HBFEM (多谐波叠加版)
            % 原理: 损耗功率是标量，各个频率分量的损耗近似线性叠加 (正交假设)
            
            base_freq = aft.BaseFreq;
            harmonics = aft.Harmonics;
            spaceName = space.toString();
            
            Vols = obj.computeElementVolumes();
            P_total = 0;
            
            fprintf('   [Loss] Calculating HBFEM Iron Loss...\n');
            
            for k = 1:length(harmonics)
                h = harmonics(k);
                if h == 0, continue; end % DC 通常无磁滞/涡流损耗(或极小)
                
                freq = h * base_freq;
                A_k = A_sol(:, k);
                
                % 计算该谐波下的 B 场
                B_elems = obj.PostProcessor.computeElementB(A_k, spaceName);
                B_mag = obj.PostProcessor.computeMagnitude(B_elems); % [Ne x 1]
                
                % 应用 Steinmetz 公式
                p_dens = k_h * freq * (B_mag .^ alpha) + k_e * (freq^2) * (B_mag .^ 2);
                
                P_k = sum(p_dens .* Vols);
                P_total = P_total + P_k;
                
                % fprintf('      H%d (%.0f Hz): %.4e W\n', h, freq, P_k);
            end
        end
        
        function vols = computeElementVolumes(obj)
            P = obj.Assembler.Mesh.P;
            T = obj.Assembler.Mesh.T;
            v1 = P(:, T(1,:)); v2 = P(:, T(2,:)); v3 = P(:, T(3,:)); v4 = P(:, T(4,:));
            d1 = v2 - v1; d2 = v3 - v1; d3 = v4 - v1;
            c = cross(d1, d2, 1);
            dot_prod = sum(c .* d3, 1);
            vols = abs(dot_prod)' / 6.0;
        end
        
        function P_loss = computeBertottiLoss_TimeStep(obj, dA_dt_col, matLib)
            % COMPUTEBERTOTTILOSS_TIMESTEP 计算单步瞬时损耗
            % 基于 Bertotti 分离模型: P_total = P_classic + P_excess + P_hyst
            
            % 1. 计算 dB/dt 场
            %    B = curl(A) => dB/dt = curl(dA/dt)
            %    [修复] 显式指定 'Nedelec_P1' 空间，防止自动匹配到 Lagrange 空间导致索引越界
            dB_dt_elems = obj.PostProcessor.computeElementB(dA_dt_col, 'Nedelec_P1'); % [3 x Ne]
            
            % 2. 计算模长 |dB/dt|
            dB_dt_mag = obj.PostProcessor.computeMagnitude(dB_dt_elems); % [Ne x 1]
            
            % 3. 获取单元体积
            Vols = obj.computeElementVolumes();
            
            % 4. 遍历材料区域进行积分
            P_classic = 0;
            P_excess = 0;
            
            mesh = obj.Assembler.Mesh;
            tags = mesh.RegionTags;
            uTags = unique(tags);
            
            for k = 1:length(uTags)
                tag = uTags(k);
                
                % 获取材料数据
                if isa(matLib, 'containers.Map') && matLib.isKey(tag)
                    mat = matLib(tag);
                elseif isnumeric(matLib) && tag <= length(matLib)
                    mat = matLib(tag);
                else
                    continue; 
                end
                
                % 仅对定义了叠片参数的材料计算
                if isfield(mat, 'Lam_d') && isfield(mat, 'Lam_sigma')
                    mask = (tags == tag);
                    if ~any(mask), continue; end
                    
                    vol_sub = Vols(mask);
                    dB_dt_sub = dB_dt_mag(mask);
                    
                    % === A. 经典涡流损耗 (Classical Eddy Current Loss) ===
                    % 公式: P = (sigma * d^2 / 12) * (dB/dt)^2 * Vol
                    coeff_cl = (mat.Lam_sigma * mat.Lam_d^2) / 12.0;
                    p_dens_cl = coeff_cl * (dB_dt_sub .^ 2);
                    P_classic = P_classic + sum(p_dens_cl .* vol_sub);
                    
                    % === B. 异常损耗 (Excess Loss) ===
                    if isfield(mat, 'Loss_Kexc') && ~isempty(mat.Loss_Kexc) && mat.Loss_Kexc > 0
                        p_dens_exc = mat.Loss_Kexc * (dB_dt_sub .^ 1.5);
                        P_excess = P_excess + sum(p_dens_exc .* vol_sub);
                    end
                end
            end
            
            P_loss.Classic = P_classic;
            P_loss.Excess = P_excess;
            P_loss.Total = P_classic + P_excess;
        end
    end
end