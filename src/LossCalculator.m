classdef LossCalculator < handle
    % LOSSCALCULATOR 损耗计算模块 (v1.1 - Fix Interface)
    % 
    % 修复:
    %   1. computeOhmicLoss_* 方法现在接收 'space' 参数，
    %      以正确调用 Assembler.assembleMassWeighted。
    
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
            
            % [Fix] 传入 space 对象，而非 DoF 数组
            M_sigma = obj.Assembler.assembleMassWeighted(space, sigmaMap);
            
            % 矩阵运算计算总损耗
            % A_phasor 是复数列向量, result = A' * M * A
            term = A_phasor' * M_sigma * A_phasor;
            
            P_ohmic = 0.5 * (omega^2) * real(term);
        end
        
        function P_ohmic = computeOhmicLoss_HBFEM(obj, X_harmonics, sigmaMap, aft, space)
            % COMPUTEOHMICLOSS_HBFEM 计算 HBFEM 下的总涡流损耗
            % Formula: P = Sum( 1/2 * (k*w)^2 * A_k' * M_sigma * A_k )
            
            base_omega = 2 * pi * aft.BaseFreq;
            harmonics = aft.Harmonics;
            
            % [Fix] 传入 space 对象
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
            % COMPUTEIRONLOSS_STEINMETZ 计算铁损
            
            % 1. 计算每个单元的 B 场模值
            % 传入 space.toString() 明确指定空间
            spaceName = 'Nedelec_P1';
            if nargin >= 7, spaceName = space.toString(); end
            
            B_elems = obj.PostProcessor.computeElementB(A_sol, spaceName);
            B_mag = obj.PostProcessor.computeMagnitude(B_elems); 
            
            % 2. 计算每个单元的体积
            Vols = obj.computeElementVolumes(); 
            
            P_iron = 0;
            
            % 3. 逐单元积分
            if size(B_mag, 2) == 1
                B_peak = B_mag;
                p_dens = k_h * freq * (B_peak .^ alpha) + k_e * (freq^2) * (B_peak .^ 2);
                P_iron = sum(p_dens .* Vols);
            else 
                warning('Multi-harmonic input detected. Using Column 1 only.');
                B_peak = B_mag(:, 1);
                p_dens = k_h * freq * (B_peak .^ alpha) + k_e * (freq^2) * (B_peak .^ 2);
                P_iron = sum(p_dens .* Vols);
            end
        end
        
        function vols = computeElementVolumes(obj)
            P = obj.Assembler.Mesh.P;
            T = obj.Assembler.Mesh.T;
            
            v1 = P(:, T(1,:));
            v2 = P(:, T(2,:));
            v3 = P(:, T(3,:));
            v4 = P(:, T(4,:));
            
            d1 = v2 - v1;
            d2 = v3 - v1;
            d3 = v4 - v1;
            
            c = cross(d1, d2, 1);
            dot_prod = sum(c .* d3, 1);
            
            vols = abs(dot_prod)' / 6.0;
        end
    end
end