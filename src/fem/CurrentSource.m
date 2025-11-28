classdef CurrentSource < handle
    % CURRENTSOURCE 定义电流源项 J
    % 目前支持: 基于区域 ID 的常数向量 J
    
    properties
        RegionID    % 作用的物理区域 ID (整数)
        J_Vector    % 电流密度向量 [Jx, Jy, Jz] (对于 Domain 类型)
        Type        % 'Domain'
    end
    
    methods
        function obj = CurrentSource(regionID, J_vec)
            % 构造函数
            % obj = CurrentSource(5, [0, 0, 1e6]); % 在区域5施加 z方向 1e6 A/m^2
            obj.Type = 'Domain';
            obj.RegionID = regionID;
            
            if length(J_vec) ~= 3
                error('J_Vector must be a 3-element vector.');
            end
            obj.J_Vector = reshape(J_vec, 1, 3); % 确保是行向量
        end
    end
end