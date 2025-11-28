classdef Winding < handle
    % WINDING 线圈绕组类 (v3.0)
    % 定义线圈的几何区域、匝数、电阻和方向
    
    properties
        Name            % 名称 (String)
        RegionID        % 对应的网格区域 ID (Integer)
        Turns           % 匝数 N (Float)
        Resistance      % 直流电阻 R (Ohm)
        CrossSectionArea % 导线束总截面积 S (m^2)
        Direction       % 电流方向向量 [x, y, z] (假设线圈段是直的，或由网格几何决定)
    end
    
    methods
        function obj = Winding(name, id, N, R, S, dir)
            if nargin > 0
                obj.Name = name;
                obj.RegionID = id;
                obj.Turns = N;
                obj.Resistance = R;
                obj.CrossSectionArea = S;
                if norm(dir) > 0
                    obj.Direction = reshape(dir, 1, 3) / norm(dir); % 归一化方向
                else
                    obj.Direction = [0 0 0];
                end
            end
        end
        
        function J_eff = getEffectiveCurrentDensity(obj)
            % 获取单位电流下的等效电流密度向量
            % J_eff = (N / S) * direction
            % 物理意义: 当线圈通入 1A 电流时，区域内的电流密度 J
            if obj.CrossSectionArea <= 0
                error('Winding Area must be positive.');
            end
            J_eff = (obj.Turns / obj.CrossSectionArea) * obj.Direction;
        end
    end
end