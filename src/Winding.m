classdef Winding < handle
    % WINDING 线圈绕组类 (v3.1 - Spatially Varying Direction)
    
    properties
        Name            
        RegionID        
        Turns           
        Resistance      
        CrossSectionArea 
        Direction       % 全局统一方向 (备用)
        
        % [New] 单元级方向场 (优先使用)
        % [3 x NumElements] 矩阵。如果非空，内核将优先读取此数据。
        DirectionField  
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
                    obj.Direction = reshape(dir, 1, 3) / norm(dir); 
                else
                    obj.Direction = [0 0 0];
                end
                obj.DirectionField = []; % 默认为空
            end
        end
        
        function setDirectionField(obj, field)
            % 设置空间变化的电流方向
            obj.DirectionField = field;
        end
    end
end