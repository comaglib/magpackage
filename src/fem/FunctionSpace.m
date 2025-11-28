classdef FunctionSpace
    % FUNCTIONSPACE 定义有限元离散空间
    % 例如: Lagrange P1, Nedelec P1
    
    properties
        Type  % 'Lagrange' or 'Nedelec'
        Order % 1 (目前仅支持一阶)
    end
    
    methods
        function obj = FunctionSpace(type, order)
            if nargin > 0
                obj.Type = type;
                obj.Order = order;
            end
        end
        
        function str = toString(obj)
            str = sprintf('%s_P%d', obj.Type, obj.Order);
        end
    end
end