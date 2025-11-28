classdef FemConfig
    % FEMCONFIG 集中管理求解器的底层配置参数
    % 包括并行切片大小、积分精度、容差等
    
    properties
        % 并行组装时的切片大小 (根据内存和CPU核数调整)
        % 较大的值减少调度开销，较小的值平衡负载 50000
        AssemblyChunkSize = 2000
        
        % 默认积分阶数 (1=中心点, 2=4点)
        DefaultQuadratureOrder = 2
        
        % 线性求解器容差
        LinearTolerance = 1e-6
        
        % 调试模式
        DebugMode = false
    end
    
    methods (Static)
        function obj = Default()
            % 快速获取默认配置
            obj = FemConfig();
        end
    end
end