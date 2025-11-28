classdef Assembler < handle
    % ASSEMBLER 有限元组装器 (v5.0 - Winding Support)
    
    properties
        Mesh
        DofHandler
        Config
    end
    
    methods
        function obj = Assembler(mesh, dofHandler, config)
            obj.Mesh = mesh;
            obj.DofHandler = dofHandler;
            if nargin < 3
                obj.Config = FemConfig.Default();
            else
                obj.Config = config;
            end
        end
        
        function K = assembleStiffness(obj, space, materialMap)
            packedData = obj.preparePackedData(space);
            elemTags = obj.Mesh.RegionTags;
            try
                Nu_vec = materialMap(elemTags);
            catch
                error('MaterialMap Error (Stiffness).');
            end
            packedData.Nu = Nu_vec(:);
            
            if strcmpi(space.Type, 'Nedelec')
                [I, J, V] = assemble_curl_curl_kernel(packedData, obj.Config);
                K = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space.');
            end
        end
        
        function M = assembleMass(obj, space)
            packedData = obj.preparePackedData(space);
            if strcmpi(space.Type, 'Nedelec')
                [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                M = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space.');
            end
        end
        
        function M = assembleMassWeighted(obj, space, coeffMap)
            packedData = obj.preparePackedData(space);
            elemTags = obj.Mesh.RegionTags;
            try
                Coeff_vec = coeffMap(elemTags);
            catch
                error('MaterialMap Error (Mass).');
            end
            packedData.Coeff = Coeff_vec(:);
            
            if strcmpi(space.Type, 'Nedelec')
                [I, J, V] = assemble_mass_kernel(packedData, obj.Config);
                M = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
            else
                error('Unsupported space.');
            end
        end
        
        function F = assembleSource(obj, space, sourceMap)
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            F = assemble_source_kernel(packedData, sourceMap, obj.Config);
        end
        
        % [新增] 组装绕组耦合向量
        function C = assembleWinding(obj, space, windingObj)
            % ASSEMBLEWINDING 组装绕组耦合向量 C
            % packedData = obj.preparePackedData(space);
            % packedData.RegionTags = obj.Mesh.RegionTags;
            % 为了减少重复代码，直接调用 preparePackedData
            
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            
            C = assemble_winding_kernel(packedData, windingObj, obj.Config);
        end
        
        function [J_mat, R_mat] = assembleJacobian(obj, space, solutionA, matLibData, calcJ)
            if nargin < 5, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            if strcmpi(space.Type, 'Nedelec')
                [I, J, V, R_vec] = assemble_jacobian_kernel(packedData, solutionA, matLibData, obj.Config, calcJ);
                if calcJ
                    J_mat = sparse(I, J, V, obj.DofHandler.NumGlobalDofs, obj.DofHandler.NumGlobalDofs);
                else
                    J_mat = [];
                end
                R_mat = R_vec; 
            else
                error('Unsupported space.');
            end
        end
        
        function [J_triplets, Res_mat] = assembleHBFEM(obj, space, solHarmonics, aftObj, matLibData, calcJ)
            if nargin < 6, calcJ = true; end
            packedData = obj.preparePackedData(space);
            packedData.RegionTags = obj.Mesh.RegionTags;
            if strcmpi(space.Type, 'Nedelec')
                [I, J, V, R_mat] = assemble_hbfem_kernel(packedData, solHarmonics, aftObj, matLibData, obj.Config, calcJ);
                if calcJ
                    J_triplets.I = I; J_triplets.J = J; J_triplets.V = V;
                else
                    J_triplets = [];
                end
                Res_mat = R_mat;
            else
                error('Unsupported space.');
            end
        end
    end
    
    methods (Access = private)
        function packedData = preparePackedData(obj, space)
            packedData = obj.DofHandler.packForKernel(space);
            packedData.P = obj.Mesh.P;
            packedData.T = obj.Mesh.T;
        end
    end
end