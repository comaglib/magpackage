/**
 * assemble_winding_kernel_mex.cpp
 * 功能: 组装线圈电流源向量 F = (Scale * Dir * N_i)
 * 输入: P, T, Dofs, Signs, RegionTags, CoilID, DirField, UniformDir, Scale, Qw, RefVal
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include "MexElemUtils.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 11) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 11 inputs.");

    // Inputs
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* Dofs = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* Tags = mxGetPr(prhs[4]);
    double  CoilID = mxGetScalar(prhs[5]);
    
    // DirField: [3 x NumElems] (Optional, can be empty)
    double* DirField = NULL;
    if (!mxIsEmpty(prhs[6])) DirField = mxGetPr(prhs[6]);
    
    // UniformDir: [3 x 1]
    double* UniDir = mxGetPr(prhs[7]);
    double  Scale = mxGetScalar(prhs[8]);
    
    double* Qw = mxGetPr(prhs[9]);
    double* RefVal = mxGetPr(prhs[10]);

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[9]);
    
    // Outputs
    size_t numEntries = numElems * 6;
    plhs[0] = mxCreateDoubleMatrix(numEntries, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(numEntries, 1, mxREAL); // V
    
    double* I_out = mxGetPr(plhs[0]);
    double* V_out = mxGetPr(plhs[1]);

    #pragma omp parallel for
    for (long e = 0; e < numElems; e++) {
        
        // Output pointer for this element
        long base_idx = e * 6;
        
        // Check Region ID
        if ((int)Tags[e] != (int)CoilID) {
            // Fill zeros if not in coil
            for(int i=0; i<6; i++) {
                I_out[base_idx + i] = 0; // or 1, doesn't matter as V is 0
                V_out[base_idx + i] = 0; 
            }
            continue;
        }

        // 1. Geometry
        double p_local[12];
        for (int n = 0; n < 4; n++) {
            long node_idx = (long)T[n + e*4] - 1; 
            p_local[n*3 + 0] = P[node_idx*3 + 0];
            p_local[n*3 + 1] = P[node_idx*3 + 1];
            p_local[n*3 + 2] = P[node_idx*3 + 2];
        }

        Mat3x3 J, J_inv;
        double detJ = MexElemUtils::compute_jacobian_3d(p_local, J);
        double absDetJ = std::abs(detJ);
        MexElemUtils::compute_inverse_3x3(J, detJ, J_inv);

        // Determine Current Direction J_dir
        double J_phys[3];
        if (DirField != NULL) {
            J_phys[0] = DirField[0 + e*3];
            J_phys[1] = DirField[1 + e*3];
            J_phys[2] = DirField[2 + e*3];
        } else {
            J_phys[0] = UniDir[0];
            J_phys[1] = UniDir[1];
            J_phys[2] = UniDir[2];
        }
        
        // Apply Scaling (N_turns * Current / Area)
        J_phys[0] *= Scale;
        J_phys[1] *= Scale;
        J_phys[2] *= Scale;

        // Optimization: J_trans = J^{-1} * J_phys
        double J_trans[3];
        J_inv.multVec(J_phys, J_trans);

        double F_local[6] = {0};
        
        for (int q = 0; q < numQP; q++) {
            double w = Qw[q] * absDetJ;
            long q_offset = q * 18;

            for (int i = 0; i < 6; i++) {
                double ref_n[3];
                // Nedelec Ref Basis [6 x 3 x Nq]
                ref_n[0] = RefVal[i + 0  + q_offset];
                ref_n[1] = RefVal[i + 6  + q_offset];
                ref_n[2] = RefVal[i + 12 + q_offset];
                
                double val = J_trans[0]*ref_n[0] + J_trans[1]*ref_n[1] + J_trans[2]*ref_n[2];
                F_local[i] += val * w;
            }
        }

        // Fill Output
        double signs[6];
        for(int k=0; k<6; k++) signs[k] = Signs[k + e*6];
        long dofs[6];
        for(int k=0; k<6; k++) dofs[k] = (long)Dofs[k + e*6];

        for (int i = 0; i < 6; i++) {
            I_out[base_idx + i] = (double)dofs[i];
            V_out[base_idx + i] = F_local[i] * signs[i];
        }
    }
}