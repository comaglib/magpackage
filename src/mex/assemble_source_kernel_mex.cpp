/**
 * assemble_source_kernel_mex.cpp
 * 功能: 组装 Nedelec 单元的源向量 F = (J_s, Ni)
 * * 修正日志:
 * 1. [Fix Indexing] 修正 RefVal 索引逻辑。
 * 2. [Fix Weight] 移除 / 6.0 因子。
 * 3. [Compat] 使用 mxGetPr。
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include "MexElemUtils.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 7) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 7 inputs.");

    // Inputs
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* Dofs = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* Js = mxGetPr(prhs[4]); // [3 x NumElems]
    double* Qw = mxGetPr(prhs[5]);
    double* RefVal = mxGetPr(prhs[6]); // [6 x 3 x NumQP]

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[5]);
    
    // Check Js dimensions
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != numElems) {
        mexErrMsgIdAndTxt("MagPackage:Dim", "Js must be [3 x NumElems]");
    }

    // Output: Sparse triplet format (I, V) but logic is element-wise dense
    // Each element contributes to 6 DOFs
    size_t numEntries = numElems * 6;
    plhs[0] = mxCreateDoubleMatrix(numEntries, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(numEntries, 1, mxREAL); // V
    
    double* I_out = mxGetPr(plhs[0]);
    double* V_out = mxGetPr(plhs[1]);

    #pragma omp parallel for
    for (long e = 0; e < numElems; e++) {
        
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

        // Get Current Density for this element (Physical Space)
        double J_phys[3];
        J_phys[0] = Js[0 + e*3];
        J_phys[1] = Js[1 + e*3];
        J_phys[2] = Js[2 + e*3];

        // Optimization: Pre-transform J_phys
        // Integrand: J_phys . N_phys
        // N_phys = J^{-T} * N_ref
        // Dot Product property: A . (M * B) = (M^T * A) . B
        // J_phys . (J^{-T} * N_ref) = ((J^{-T})^T * J_phys) . N_ref
        //                           = (J^{-1} * J_phys) . N_ref
        double J_trans[3];
        J_inv.multVec(J_phys, J_trans);

        double F_local[6] = {0};
        
        // Skip if source is zero
        if (std::abs(J_phys[0]) > 1e-12 || std::abs(J_phys[1]) > 1e-12 || std::abs(J_phys[2]) > 1e-12) {
            
            for (int q = 0; q < numQP; q++) {
                // [Fix Weight]
                double w = Qw[q] * absDetJ;
                long q_offset = q * 18;

                for (int i = 0; i < 6; i++) {
                    double ref_n[3];
                    // [Fix Indexing]
                    // Layout: [6 rows (basis) x 3 cols (xyz) x Nq]
                    ref_n[0] = RefVal[i + 0  + q_offset];
                    ref_n[1] = RefVal[i + 6  + q_offset];
                    ref_n[2] = RefVal[i + 12 + q_offset];
                    
                    // Dot product: J_trans . N_ref
                    double val = J_trans[0]*ref_n[0] + J_trans[1]*ref_n[1] + J_trans[2]*ref_n[2];
                    F_local[i] += val * w;
                }
            }
        }

        // Fill Output
        double signs[6];
        for(int k=0; k<6; k++) signs[k] = Signs[k + e*6];
        long dofs[6];
        for(int k=0; k<6; k++) dofs[k] = (long)Dofs[k + e*6];

        long base_idx = e * 6;
        for (int i = 0; i < 6; i++) {
            I_out[base_idx + i] = (double)dofs[i];
            V_out[base_idx + i] = F_local[i] * signs[i];
        }
    }
}