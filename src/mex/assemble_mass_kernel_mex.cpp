/**
 * assemble_mass_kernel_mex.cpp
 * 功能: 组装 Nedelec 单元的质量矩阵 M = (sigma * Ni, Nj)
 * * 修正日志:
 * 1. [Fix Indexing] 修正 RefVal (基函数值) 的索引逻辑，匹配 [6 x 3 x Nq] 布局。
 * 2. [Fix Weight] 移除冗余的 / 6.0 因子。
 * 3. [Compat] 使用 mxGetPr。
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include "MexElemUtils.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Check Inputs
    if (nrhs != 7) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 7 inputs.");

    // Inputs: P, T, Dofs, Signs, Coeffs(Optional), Qw, RefVal
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* Dofs = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* Coeffs = mxGetPr(prhs[4]); // Can be NULL/Empty if logic handles it
    double* Qw = mxGetPr(prhs[5]);
    double* RefVal = mxGetPr(prhs[6]); // [6 x 3 x NumQP]

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[5]);
    
    // Check if Coeffs is provided and has correct size
    bool hasCoeffs = (mxGetNumberOfElements(prhs[4]) == numElems);

    // Outputs
    size_t numTriplets = numElems * 36;
    plhs[0] = mxCreateDoubleMatrix(numTriplets, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(numTriplets, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(numTriplets, 1, mxREAL);
    
    double* I_out = mxGetPr(plhs[0]);
    double* J_out = mxGetPr(plhs[1]);
    double* V_out = mxGetPr(plhs[2]);

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
        
        // Compute J^{-1} for Piola transform: N = J^{-T} * N_ref
        MexElemUtils::compute_inverse_3x3(J, detJ, J_inv);

        double M_local[6][6] = {0};
        double sigma = hasCoeffs ? Coeffs[e] : 1.0;
        
        // Skip calculation if sigma is effectively zero
        if (std::abs(sigma) > 1e-12) {
            
            for (int q = 0; q < numQP; q++) {
                // [Fix Weight] No / 6.0
                double w = Qw[q] * absDetJ;
                
                // Transform Basis Functions
                double N_phy[3][6]; 
                long q_offset = q * 18; // 6 rows * 3 cols * q

                for (int i = 0; i < 6; i++) {
                    double ref_n[3];
                    // [Fix Indexing]
                    // Layout: [6 rows (basis) x 3 cols (xyz) x Nq]
                    // Index = row + 6*col + 18*q
                    ref_n[0] = RefVal[i + 0  + q_offset];
                    ref_n[1] = RefVal[i + 6  + q_offset];
                    ref_n[2] = RefVal[i + 12 + q_offset];
                    
                    // Piola Transform: N_phy = J^{-T} * ref_n
                    // Equivalent to: (J^{-1})^T * ref_n
                    // Using helper: multVecTrans(input, output) using J_inv
                    double temp[3];
                    J_inv.multVecTrans(ref_n, temp);
                    
                    N_phy[0][i] = temp[0];
                    N_phy[1][i] = temp[1];
                    N_phy[2][i] = temp[2];
                }
                
                // Local Assembly
                for (int i = 0; i < 6; i++) {
                    for (int j = 0; j < 6; j++) {
                        double dot = N_phy[0][i]*N_phy[0][j] + 
                                     N_phy[1][i]*N_phy[1][j] + 
                                     N_phy[2][i]*N_phy[2][j];
                        M_local[i][j] += dot * sigma * w;
                    }
                }
            }
        }

        // Output Filling
        double signs[6];
        for(int k=0; k<6; k++) signs[k] = Signs[k + e*6];
        long dofs[6];
        for(int k=0; k<6; k++) dofs[k] = (long)Dofs[k + e*6];

        long base_idx = e * 36;
        int count = 0;
        for (int c = 0; c < 6; c++) {
            for (int r = 0; r < 6; r++) {
                I_out[base_idx + count] = (double)dofs[r];
                J_out[base_idx + count] = (double)dofs[c];
                V_out[base_idx + count] = M_local[r][c] * signs[r] * signs[c];
                count++;
            }
        }
    }
}