/**
 * assemble_scalar_laplacian_kernel_mex.cpp
 * 功能: 组装标量 Laplacian 矩阵 K = (sigma * grad phi_i, grad phi_j)
 * 输入: P, T, Dofs, Coeffs, Qw, GradRef
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include "MexElemUtils.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 6 inputs.");

    // Inputs
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* Dofs = mxGetPr(prhs[2]); // [4 x NumElems]
    double* Coeffs = mxGetPr(prhs[3]); // [NumElems]
    double* Qw = mxGetPr(prhs[4]);
    double* GradRef = mxGetPr(prhs[5]); // [4 x 3 x Nq]

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[4]);
    
    // Outputs
    // 4 nodes per tet -> 16 entries per elem
    size_t numTriplets = numElems * 16;
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
        MexElemUtils::compute_inverse_3x3(J, detJ, J_inv);

        double K_local[4][4] = {0};
        double sigma = Coeffs[e];
        
        if (std::abs(sigma) > 1e-12) {
            for (int q = 0; q < numQP; q++) {
                double w = Qw[q] * absDetJ; // Weight
                long q_offset = q * 12; // 4 rows * 3 cols

                // Transform Gradients: G_phy = J^{-T} * G_ref
                double G_phy[3][4]; 
                for (int i = 0; i < 4; i++) {
                    double ref_g[3];
                    // Layout: [4 x 3 x Nq] -> i + 4*c + 12*q
                    ref_g[0] = GradRef[i + 0 + q_offset];
                    ref_g[1] = GradRef[i + 4 + q_offset];
                    ref_g[2] = GradRef[i + 8 + q_offset];
                    
                    double temp[3];
                    J_inv.multVecTrans(ref_g, temp);
                    G_phy[0][i] = temp[0];
                    G_phy[1][i] = temp[1];
                    G_phy[2][i] = temp[2];
                }
                
                // Assembly
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        double dot = G_phy[0][i]*G_phy[0][j] + 
                                     G_phy[1][i]*G_phy[1][j] + 
                                     G_phy[2][i]*G_phy[2][j];
                        K_local[i][j] += dot * sigma * w;
                    }
                }
            }
        }

        // Fill Outputs
        long dofs[4];
        for(int k=0; k<4; k++) dofs[k] = (long)Dofs[k + e*4];

        long base_idx = e * 16;
        int count = 0;
        for (int c = 0; c < 4; c++) {
            for (int r = 0; r < 4; r++) {
                I_out[base_idx + count] = (double)dofs[r];
                J_out[base_idx + count] = (double)dofs[c];
                V_out[base_idx + count] = K_local[r][c];
                count++;
            }
        }
    }
}