/**
 * assemble_coupling_kernel_mex.cpp
 * 功能: 组装 A-V 耦合矩阵 C = (sigma * N_i, grad phi_j)
 * 输入: P, T, DofsE(Edge), Signs, DofsN(Node), Coeffs, Qw, NedRef, GradRef
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include "MexElemUtils.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 9) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 9 inputs.");

    // Inputs
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* DofsE = mxGetPr(prhs[2]); // Edge DoFs (Rows, 6 per elem)
    double* Signs = mxGetPr(prhs[3]); // Edge Signs
    double* DofsN = mxGetPr(prhs[4]); // Node DoFs (Cols, 4 per elem)
    double* Coeffs = mxGetPr(prhs[5]); // Sigma
    double* Qw = mxGetPr(prhs[6]);
    double* NedRef = mxGetPr(prhs[7]);  // [6 x 3 x Nq]
    double* GradRef = mxGetPr(prhs[8]); // [4 x 3 x Nq] (Usually constant for P1, but support Nq)

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[6]);
    
    // Output: Sparse Triplet (I, J, V)
    // Local matrix is 6 (Edge) x 4 (Node) = 24 entries per elem
    size_t numTriplets = numElems * 24;
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

        double C_local[6][4] = {0}; // 6 Edges x 4 Nodes
        double sigma = Coeffs[e];
        
        if (std::abs(sigma) > 1e-12) {
            for (int q = 0; q < numQP; q++) {
                double w = Qw[q] * absDetJ;
                long q_offset = q * 18; // 6*3
                long q_offset_grad = q * 12; // 4*3

                // Transform Nedelec Basis (Edge) -> N_phy = J^{-T} * N_ref
                double N_phy[3][6]; 
                for (int i = 0; i < 6; i++) {
                    double ref_n[3];
                    ref_n[0] = NedRef[i + 0  + q_offset];
                    ref_n[1] = NedRef[i + 6  + q_offset];
                    ref_n[2] = NedRef[i + 12 + q_offset];
                    
                    double temp[3];
                    J_inv.multVecTrans(ref_n, temp);
                    N_phy[0][i] = temp[0];
                    N_phy[1][i] = temp[1];
                    N_phy[2][i] = temp[2];
                }

                // Transform Gradient Basis (Node) -> Grad_phy = J^{-T} * Grad_ref
                double G_phy[3][4];
                for (int i = 0; i < 4; i++) {
                    double ref_g[3];
                    // GradRef Layout: [4 rows x 3 cols x Nq]
                    ref_g[0] = GradRef[i + 0 + q_offset_grad];
                    ref_g[1] = GradRef[i + 4 + q_offset_grad];
                    ref_g[2] = GradRef[i + 8 + q_offset_grad];

                    double temp[3];
                    J_inv.multVecTrans(ref_g, temp);
                    G_phy[0][i] = temp[0];
                    G_phy[1][i] = temp[1];
                    G_phy[2][i] = temp[2];
                }
                
                // Assembly: (sigma * N_i) . (Grad phi_j)
                for (int i = 0; i < 6; i++) {   // Edge Row
                    for (int j = 0; j < 4; j++) { // Node Col
                        double dot = N_phy[0][i]*G_phy[0][j] + 
                                     N_phy[1][i]*G_phy[1][j] + 
                                     N_phy[2][i]*G_phy[2][j];
                        C_local[i][j] += dot * sigma * w;
                    }
                }
            }
        }

        double signs[6];
        for(int k=0; k<6; k++) signs[k] = Signs[k + e*6];
        
        long row_dofs[6];
        for(int k=0; k<6; k++) row_dofs[k] = (long)DofsE[k + e*6];
        
        long col_dofs[4];
        for(int k=0; k<4; k++) col_dofs[k] = (long)DofsN[k + e*4];

        long base_idx = e * 24;
        int count = 0;
        for (int c = 0; c < 4; c++) {     // Node Cols
            for (int r = 0; r < 6; r++) { // Edge Rows
                I_out[base_idx + count] = (double)row_dofs[r];
                J_out[base_idx + count] = (double)col_dofs[c];
                
                // Sign only applies to Edge functions
                V_out[base_idx + count] = C_local[r][c] * signs[r];
                count++;
            }
        }
    }
}