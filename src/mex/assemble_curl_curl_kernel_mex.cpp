/**
 * assemble_curl_curl_kernel_mex.cpp (Final Fixed Version)
 * 功能: 组装 Nedelec 单元的旋度-旋度 (Curl-Curl) 矩阵
 * * 修复日志:
 * 1. [Critical Fix] 修正 RefCurl 索引逻辑。
 * MATLAB 输入为 [6 x 3 x Nq] (列优先: Basis -> Component -> QP)。
 * 原代码错误地假设为 [3 x 6]，导致数据读取错位。
 * 2. [Verified] 移除了错误的 1/6 权重因子。
 * 3. [Compat] 保持 mxGetPr 以兼容 MinGW。
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
    double* Nu = mxGetPr(prhs[4]);
    double* Qw = mxGetPr(prhs[5]);
    double* RefCurl = mxGetPr(prhs[6]); // [6 x 3 x NumQP]

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[5]);
    
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
        
        // 1. Get Geometry
        double p_local[12];
        for (int n = 0; n < 4; n++) {
            long node_idx = (long)T[n + e*4] - 1; 
            p_local[n*3 + 0] = P[node_idx*3 + 0];
            p_local[n*3 + 1] = P[node_idx*3 + 1];
            p_local[n*3 + 2] = P[node_idx*3 + 2];
        }

        Mat3x3 J_mat;
        double detJ = MexElemUtils::compute_jacobian_3d(p_local, J_mat);
        double absDetJ = std::abs(detJ);
        double invDetJ = 1.0 / detJ;

        double K_local[6][6] = {0};
        double nu_val = Nu[e];
        
        // 2. Integration Loop
        for (int q = 0; q < numQP; q++) {
            double w = Qw[q] * absDetJ; // Correct weight (no / 6.0)
            
            double B_real[3][6]; 
            
            // RefCurl Memory Layout: [6 rows x 3 cols x NumQP]
            // Column-Major Index: i + 6*c + 18*q
            long q_offset = q * 18;
            
            for (int i = 0; i < 6; i++) {
                double ref_c[3];
                // [Corrected Indexing]
                // Col 0 (x): i + 0*6
                // Col 1 (y): i + 1*6
                // Col 2 (z): i + 2*6
                ref_c[0] = RefCurl[i + 0  + q_offset];
                ref_c[1] = RefCurl[i + 6  + q_offset];
                ref_c[2] = RefCurl[i + 12 + q_offset];
                
                // Piola Transform: (1/detJ) * J * curl_ref
                double temp[3];
                J_mat.multVec(ref_c, temp);
                
                B_real[0][i] = temp[0] * invDetJ;
                B_real[1][i] = temp[1] * invDetJ;
                B_real[2][i] = temp[2] * invDetJ;
            }

            // Local Stiffness Assembly
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    double dot_val = B_real[0][i]*B_real[0][j] + 
                                     B_real[1][i]*B_real[1][j] + 
                                     B_real[2][i]*B_real[2][j];
                    K_local[i][j] += dot_val * nu_val * w;
                }
            }
        }

        // 3. Fill Output
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
                V_out[base_idx + count] = K_local[r][c] * signs[r] * signs[c];
                count++;
            }
        }
    }
}