/**
 * assemble_jacobian_kernel_mex.cpp (Update for V7 Utils)
 * 功能: 适配新的 MexMaterialUtils 构造函数 (接收 explicit coe_start)
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <vector>
#include "MexElemUtils.hpp"
#include "MexMaterialUtils.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Updated Expected Args: 18 inputs (Added M_coe_start at index 15)
    if (nrhs != 18) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 18 inputs.");

    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* Dofs = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* Tags = mxGetPr(prhs[4]);
    double* SolA = mxGetPr(prhs[5]);
    double* Qw = mxGetPr(prhs[6]);
    double* RefCurl = mxGetPr(prhs[7]);

    double* M_lin = mxGetPr(prhs[8]);
    double* M_isNon = mxGetPr(prhs[9]);
    double* M_maxB = mxGetPr(prhs[10]);
    double* M_brk = mxGetPr(prhs[11]);
    double* M_coe = mxGetPr(prhs[12]);
    double* M_start = mxGetPr(prhs[13]);
    double* M_cnt = mxGetPr(prhs[14]);
    double* M_coe_start = mxGetPr(prhs[15]); // NEW input
    
    double  TimeStep = mxGetScalar(prhs[16]); 
    bool    CalcJ = (mxGetScalar(prhs[17]) > 0.5);

    size_t numElems = mxGetN(prhs[1]);
    size_t numQP = mxGetNumberOfElements(prhs[6]);
    int maxTag = (int)mxGetNumberOfElements(prhs[8]) - 1;

    // Use new constructor
    MaterialEvaluator matEval(M_lin, M_isNon, M_maxB, M_brk, M_coe, M_start, M_cnt, M_coe_start, maxTag);

    size_t numJ = CalcJ ? (numElems * 36) : 0;
    size_t numR = numElems * 6;

    plhs[0] = mxCreateDoubleMatrix(numJ, 1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(numJ, 1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(numJ, 1, mxREAL); 
    plhs[3] = mxCreateDoubleMatrix(numR, 1, mxREAL); 
    plhs[4] = mxCreateDoubleMatrix(numR, 1, mxREAL); 

    double* I_out = mxGetPr(plhs[0]);
    double* J_out = mxGetPr(plhs[1]);
    double* V_out = mxGetPr(plhs[2]);
    double* R_idx = mxGetPr(plhs[3]);
    double* R_val = mxGetPr(plhs[4]);

    #pragma omp parallel for
    for (long e = 0; e < numElems; e++) {
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

        double signs[6];
        long dofs[6];
        double a_local[6];

        for(int k=0; k<6; k++) {
            signs[k] = Signs[k + e*6];
            dofs[k] = (long)Dofs[k + e*6];
            long global_idx = dofs[k] - 1;
            if (global_idx >= 0) {
                a_local[k] = SolA[global_idx] * signs[k];
            } else {
                a_local[k] = 0.0;
            }
        }

        double K_local[6][6] = {0};
        double R_local[6]    = {0};
        int tag = (int)Tags[e];

        for (int q = 0; q < numQP; q++) {
            double w = Qw[q] * absDetJ; 
            long q_offset = q * 18;

            double B[3] = {0, 0, 0};
            double CurlNi[6][3]; 

            for (int i = 0; i < 6; i++) {
                double ref_c[3];
                ref_c[0] = RefCurl[i + 0  + q_offset];
                ref_c[1] = RefCurl[i + 6  + q_offset];
                ref_c[2] = RefCurl[i + 12 + q_offset];
                double temp[3];
                J_mat.multVec(ref_c, temp);
                CurlNi[i][0] = temp[0] * invDetJ;
                CurlNi[i][1] = temp[1] * invDetJ;
                CurlNi[i][2] = temp[2] * invDetJ;
                B[0] += a_local[i] * CurlNi[i][0];
                B[1] += a_local[i] * CurlNi[i][1];
                B[2] += a_local[i] * CurlNi[i][2];
            }

            double B_sq = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
            double dnu = 0.0;
            double nu = matEval.evaluate(tag, B_sq, dnu);
            
            for (int i = 0; i < 6; i++) {
                double B_dot_Ci = B[0]*CurlNi[i][0] + B[1]*CurlNi[i][1] + B[2]*CurlNi[i][2];
                R_local[i] += nu * B_dot_Ci * w;
                if (CalcJ) {
                    for (int j = 0; j < 6; j++) {
                        double Ci_dot_Cj = CurlNi[i][0]*CurlNi[j][0] + 
                                           CurlNi[i][1]*CurlNi[j][1] + 
                                           CurlNi[i][2]*CurlNi[j][2];
                        double B_dot_Cj = B[0]*CurlNi[j][0] + B[1]*CurlNi[j][1] + B[2]*CurlNi[j][2];
                        double val = nu * Ci_dot_Cj + 2.0 * dnu * B_dot_Ci * B_dot_Cj;
                        K_local[i][j] += val * w;
                    }
                }
            }
        }

        long r_base = e * 6;
        for (int i = 0; i < 6; i++) {
            R_idx[r_base + i] = (double)dofs[i];
            R_val[r_base + i] = R_local[i] * signs[i];
        }

        if (CalcJ) {
            long j_base = e * 36;
            int count = 0;
            for (int c = 0; c < 6; c++) {
                for (int r = 0; r < 6; r++) {
                    I_out[j_base + count] = (double)dofs[r];
                    J_out[j_base + count] = (double)dofs[c];
                    V_out[j_base + count] = K_local[r][c] * signs[r] * signs[c];
                    count++;
                }
            }
        }
    }
}