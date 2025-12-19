#include "mex.h"
#include <omp.h>
#include <cmath>
#include "MexElemUtils.hpp"
#include "MexMaterialUtils.hpp" 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // [修复点]: 增加到 17 个输入参数
    if (nrhs != 17) mexErrMsgIdAndTxt("MagPackage:Args", "Expected 17 inputs.");

    // 基础输入
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* Dofs = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* Tags = mxGetPr(prhs[4]);
    double* SolA = mxGetPr(prhs[5]);
    double* Qw = mxGetPr(prhs[6]);
    double* RefCurl = mxGetPr(prhs[7]);

    // 材质打包数据
    double* M_lin = mxGetPr(prhs[8]);
    double* M_isNon = mxGetPr(prhs[9]);
    double* M_maxB = mxGetPr(prhs[10]);
    double* M_brk = mxGetPr(prhs[11]);
    double* M_coe = mxGetPr(prhs[12]);
    double* M_start = mxGetPr(prhs[13]);
    double* M_cnt = mxGetPr(prhs[14]);
    double* M_coe_start = mxGetPr(prhs[15]); // [新增参数]

    bool CalcJ = (bool)mxGetScalar(prhs[16]); // [索引后移]

    size_t numElems = mxGetN(prhs[1]);
    size_t n_q = mxGetM(prhs[6]) * mxGetN(prhs[6]);

    MexMaterialData matData;
    matData.linearNu = M_lin;
    matData.isNonlinear = M_isNon;
    matData.maxBSq = M_maxB;
    matData.breaks = M_brk;
    matData.coefs = M_coe;
    matData.splineStart = M_start;
    matData.splineCount = M_cnt;
    matData.splineCoeStart = M_coe_start; // [新增挂载]

    // 输出分配
    size_t nnz = CalcJ ? (numElems * 36) : 0;
    plhs[0] = mxCreateDoubleMatrix(nnz, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(nnz, 1, mxREAL); // J
    plhs[2] = mxCreateDoubleMatrix(nnz, 1, mxREAL); // V
    plhs[3] = mxCreateDoubleMatrix(numElems * 6, 1, mxREAL); // R_idx
    plhs[4] = mxCreateDoubleMatrix(numElems * 6, 1, mxREAL); // R_val

    double* I_out = mxGetPr(plhs[0]);
    double* J_out = mxGetPr(plhs[1]);
    double* V_out = mxGetPr(plhs[2]);
    double* R_idx = mxGetPr(plhs[3]);
    double* R_val = mxGetPr(plhs[4]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < numElems; e++) {
        double nodes[4][3];
        for (int i = 0; i < 4; i++) {
            int nid = (int)T[i + e*4] - 1;
            nodes[i][0] = P[nid*3 + 0];
            nodes[i][1] = P[nid*3 + 1];
            nodes[i][2] = P[nid*3 + 2];
        }

        Mat3x3 J_mat;
        double detJ = MexElemUtils::compute_jacobian_3d(&nodes[0][0], J_mat);
        double absDetJ = std::abs(detJ);
        double invDetJ = 1.0 / detJ;

        int dofs[6];
        double signs[6];
        double A_local[6];
        for (int i = 0; i < 6; i++) {
            dofs[i] = (int)Dofs[i + e*6];
            signs[i] = Signs[i + e*6];
            int gid = dofs[i] - 1; 
            if (gid >= 0) A_local[i] = SolA[gid] * signs[i];
            else          A_local[i] = 0.0;
        }

        double CurlNi[6][3];
        double R_local[6] = {0};
        double K_local[6][6] = {0};
        int tag = (int)Tags[e];

        for (size_t q = 0; q < n_q; q++) {
            double w = Qw[q] * absDetJ;
            for (int i = 0; i < 6; i++) {
                double c_ref[3];
                c_ref[0] = RefCurl[i + 0*6 + q*18];
                c_ref[1] = RefCurl[i + 1*6 + q*18];
                c_ref[2] = RefCurl[i + 2*6 + q*18];

                double c_phy[3];
                J_mat.multVec(c_ref, c_phy);
                CurlNi[i][0] = c_phy[0] * invDetJ;
                CurlNi[i][1] = c_phy[1] * invDetJ;
                CurlNi[i][2] = c_phy[2] * invDetJ;
            }

            double B[3] = {0, 0, 0};
            for (int i = 0; i < 6; i++) {
                B[0] += A_local[i] * CurlNi[i][0];
                B[1] += A_local[i] * CurlNi[i][1];
                B[2] += A_local[i] * CurlNi[i][2];
            }
            double B_sq = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];

            double nu, dnu;
            evaluate_nu_derivative(B_sq, tag, matData, nu, dnu);
            
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
            long k_base = e * 36;
            int idx = 0;
            for (int j = 0; j < 6; j++) {
                for (int i = 0; i < 6; i++) {
                    I_out[k_base + idx] = (double)dofs[i];
                    J_out[k_base + idx] = (double)dofs[j];
                    V_out[k_base + idx] = K_local[i][j] * signs[i] * signs[j];
                    idx++;
                }
            }
        }
    }
}