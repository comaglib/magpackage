#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>
#include "MexMaterialUtils.hpp"

// 复用 Jacobian 计算
void compute_jacobian_utils(const double* P, const double* T_col, 
                           double* J_mat, double& detJ, double& invDetJ) {
    int n1 = (int)T_col[0] - 1; int n2 = (int)T_col[1] - 1;
    int n3 = (int)T_col[2] - 1; int n4 = (int)T_col[3] - 1;
    double x1 = P[3*n1], y1 = P[3*n1+1], z1 = P[3*n1+2];
    double x2 = P[3*n2], y2 = P[3*n2+1], z2 = P[3*n2+2];
    double x3 = P[3*n3], y3 = P[3*n3+1], z3 = P[3*n3+2];
    double x4 = P[3*n4], y4 = P[3*n4+1], z4 = P[3*n4+2];
    
    J_mat[0] = x2-x1; J_mat[1] = x3-x1; J_mat[2] = x4-x1;
    J_mat[3] = y2-y1; J_mat[4] = y3-y1; J_mat[5] = y4-y1;
    J_mat[6] = z2-z1; J_mat[7] = z3-z1; J_mat[8] = z4-z1;

    detJ = J_mat[0]*(J_mat[4]*J_mat[8] - J_mat[5]*J_mat[7]) - 
           J_mat[1]*(J_mat[3]*J_mat[8] - J_mat[5]*J_mat[6]) + 
           J_mat[2]*(J_mat[3]*J_mat[7] - J_mat[4]*J_mat[6]);
    invDetJ = 1.0 / detJ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Inputs: 
    // 0: P, 1: T, 2: Dofs, 3: Signs, 4: Tags
    // 5: SolutionA (Global Vector)
    // 6: Qw, 7: CurlRef
    // 8-15: Material Data Arrays (LinearNu, IsNonlinear, MaxBSq, Breaks, Coefs, Starts, Counts, NumTags)
    // 16: RequestJacobian (scalar)
    
    if (nrhs < 17) mexErrMsgIdAndTxt("MagPackage:Jacobian", "Need 17 inputs.");

    double *P = mxGetPr(prhs[0]);
    double *T = mxGetPr(prhs[1]); size_t n_elems = mxGetN(prhs[1]);
    double *Dofs = mxGetPr(prhs[2]);
    double *Signs = mxGetPr(prhs[3]);
    double *Tags = mxGetPr(prhs[4]);
    double *SolA = mxGetPr(prhs[5]); // Global solution
    
    double *Qw = mxGetPr(prhs[6]); size_t n_q = mxGetM(prhs[6]) * mxGetN(prhs[6]);
    double *CurlRef = mxGetPr(prhs[7]);
    
    // Construct Material Helper
    MexMaterialData matData;
    matData.linearNu    = mxGetPr(prhs[8]);
    matData.isNonlinear = mxGetPr(prhs[9]);
    matData.maxBSq      = mxGetPr(prhs[10]);
    matData.breaks      = mxGetPr(prhs[11]);
    matData.coefs       = mxGetPr(prhs[12]);
    matData.splineStart = mxGetPr(prhs[13]);
    matData.splineCount = mxGetPr(prhs[14]);
    // input 15 is numTags, unused directly
    
    bool requestJ = (bool)mxGetScalar(prhs[16]);

    // Alloc Outputs
    // K: 36 entries, R: 6 entries
    size_t n_K = requestJ ? (n_elems * 36) : 0;
    size_t n_R = n_elems * 6;
    
    plhs[0] = mxCreateDoubleMatrix(n_K, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(n_K, 1, mxREAL); // J
    plhs[2] = mxCreateDoubleMatrix(n_K, 1, mxREAL); // V
    plhs[3] = mxCreateDoubleMatrix(n_R, 1, mxREAL); // R_idx
    plhs[4] = mxCreateDoubleMatrix(n_R, 1, mxREAL); // R_val

    double *I_out = mxGetPr(plhs[0]);
    double *J_out = mxGetPr(plhs[1]);
    double *V_out = mxGetPr(plhs[2]);
    double *R_idx = mxGetPr(plhs[3]);
    double *R_val = mxGetPr(plhs[4]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        // 1. Get Element Data
        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        
        double dofs[6], s[6], A_local[6];
        for(int k=0; k<6; k++) {
            dofs[k] = Dofs[k + e*6];
            s[k] = Signs[k + e*6];
            // Gather solution: A_local = s * A_global
            int global_id = (int)dofs[k] - 1; // 0-based
            A_local[k] = SolA[global_id] * s[k];
        }
        
        int tag = (int)Tags[e];
        
        // 2. Jacobian
        double J_mat[9], detJ, invDetJ;
        compute_jacobian_utils(P, T_local, J_mat, detJ, invDetJ);
        double absDetJ = std::abs(detJ);

        double Ke[36] = {0};
        double Re[6] = {0};
        
        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ;
            
            // Calc B = Curl(A)
            // C_phy = (C_ref * J_mat') * invDetJ
            // C_phy is 6x3. We compute B = C_phy' * A_local (3x1 vector)
            double B[3] = {0};
            double C_phy[18]; // Store for J calculation
            
            for (int i = 0; i < 6; ++i) {
                // Fetch CurlRef (row i)
                double cr_x = CurlRef[i + 0*6 + q*18];
                double cr_y = CurlRef[i + 1*6 + q*18];
                double cr_z = CurlRef[i + 2*6 + q*18];
                
                // Piola
                double cp_x = (cr_x * J_mat[0] + cr_y * J_mat[1] + cr_z * J_mat[2]) * invDetJ;
                double cp_y = (cr_x * J_mat[3] + cr_y * J_mat[4] + cr_z * J_mat[5]) * invDetJ;
                double cp_z = (cr_x * J_mat[6] + cr_y * J_mat[7] + cr_z * J_mat[8]) * invDetJ;
                
                C_phy[i+0*6] = cp_x;
                C_phy[i+1*6] = cp_y;
                C_phy[i+2*6] = cp_z;
                
                double val = A_local[i];
                B[0] += cp_x * val;
                B[1] += cp_y * val;
                B[2] += cp_z * val;
            }
            
            double B_sq = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
            
            // Material Eval
            double nu, dnu_db2;
            evaluate_nu_derivative(B_sq, tag, matData, nu, dnu_db2);
            
            // Residual contribution: w * (C_phy * (nu * B))
            // H = nu * B
            double Hx = nu * B[0];
            double Hy = nu * B[1];
            double Hz = nu * B[2];
            
            for (int i=0; i<6; ++i) {
                Re[i] += w * (C_phy[i+0*6]*Hx + C_phy[i+1*6]*Hy + C_phy[i+2*6]*Hz);
            }
            
            if (requestJ) {
                // Jacobian: w * C_phy * M_diff * C_phy'
                // M_diff = nu*I + 2*dnu * (B * B')
                // Let's compute Temp = C_phy * M_diff
                // M_diff is symmetric 3x3
                
                double two_dnu = 2.0 * dnu_db2;
                
                // Temp_i (row i of C*M)
                // Temp_i = nu * C_i + 2*dnu * (C_i . B) * B^T
                
                for (int i=0; i<6; ++i) {
                    double ci_x = C_phy[i+0*6];
                    double ci_y = C_phy[i+1*6];
                    double ci_z = C_phy[i+2*6];
                    
                    double dot_CB = ci_x*B[0] + ci_y*B[1] + ci_z*B[2];
                    double scalar = two_dnu * dot_CB;
                    
                    // Temp_vector = nu*Ci + scalar*B
                    double tx = nu*ci_x + scalar*B[0];
                    double ty = nu*ci_y + scalar*B[1];
                    double tz = nu*ci_z + scalar*B[2];
                    
                    // Ke_ij += w * (Temp_i . C_j)
                    for (int j=0; j<6; ++j) {
                        double cj_x = C_phy[j+0*6];
                        double cj_y = C_phy[j+1*6];
                        double cj_z = C_phy[j+2*6];
                        
                        Ke[i + j*6] += w * (tx*cj_x + ty*cj_y + tz*cj_z);
                    }
                }
            }
        }
        
        // Store Re
        size_t off_R = e * 6;
        for(int k=0; k<6; k++) {
            R_idx[off_R + k] = dofs[k];
            R_val[off_R + k] = Re[k] * s[k]; // Apply sign
        }
        
        // Store Ke
        if (requestJ) {
            size_t off_K = e * 36;
            int idx = 0;
            for (int c = 0; c < 6; ++c) {
                for (int r = 0; r < 6; ++r) {
                    V_out[off_K + idx] = Ke[r + c*6] * s[r] * s[c]; // Signs
                    I_out[off_K + idx] = dofs[r];
                    J_out[off_K + idx] = dofs[c];
                    idx++;
                }
            }
        }
    }
}