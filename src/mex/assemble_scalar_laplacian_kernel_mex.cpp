#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>

// 复用 Jacobian 计算逻辑
void compute_jacobian_and_invJ(const double* P, const double* T_col, 
                               double& detJ, double* invJ) {
    int n1 = (int)T_col[0] - 1;
    int n2 = (int)T_col[1] - 1;
    int n3 = (int)T_col[2] - 1;
    int n4 = (int)T_col[3] - 1;

    double x1 = P[3*n1], y1 = P[3*n1+1], z1 = P[3*n1+2];
    double x2 = P[3*n2], y2 = P[3*n2+1], z2 = P[3*n2+2];
    double x3 = P[3*n3], y3 = P[3*n3+1], z3 = P[3*n3+2];
    double x4 = P[3*n4], y4 = P[3*n4+1], z4 = P[3*n4+2];

    double d1x = x2-x1, d1y = y2-y1, d1z = z2-z1;
    double d2x = x3-x1, d2y = y3-y1, d2z = z3-z1;
    double d3x = x4-x1, d3y = y4-y1, d3z = z4-z1;

    detJ = d1x*(d2y*d3z - d2z*d3y) - d1y*(d2x*d3z - d2z*d3x) + d1z*(d2x*d3y - d2y*d3x);
    double invDetJ = 1.0 / detJ;

    invJ[0] = (d2y*d3z - d2z*d3y) * invDetJ;
    invJ[1] = (d2z*d3x - d2x*d3z) * invDetJ;
    invJ[2] = (d2x*d3y - d2y*d3x) * invDetJ;
    invJ[3] = (d3y*d1z - d3z*d1y) * invDetJ;
    invJ[4] = (d3z*d1x - d3x*d1z) * invDetJ;
    invJ[5] = (d3x*d1y - d3y*d1x) * invDetJ;
    invJ[6] = (d1y*d2z - d1z*d2y) * invDetJ;
    invJ[7] = (d1z*d2x - d1x*d2z) * invDetJ;
    invJ[8] = (d1x*d2y - d1y*d2x) * invDetJ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Inputs: P, T, Dofs, Coeff, Qw, GradRef
    // Note: Lagrange P1 has 4 nodes per element. Dofs is [4 x Ne].
    if (nrhs < 6) mexErrMsgIdAndTxt("MagPackage:ScalarLaplacian", "Need 6 inputs.");

    double *P = mxGetPr(prhs[0]);
    double *T = mxGetPr(prhs[1]); 
    size_t n_elems = mxGetN(prhs[1]);
    
    double *Dofs = mxGetPr(prhs[2]);
    
    // Coeff can be empty or scalar or vector
    double *Coeff = nullptr;
    bool has_coeff = false;
    if (!mxIsEmpty(prhs[3])) {
        Coeff = mxGetPr(prhs[3]);
        has_coeff = true;
    }
    
    double *Qw = mxGetPr(prhs[4]);
    size_t n_q = mxGetM(prhs[4]) * mxGetN(prhs[4]);
    
    double *GradRef = mxGetPr(prhs[5]); // [4 x 3 x Nq]

    // Output: I, J, V (16 entries per elem)
    size_t n_triplets = n_elems * 16;
    plhs[0] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL);
    
    double *I_out = mxGetPr(plhs[0]);
    double *J_out = mxGetPr(plhs[1]);
    double *V_out = mxGetPr(plhs[2]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        double c_val = has_coeff ? Coeff[e] : 1.0;
        if (c_val == 0.0) {
            // Fill zeros
            size_t offset = e * 16;
            double dofs[4];
            for(int k=0; k<4; k++) dofs[k] = Dofs[k + e*4];
            int idx = 0;
            for (int c = 0; c < 4; ++c) {
                for (int r = 0; r < 4; ++r) {
                    I_out[offset + idx] = dofs[r];
                    J_out[offset + idx] = dofs[c];
                    V_out[offset + idx] = 0.0;
                    idx++;
                }
            }
            continue;
        }

        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        
        double detJ;
        double invJ[9]; 
        compute_jacobian_and_invJ(P, T_local, detJ, invJ);
        double absDetJ = std::abs(detJ);

        double Ke[16] = {0};

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ * c_val;
            
            double G_phy[12]; // 4 nodes x 3 coords
            
            // Grad_phy = Grad_ref * J^{-1}
            // G_ref row i: [gr_x, gr_y, gr_z]
            // invJ is Row-Major J^{-1}
            for (int i = 0; i < 4; ++i) {
                double gr_x = GradRef[i + 0*4 + q*12];
                double gr_y = GradRef[i + 1*4 + q*12];
                double gr_z = GradRef[i + 2*4 + q*12];
                
                // Row vector mul: [gx, gy, gz] * [invJ]
                G_phy[i+0*4] = gr_x * invJ[0] + gr_y * invJ[3] + gr_z * invJ[6];
                G_phy[i+1*4] = gr_x * invJ[1] + gr_y * invJ[4] + gr_z * invJ[7];
                G_phy[i+2*4] = gr_x * invJ[2] + gr_y * invJ[5] + gr_z * invJ[8];
            }
            
            // Ke += w * (G_phy * G_phy')
            for (int c = 0; c < 4; ++c) {
                for (int r = 0; r < 4; ++r) {
                    double dot_val = G_phy[r + 0*4] * G_phy[c + 0*4] + 
                                     G_phy[r + 1*4] * G_phy[c + 1*4] + 
                                     G_phy[r + 2*4] * G_phy[c + 2*4];
                    Ke[r + c*4] += w * dot_val;
                }
            }
        }
        
        double dofs[4];
        for(int k=0; k<4; k++) dofs[k] = Dofs[k + e*4];
        
        size_t offset = e * 16;
        int idx = 0;
        for (int c = 0; c < 4; ++c) {
            for (int r = 0; r < 4; ++r) {
                V_out[offset + idx] = Ke[r + c*4]; // No signs for scalar Lagrange
                I_out[offset + idx] = dofs[r];
                J_out[offset + idx] = dofs[c];
                idx++;
            }
        }
    }
}