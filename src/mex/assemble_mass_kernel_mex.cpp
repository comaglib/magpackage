#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>

// 辅助：计算 Jacobian 和 inv(J)'
// N_phy = N_ref * iJt'
void compute_jacobian_and_iJt_trans(const double* P, const double* T_col, 
                                    double& detJ, double* iJt_trans) {
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

    // detJ
    detJ = d1x*(d2y*d3z - d2z*d3y) - d1y*(d2x*d3z - d2z*d3x) + d1z*(d2x*d3y - d2y*d3x);
    double invDetJ = 1.0 / detJ;

    // iJt (inv(J)^T) 的列是 J^{-T} 的列 = grad_xi / detJ
    // cross(d2, d3)
    iJt_trans[0] = (d2y*d3z - d2z*d3y) * invDetJ;
    iJt_trans[1] = (d2z*d3x - d2x*d3z) * invDetJ;
    iJt_trans[2] = (d2x*d3y - d2y*d3x) * invDetJ;

    // cross(d3, d1)
    iJt_trans[3] = (d3y*d1z - d3z*d1y) * invDetJ;
    iJt_trans[4] = (d3z*d1x - d3x*d1z) * invDetJ;
    iJt_trans[5] = (d3x*d1y - d3y*d1x) * invDetJ;

    // cross(d1, d2)
    iJt_trans[6] = (d1y*d2z - d1z*d2y) * invDetJ;
    iJt_trans[7] = (d1z*d2x - d1x*d2z) * invDetJ;
    iJt_trans[8] = (d1x*d2y - d1y*d2x) * invDetJ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Inputs: P, T, Dofs, Signs, Coeff, Qw, ValRef
    if (nrhs < 7) mexErrMsgIdAndTxt("MagPackage:Mass", "Need 7 inputs.");

    double *P = mxGetPr(prhs[0]);
    double *T = mxGetPr(prhs[1]); 
    size_t n_elems = mxGetN(prhs[1]);
    
    double *Dofs = mxGetPr(prhs[2]);
    double *Signs = mxGetPr(prhs[3]);
    
    // Coeff 可能为空
    double *Coeff = nullptr;
    bool has_coeff = false;
    if (!mxIsEmpty(prhs[4])) {
        Coeff = mxGetPr(prhs[4]);
        has_coeff = true;
    }
    
    double *Qw = mxGetPr(prhs[5]);
    size_t n_q = mxGetM(prhs[5]) * mxGetN(prhs[5]);
    
    double *ValRef = mxGetPr(prhs[6]); // [6 x 3 x Nq]

    // Output: I, J, V
    size_t n_triplets = n_elems * 36;
    plhs[0] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL);
    
    double *I_out = mxGetPr(plhs[0]);
    double *J_out = mxGetPr(plhs[1]);
    double *V_out = mxGetPr(plhs[2]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        // Coeff check
        double c_val = has_coeff ? Coeff[e] : 1.0;
        if (c_val == 0.0) {
            // Fill with zeros (or handle sparsity properly, here we just fill zeros to keep index alignment simple)
            size_t offset = e * 36;
            // Need to fill I, J with valid dofs to avoid sparse error, but V=0
            double dofs[6];
            for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
            int idx = 0;
            for (int c = 0; c < 6; ++c) {
                for (int r = 0; r < 6; ++r) {
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
        double iJt_trans[9]; // inv(J)'
        compute_jacobian_and_iJt_trans(P, T_local, detJ, iJt_trans);
        double absDetJ = std::abs(detJ);

        double Me[36] = {0};

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ * c_val;
            
            // ValRef: [6x3] at q
            // N_phy = N_ref * iJt_trans' 
            // N_ref row i: [nr_x, nr_y, nr_z]
            // iJt_trans = [j1x j1y j1z; j2x ...; j3x ...] (row major storage of matrix)
            // But wait, iJt_trans calculated above is actually inv(J)' flattened.
            // Let M = iJt_trans. 
            // N_phy = N_ref * M
            
            double N_phy[18]; 

            for (int i = 0; i < 6; ++i) {
                double nr_x = ValRef[i + 0*6 + q*18];
                double nr_y = ValRef[i + 1*6 + q*18];
                double nr_z = ValRef[i + 2*6 + q*18];
                
                // M is 3x3. M_00=iJt_trans[0], M_01=iJt_trans[1]...
                // N_phy_x = nr_x * M00 + nr_y * M10 + nr_z * M20
                // iJt_trans in my code: 
                // [0,1,2] is col 1 of inv(J)^T (which is row 1 of inv(J))
                // Code above computes iJt_trans such that it IS inv(J)^T.
                // So N_phy = N_ref * iJt_trans
                
                N_phy[i+0*6] = nr_x * iJt_trans[0] + nr_y * iJt_trans[3] + nr_z * iJt_trans[6];
                N_phy[i+1*6] = nr_x * iJt_trans[1] + nr_y * iJt_trans[4] + nr_z * iJt_trans[7];
                N_phy[i+2*6] = nr_x * iJt_trans[2] + nr_y * iJt_trans[5] + nr_z * iJt_trans[8];
            }
            
            // Me += w * N_phy * N_phy'
            for (int c = 0; c < 6; ++c) {
                for (int r = 0; r < 6; ++r) {
                    double dot_val = N_phy[r + 0*6] * N_phy[c + 0*6] + 
                                     N_phy[r + 1*6] * N_phy[c + 1*6] + 
                                     N_phy[r + 2*6] * N_phy[c + 2*6];
                    Me[r + c*6] += w * dot_val;
                }
            }
        }
        
        // Fill
        double s[6];
        for(int k=0; k<6; k++) s[k] = Signs[k + e*6];
        double dofs[6];
        for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
        
        size_t offset = e * 36;
        int idx = 0;
        for (int c = 0; c < 6; ++c) {
            for (int r = 0; r < 6; ++r) {
                V_out[offset + idx] = Me[r + c*6] * s[r] * s[c];
                I_out[offset + idx] = dofs[r];
                J_out[offset + idx] = dofs[c];
                idx++;
            }
        }
    }
}