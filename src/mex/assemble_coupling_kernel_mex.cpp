/**
 * assemble_coupling_kernel_mex.cpp
 * * 功能: 组装 A-V 耦合矩阵
 * 公式: C_ij = \int sigma * (N_i . grad(phi_j)) dV
 * * 输入:
 * 0: P        (3 x NumNodes)
 * 1: T        (4 x NumElems)
 * 2: DofsE    (6 x NumElems) - Nedelec DoFs
 * 3: Signs    (6 x NumElems)
 * 4: DofsN    (4 x NumElems) - Lagrange DoFs
 * 5: Coeffs   (NumElems x 1) - Sigma
 * 6: Qw       (n_q x 1)
 * 7: NedRef   (6 x 3 x n_q)
 * 8: GradRef  (4 x 3 x n_q) - Lagrange Gradients
 * * 输出:
 * 0: I (Row Indices - Edge)
 * 1: J (Col Indices - Node)
 * 2: V (Values)
 */

#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>

void compute_jacobian_and_invJ(const double* P, const double* T_col, 
                               double& detJ, double* invJ) {
    int n1 = (int)T_col[0] - 1; int n2 = (int)T_col[1] - 1;
    int n3 = (int)T_col[2] - 1; int n4 = (int)T_col[3] - 1;

    double x1 = P[3*n1], y1 = P[3*n1+1], z1 = P[3*n1+2];
    double x2 = P[3*n2], y2 = P[3*n2+1], z2 = P[3*n2+2];
    double x3 = P[3*n3], y3 = P[3*n3+1], z3 = P[3*n3+2];
    double x4 = P[3*n4], y4 = P[3*n4+1], z4 = P[3*n4+2];

    double d1x = x2-x1, d1y = y2-y1, d1z = z2-z1;
    double d2x = x3-x1, d2y = y3-y1, d2z = z3-z1;
    double d3x = x4-x1, d3y = y4-y1, d3z = z4-z1;

    detJ = d1x*(d2y*d3z - d2z*d3y) - d1y*(d2x*d3z - d2z*d3x) + d1z*(d2x*d3y - d2y*d3x);
    double invDet = 1.0 / detJ;

    invJ[0] = (d2y*d3z - d2z*d3y) * invDet; invJ[3] = (d3y*d1z - d3z*d1y) * invDet; invJ[6] = (d2x*d3y - d2y*d3x) * invDet;
    invJ[1] = (d3z*d1y - d1z*d3y) * invDet; invJ[4] = (d1x*d3z - d1z*d3x) * invDet; invJ[7] = (d3x*d1y - d1x*d3y) * invDet;
    invJ[2] = (d1z*d2y - d1y*d2z) * invDet; invJ[5] = (d1y*d2x - d1x*d2y) * invDet; invJ[8] = (d1x*d2z - d1z*d2x) * invDet;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 9) mexErrMsgTxt("Args: P, T, DofsE, Signs, DofsN, Coeffs, Qw, NedRef, GradRef");

    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* DofsE = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* DofsN = mxGetPr(prhs[4]);
    double* Coeffs = mxGetPr(prhs[5]);
    double* Qw = mxGetPr(prhs[6]);
    double* NedRef = mxGetPr(prhs[7]);
    double* GradRef = mxGetPr(prhs[8]);

    size_t numElems = mxGetN(prhs[1]);
    size_t n_q = mxGetM(prhs[6]);

    // Output: 6 edges * 4 nodes = 24 entries per element
    plhs[0] = mxCreateDoubleMatrix(24 * numElems, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(24 * numElems, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(24 * numElems, 1, mxREAL);

    double* I = mxGetPr(plhs[0]);
    double* J = mxGetPr(plhs[1]);
    double* V = mxGetPr(plhs[2]);

    #pragma omp parallel for
    for (size_t e = 0; e < numElems; ++e) {
        double sigma = Coeffs[e];
        if (std::abs(sigma) < 1e-15) {
             // Fill zeros
             double* de = &DofsE[e*6];
             double* dn = &DofsN[e*4];
             for(int k=0; k<24; ++k) {
                 I[e*24 + k] = de[k % 6];
                 J[e*24 + k] = dn[k / 6];
                 V[e*24 + k] = 0.0;
             }
             continue;
        }

        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        double detJ, invJ[9]; 
        compute_jacobian_and_invJ(P, T_local, detJ, invJ);
        double absDetJ = std::abs(detJ);

        // Local Coupling Matrix: [6 x 4] (Edge x Node)
        double Ce[24] = {0}; 

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ * sigma;
            size_t ned_off = q * 18; // 6x3
            size_t grad_off = q * 12; // 4x3

            // Precompute Physical Gradients (GradPhi)
            // G_phy = GradRef * invJ [4x3]
            double G_phy[12]; 
            for(int j=0; j<4; ++j) {
                double gr_x = GradRef[j + 0*4 + grad_off];
                double gr_y = GradRef[j + 1*4 + grad_off];
                double gr_z = GradRef[j + 2*4 + grad_off];
                
                G_phy[j + 0*4] = gr_x*invJ[0] + gr_y*invJ[3] + gr_z*invJ[6];
                G_phy[j + 1*4] = gr_x*invJ[1] + gr_y*invJ[4] + gr_z*invJ[7];
                G_phy[j + 2*4] = gr_x*invJ[2] + gr_y*invJ[5] + gr_z*invJ[8];
            }

            // Precompute Physical Nedelec (N_phy)
            // N_phy = NedRef * invJ [6x3]
            double N_phy[18];
            for(int i=0; i<6; ++i) {
                double nr_x = NedRef[i + 0*6 + ned_off];
                double nr_y = NedRef[i + 1*6 + ned_off];
                double nr_z = NedRef[i + 2*6 + ned_off];
                
                N_phy[i + 0*6] = nr_x*invJ[0] + nr_y*invJ[3] + nr_z*invJ[6];
                N_phy[i + 1*6] = nr_x*invJ[1] + nr_y*invJ[4] + nr_z*invJ[7];
                N_phy[i + 2*6] = nr_x*invJ[2] + nr_y*invJ[5] + nr_z*invJ[8];
            }

            // Accumulate Ce += w * (N_phy * G_phy')
            // Ce[row_edge, col_node]
            // N_phy[edge, k] * G_phy[node, k]
            for (int col = 0; col < 4; ++col) {     // Node j
                for (int row = 0; row < 6; ++row) { // Edge i
                    double dot = N_phy[row + 0*6] * G_phy[col + 0*4] + 
                                 N_phy[row + 1*6] * G_phy[col + 1*4] + 
                                 N_phy[row + 2*6] * G_phy[col + 2*4];
                    
                    // Col-major storage for 6x4 matrix: idx = row + col*6
                    Ce[row + col*6] += w * dot;
                }
            }
        }

        // Apply Signs and Store
        double* s = &Signs[e*6];
        double* de = &DofsE[e*6];
        double* dn = &DofsN[e*4];
        
        for (int col = 0; col < 4; ++col) {
            for (int row = 0; row < 6; ++row) {
                size_t global_idx = e*24 + (row + col*6);
                I[global_idx] = de[row];
                J[global_idx] = dn[col];
                V[global_idx] = Ce[row + col*6] * s[row];
            }
        }
    }
}