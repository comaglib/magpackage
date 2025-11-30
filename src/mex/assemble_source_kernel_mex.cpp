/**
 * assemble_source_kernel_mex.cpp
 * * 功能: 高性能组装电流源向量 (Source Vector Assembly)
 * 公式: F_i = \int (N_i \cdot J_s) dV
 * * 输入:
 * 0: P        (3 x NumNodes)
 * 1: T        (4 x NumElems)
 * 2: CellDofs (6 x NumElems)
 * 3: Signs    (6 x NumElems)
 * 4: Js       (3 x NumElems) - 单元电流密度向量
 * 5: Qw       (NumQuadPoints x 1)
 * 6: ValRef   (6 x 3 x NumQuadPoints) - Nedelec 基函数参考值
 * * 输出:
 * 0: I (Indices)
 * 1: V (Values)
 */

#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>

// 复用 Jacobian 计算逻辑 (Copy from assemble_winding_kernel_mex.cpp)
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

    // J = [d1x d2x d3x; d1y d2y d3y; d1z d2z d3z]
    // detJ = d1 . (d2 x d3)
    detJ = d1x*(d2y*d3z - d2z*d3y) - 
           d1y*(d2x*d3z - d2z*d3x) + 
           d1z*(d2x*d3y - d2y*d3x);

    double invDet = 1.0 / detJ;

    // InvJ = (1/det) * Adj(J)
    // Row 0
    invJ[0] = (d2y*d3z - d2z*d3y) * invDet;
    invJ[3] = (d3y*d1z - d3z*d1y) * invDet; // Transposed for col-major? 
    // Wait, standard inverse formula Adj(J)_ji.
    // Let's stick to the implementation verified in assemble_winding
    // Adj(0,0) = +(d2y*d3z - d2z*d3y) -> invJ[0]
    // Adj(0,1) = -(d2x*d3z - d2z*d3x) -> invJ[3] (Row 0, Col 1 in C is index 3 in Col-Major)
    invJ[6] = (d2x*d3y - d2y*d3x) * invDet; // Row 0, Col 2 -> index 6

    // Row 1
    invJ[1] = (d3z*d1y - d1z*d3y) * invDet; // Adj(1,0) -> Row 1, Col 0 -> Index 1
    invJ[4] = (d1x*d3z - d1z*d3x) * invDet; 
    invJ[7] = (d3x*d1y - d1x*d3y) * invDet;

    // Row 2
    invJ[2] = (d1z*d2y - d1y*d2z) * invDet;
    invJ[5] = (d1y*d2x - d1x*d2y) * invDet;
    invJ[8] = (d1x*d2z - d1z*d2x) * invDet;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs < 7) mexErrMsgTxt("Input arguments: P, T, CellDofs, Signs, Js, Qw, ValRef");

    // Inputs
    double* P = mxGetPr(prhs[0]);
    double* T = mxGetPr(prhs[1]);
    double* CellDofs = mxGetPr(prhs[2]);
    double* Signs = mxGetPr(prhs[3]);
    double* Js = mxGetPr(prhs[4]);
    double* Qw = mxGetPr(prhs[5]);
    double* ValRef = mxGetPr(prhs[6]);

    size_t numElems = mxGetN(prhs[1]);
    size_t n_q = mxGetM(prhs[5]);

    // Output Vectors (Triplets)
    // Each element contributes 6 entries
    plhs[0] = mxCreateDoubleMatrix(6 * numElems, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(6 * numElems, 1, mxREAL);

    double* I = mxGetPr(plhs[0]);
    double* V = mxGetPr(plhs[1]);

    #pragma omp parallel for
    for (size_t e = 0; e < numElems; ++e) {
        // 1. Get Geometry
        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        
        double detJ;
        double invJ[9]; 
        compute_jacobian_and_invJ(P, T_local, detJ, invJ);
        double absDetJ = std::abs(detJ);

        // 2. Get Source J for this element
        double Jx = Js[0 + e*3];
        double Jy = Js[1 + e*3];
        double Jz = Js[2 + e*3];
        
        // Skip empty source (Optimization)
        if (std::abs(Jx) < 1e-15 && std::abs(Jy) < 1e-15 && std::abs(Jz) < 1e-15) {
             // Fill zeros
             double* dofs = &CellDofs[e*6];
             for(int k=0; k<6; k++) {
                 size_t idx = e*6 + k;
                 I[idx] = dofs[k]; 
                 V[idx] = 0.0;
             }
             continue;
        }

        // 3. Integration
        double Fe[6] = {0};

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ;
            
            // Nedelec Pushforward: N_phy = N_ref * J^{-1}
            // N_phy is [6 x 3]
            // We calculate N_phy . J_vec directly
            
            size_t q_offset = q * 18; // 6*3
            
            for (int i = 0; i < 6; ++i) {
                // Get N_ref[i, :]
                double nr_x = ValRef[i + 0*6 + q_offset];
                double nr_y = ValRef[i + 1*6 + q_offset];
                double nr_z = ValRef[i + 2*6 + q_offset];
                
                // N_phy_x = nr_x * invJ[0] + nr_y * invJ[3] + nr_z * invJ[6]
                double nx = nr_x * invJ[0] + nr_y * invJ[3] + nr_z * invJ[6];
                double ny = nr_x * invJ[1] + nr_y * invJ[4] + nr_z * invJ[7];
                double nz = nr_x * invJ[2] + nr_y * invJ[5] + nr_z * invJ[8];
                
                // Dot product with Js
                double dot_val = nx * Jx + ny * Jy + nz * Jz;
                
                Fe[i] += w * dot_val;
            }
        }

        // 4. Signs and Output
        double* s = &Signs[e*6];
        double* dofs = &CellDofs[e*6];
        
        for(int k=0; k<6; k++) {
            size_t idx = e*6 + k;
            I[idx] = dofs[k];
            V[idx] = Fe[k] * s[k];
        }
    }
}