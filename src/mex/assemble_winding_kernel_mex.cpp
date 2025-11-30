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

    // invJ 存储 J^{-1} (Row-Major)
    // Row 1 = cross(d2, d3) / detJ
    invJ[0] = (d2y*d3z - d2z*d3y) * invDetJ;
    invJ[1] = (d2z*d3x - d2x*d3z) * invDetJ;
    invJ[2] = (d2x*d3y - d2y*d3x) * invDetJ;
    
    // Row 2 = cross(d3, d1) / detJ
    invJ[3] = (d3y*d1z - d3z*d1y) * invDetJ;
    invJ[4] = (d3z*d1x - d3x*d1z) * invDetJ;
    invJ[5] = (d3x*d1y - d3y*d1x) * invDetJ;
    
    // Row 3 = cross(d1, d2) / detJ
    invJ[6] = (d1y*d2z - d1z*d2y) * invDetJ;
    invJ[7] = (d1z*d2x - d1x*d2z) * invDetJ;
    invJ[8] = (d1x*d2y - d1y*d2x) * invDetJ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Inputs: 
    // 0: P, 1: T, 2: Dofs, 3: Signs, 4: Tags, 
    // 5: TargetRegion, 6: DirectionField (3xNe or empty), 7: ConstDir (3x1), 
    // 8: CurrentScale, 9: Qw, 10: ValRef
    if (nrhs < 11) mexErrMsgIdAndTxt("MagPackage:Winding", "Need 11 inputs.");

    double *P = mxGetPr(prhs[0]);
    double *T = mxGetPr(prhs[1]); 
    size_t n_elems = mxGetN(prhs[1]);
    
    double *Dofs = mxGetPr(prhs[2]);
    double *Signs = mxGetPr(prhs[3]);
    double *Tags = mxGetPr(prhs[4]);
    
    int target_region = (int)mxGetScalar(prhs[5]);
    
    double *DirField = nullptr;
    bool has_field = false;
    if (!mxIsEmpty(prhs[6])) {
        DirField = mxGetPr(prhs[6]);
        has_field = true;
    }
    
    double *ConstDir = mxGetPr(prhs[7]);
    double current_scale = mxGetScalar(prhs[8]);
    
    double *Qw = mxGetPr(prhs[9]);
    size_t n_q = mxGetM(prhs[9]) * mxGetN(prhs[9]);
    
    double *ValRef = mxGetPr(prhs[10]);

    // Output: I, V
    size_t n_entries = n_elems * 6;
    plhs[0] = mxCreateDoubleMatrix(n_entries, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_entries, 1, mxREAL);
    
    double *I_out = mxGetPr(plhs[0]);
    double *V_out = mxGetPr(plhs[1]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        // Filter by Tag
        if ((int)Tags[e] != target_region) {
            // Fill zeros
            size_t offset = e * 6;
            double dofs[6];
            for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
            for(int k=0; k<6; k++) {
                I_out[offset + k] = dofs[k];
                V_out[offset + k] = 0.0;
            }
            continue;
        }
        
        // Determine Direction Vector
        double tx, ty, tz;
        if (has_field) {
            tx = DirField[0 + e*3];
            ty = DirField[1 + e*3];
            tz = DirField[2 + e*3];
        } else {
            tx = ConstDir[0]; ty = ConstDir[1]; tz = ConstDir[2];
        }
        
        double Jx = tx * current_scale;
        double Jy = ty * current_scale;
        double Jz = tz * current_scale;

        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        
        double detJ;
        double invJ[9]; 
        compute_jacobian_and_invJ(P, T_local, detJ, invJ);
        double absDetJ = std::abs(detJ);

        double Ce[6] = {0};

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ;
            
            // N_phy = N_ref * J^{-1}
            for (int i = 0; i < 6; ++i) {
                double nr_x = ValRef[i + 0*6 + q*18];
                double nr_y = ValRef[i + 1*6 + q*18];
                double nr_z = ValRef[i + 2*6 + q*18];
                
                double nx = nr_x * invJ[0] + nr_y * invJ[3] + nr_z * invJ[6];
                double ny = nr_x * invJ[1] + nr_y * invJ[4] + nr_z * invJ[7];
                double nz = nr_x * invJ[2] + nr_y * invJ[5] + nr_z * invJ[8];
                
                // Ce += w * dot(N, J)
                Ce[i] += w * (nx * Jx + ny * Jy + nz * Jz);
            }
        }
        
        double s[6];
        for(int k=0; k<6; k++) s[k] = Signs[k + e*6];
        double dofs[6];
        for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
        
        size_t offset = e * 6;
        for(int k=0; k<6; k++) {
            I_out[offset + k] = dofs[k];
            V_out[offset + k] = Ce[k] * s[k];
        }
    }
}