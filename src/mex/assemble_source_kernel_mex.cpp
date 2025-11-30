#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>

// 复用 Mass Kernel 中的 Jacobian 计算函数
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

    detJ = d1x*(d2y*d3z - d2z*d3y) - d1y*(d2x*d3z - d2z*d3x) + d1z*(d2x*d3y - d2y*d3x);
    double invDetJ = 1.0 / detJ;

    iJt_trans[0] = (d2y*d3z - d2z*d3y) * invDetJ;
    iJt_trans[1] = (d2z*d3x - d2x*d3z) * invDetJ;
    iJt_trans[2] = (d2x*d3y - d2y*d3x) * invDetJ;
    iJt_trans[3] = (d3y*d1z - d3z*d1y) * invDetJ;
    iJt_trans[4] = (d3z*d1x - d3x*d1z) * invDetJ;
    iJt_trans[5] = (d3x*d1y - d3y*d1x) * invDetJ;
    iJt_trans[6] = (d1y*d2z - d1z*d2y) * invDetJ;
    iJt_trans[7] = (d1z*d2x - d1x*d2z) * invDetJ;
    iJt_trans[8] = (d1x*d2y - d1y*d2x) * invDetJ;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Inputs: P, T, Dofs, Signs, Tags, SourceMap, Qw, ValRef
    if (nrhs < 8) mexErrMsgIdAndTxt("MagPackage:Source", "Need 8 inputs.");

    double *P = mxGetPr(prhs[0]);
    double *T = mxGetPr(prhs[1]); 
    size_t n_elems = mxGetN(prhs[1]);
    
    double *Dofs = mxGetPr(prhs[2]);
    double *Signs = mxGetPr(prhs[3]);
    double *Tags = mxGetPr(prhs[4]);
    
    double *SrcMap = mxGetPr(prhs[5]); // [MaxTag x 3]
    size_t n_src_rows = mxGetM(prhs[5]);
    
    double *Qw = mxGetPr(prhs[6]);
    size_t n_q = mxGetM(prhs[6]) * mxGetN(prhs[6]);
    
    double *ValRef = mxGetPr(prhs[7]);

    // Output: I, V (F vector indices and values)
    // Alloc max size
    size_t n_entries = n_elems * 6;
    plhs[0] = mxCreateDoubleMatrix(n_entries, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(n_entries, 1, mxREAL); // V
    
    double *I_out = mxGetPr(plhs[0]);
    double *V_out = mxGetPr(plhs[1]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        int tag = (int)Tags[e]; // Tag is 1-based index for SrcMap
        
        // Check tag validity and if source is zero
        bool is_zero = true;
        double Jx=0, Jy=0, Jz=0;
        
        if (tag >= 1 && tag <= n_src_rows) {
            Jx = SrcMap[tag-1 + 0*n_src_rows];
            Jy = SrcMap[tag-1 + 1*n_src_rows];
            Jz = SrcMap[tag-1 + 2*n_src_rows];
            if (std::abs(Jx) > 1e-14 || std::abs(Jy) > 1e-14 || std::abs(Jz) > 1e-14) {
                is_zero = false;
            }
        }
        
        if (is_zero) {
             // Fill zero
             double dofs[6];
             for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
             size_t offset = e * 6;
             for(int k=0; k<6; k++) {
                 I_out[offset + k] = dofs[k];
                 V_out[offset + k] = 0.0;
             }
             continue;
        }

        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        
        double detJ;
        double iJt_trans[9]; 
        compute_jacobian_and_iJt_trans(P, T_local, detJ, iJt_trans);
        double absDetJ = std::abs(detJ);

        double Fe[6] = {0};

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ;
            
            // N_phy = N_ref * iJt_trans
            for (int i = 0; i < 6; ++i) {
                double nr_x = ValRef[i + 0*6 + q*18];
                double nr_y = ValRef[i + 1*6 + q*18];
                double nr_z = ValRef[i + 2*6 + q*18];
                
                double nx = nr_x * iJt_trans[0] + nr_y * iJt_trans[3] + nr_z * iJt_trans[6];
                double ny = nr_x * iJt_trans[1] + nr_y * iJt_trans[4] + nr_z * iJt_trans[7];
                double nz = nr_x * iJt_trans[2] + nr_y * iJt_trans[5] + nr_z * iJt_trans[8];
                
                // Fe += w * dot(N, J)
                Fe[i] += w * (nx * Jx + ny * Jy + nz * Jz);
            }
        }
        
        double s[6];
        for(int k=0; k<6; k++) s[k] = Signs[k + e*6];
        double dofs[6];
        for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
        
        size_t offset = e * 6;
        for(int k=0; k<6; k++) {
            I_out[offset + k] = dofs[k];
            V_out[offset + k] = Fe[k] * s[k];
        }
    }
}