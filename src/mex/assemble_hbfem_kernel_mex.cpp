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
    // 0-4: P, T, Dofs, Signs, Tags
    // 5: SolH_Real, 6: SolH_Imag (Split Complex)
    // 7: Dmat_Real, 8: Dmat_Imag (Split Complex)
    // 9: Pmat_Real, 10: Pmat_Imag (Split Complex)
    // 11: Qw, 12: CurlRef
    // 13-20: Material Data (8 arrays)
    // 21: RequestJacobian
    // 22: Scalings
    
    if (nrhs < 23) mexErrMsgIdAndTxt("MagPackage:HBFEM", "Need 23 inputs (Split Complex).");
    
    double *P = mxGetPr(prhs[0]);
    double *T = mxGetPr(prhs[1]); size_t n_elems = mxGetN(prhs[1]);
    double *Dofs = mxGetPr(prhs[2]);
    double *Signs = mxGetPr(prhs[3]);
    double *Tags = mxGetPr(prhs[4]);
    
    // Complex Solution
    double *SolH_R = mxGetPr(prhs[5]);
    double *SolH_I = mxGetPr(prhs[6]);
    size_t numDofs = mxGetM(prhs[5]);
    size_t numHarm = mxGetN(prhs[5]);
    bool has_sol_I = !mxIsEmpty(prhs[6]);
    
    // Complex Dmat
    double *Dmat_R = mxGetPr(prhs[7]);
    double *Dmat_I = mxGetPr(prhs[8]);
    size_t numTime = mxGetM(prhs[7]);
    bool has_dmat_I = !mxIsEmpty(prhs[8]);
    
    // Complex Pmat
    double *Pmat_R = mxGetPr(prhs[9]);
    double *Pmat_I = mxGetPr(prhs[10]);
    bool has_pmat_I = !mxIsEmpty(prhs[10]);
    
    double *Qw = mxGetPr(prhs[11]); size_t n_q = mxGetM(prhs[11]) * mxGetN(prhs[11]);
    double *CurlRef = mxGetPr(prhs[12]);

    MexMaterialData matData;
    matData.linearNu    = mxGetPr(prhs[13]);
    matData.isNonlinear = mxGetPr(prhs[14]);
    matData.maxBSq      = mxGetPr(prhs[15]);
    matData.breaks      = mxGetPr(prhs[16]);
    matData.coefs       = mxGetPr(prhs[17]);
    matData.splineStart = mxGetPr(prhs[18]);
    matData.splineCount = mxGetPr(prhs[19]);
    // input 20 is placeholder numTags
    
    bool requestJ = (bool)mxGetScalar(prhs[21]);
    double *Scalings = mxGetPr(prhs[22]);

    // Outputs
    // K: [I, J, V] (V is Real for DC approximation)
    size_t n_K = requestJ ? (n_elems * 36 * numHarm) : 0;
    size_t n_R = n_elems * 6 * numHarm;
    
    plhs[0] = mxCreateDoubleMatrix(n_K, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_K, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_K, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(n_R, 1, mxREAL); // R_rows
    plhs[4] = mxCreateDoubleMatrix(n_R, 1, mxREAL); // R_cols
    plhs[5] = mxCreateDoubleMatrix(n_R, 1, mxREAL); // R_val_Real
    plhs[6] = mxCreateDoubleMatrix(n_R, 1, mxREAL); // R_val_Imag (New Output)

    double *I_out = mxGetPr(plhs[0]);
    double *J_out = mxGetPr(plhs[1]);
    double *V_out = mxGetPr(plhs[2]);
    double *R_rows = mxGetPr(plhs[3]);
    double *R_cols = mxGetPr(plhs[4]);
    double *R_val_R = mxGetPr(plhs[5]);
    double *R_val_I = mxGetPr(plhs[6]);

    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4];
        
        // 1. Fetch Solution A_local [6 x NumHarm] (Complex)
        // Store Real and Imag separate
        std::vector<double> A_loc_R(6 * numHarm);
        std::vector<double> A_loc_I(6 * numHarm);
        
        double dofs[6], s[6];
        for(int k=0; k<6; k++) {
            dofs[k] = Dofs[k + e*6];
            s[k] = Signs[k + e*6];
            int global_id = (int)dofs[k] - 1;
            
            for (int h=0; h<numHarm; ++h) {
                int idx = global_id + h*numDofs;
                A_loc_R[k + h*6] = SolH_R[idx] * s[k];
                if (has_sol_I) A_loc_I[k + h*6] = SolH_I[idx] * s[k];
                else           A_loc_I[k + h*6] = 0.0;
            }
        }
        
        // 2. Compute A_time [6 x NumTime] (Real Physical Field)
        // A(t) = Real( A_phasor * Dmat^T )
        // Dmat^T = (Dr + j Di)^T = Dr^T + j Di^T
        // A * D^T = (Ar + j Ai)(Dr^T + j Di^T)
        //         = (Ar Dr^T - Ai Di^T) + j(...)
        // We only need Real part: Ar * Dr^T - Ai * Di^T
        
        std::vector<double> A_time(6 * numTime);
        for (int t=0; t<numTime; ++t) {
            for (int k=0; k<6; ++k) {
                double sum = 0;
                for (int h=0; h<numHarm; ++h) {
                    // Dmat(Row=t, Col=h)
                    double dr = Dmat_R[t + h*numTime];
                    double ar = A_loc_R[k + h*6];
                    
                    double term = ar * dr;
                    
                    if (has_sol_I && has_dmat_I) {
                        double di = Dmat_I[t + h*numTime];
                        double ai = A_loc_I[k + h*6];
                        term -= ai * di; 
                    }
                    sum += term;
                }
                A_time[k + t*6] = sum;
            }
        }
        
        int tag = (int)Tags[e];
        double J_mat[9], detJ, invDetJ;
        compute_jacobian_utils(P, T_local, J_mat, detJ, invDetJ);
        double absDetJ = std::abs(detJ);
        
        std::vector<double> Re_harm_R(6 * numHarm, 0.0);
        std::vector<double> Re_harm_I(6 * numHarm, 0.0);
        std::vector<double> Ke_base(36, 0.0);

        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ;
            
            double C_phy[18]; 
             for (int i = 0; i < 6; ++i) {
                double cr_x = CurlRef[i + 0*6 + q*18];
                double cr_y = CurlRef[i + 1*6 + q*18];
                double cr_z = CurlRef[i + 2*6 + q*18];
                C_phy[i+0*6] = (cr_x * J_mat[0] + cr_y * J_mat[1] + cr_z * J_mat[2]) * invDetJ;
                C_phy[i+1*6] = (cr_x * J_mat[3] + cr_y * J_mat[4] + cr_z * J_mat[5]) * invDetJ;
                C_phy[i+2*6] = (cr_x * J_mat[6] + cr_y * J_mat[7] + cr_z * J_mat[8]) * invDetJ;
            }
            
            double sum_nu = 0;
            
            for (int t=0; t<numTime; ++t) {
                double Bx=0, By=0, Bz=0;
                for(int i=0; i<6; ++i) {
                    double val = A_time[i + t*6];
                    Bx += C_phy[i+0*6] * val;
                    By += C_phy[i+1*6] * val;
                    Bz += C_phy[i+2*6] * val;
                }
                double B_sq = Bx*Bx + By*By + Bz*Bz;
                
                double nu, dnu;
                evaluate_nu_derivative(B_sq, tag, matData, nu, dnu);
                sum_nu += nu;
                
                double Hx = nu * Bx;
                double Hy = nu * By;
                double Hz = nu * Bz;
                
                // Add to Residual (Frequency Domain)
                // H_freq = H_time * Pmat^T
                // Pmat^T = Pr^T + j Pi^T
                // H_freq = Ht * Pr^T + j * Ht * Pi^T
                // Re_R += w * (Ci.H) * Pr
                // Re_I += w * (Ci.H) * Pi
                
                for (int h=0; h<numHarm; ++h) {
                    double pr = Pmat_R[h + t*numHarm]; 
                    double dot = 0;
                    // Optimization: Precalc dot product
                    for (int i=0; i<6; ++i) {
                         dot = C_phy[i+0*6]*Hx + C_phy[i+1*6]*Hy + C_phy[i+2*6]*Hz;
                         
                         // Accumulate Real
                         Re_harm_R[i + h*6] += w * dot * pr;
                         
                         // Accumulate Imag
                         if (has_pmat_I) {
                             double pi = Pmat_I[h + t*numHarm];
                             Re_harm_I[i + h*6] += w * dot * pi;
                         }
                    }
                }
            }
            
            if (requestJ) {
                double nu_dc = sum_nu / numTime;
                for (int c=0; c<6; ++c) {
                    for (int r=0; r<6; ++r) {
                        double dot = C_phy[r+0*6]*C_phy[c+0*6] + 
                                     C_phy[r+1*6]*C_phy[c+1*6] + 
                                     C_phy[r+2*6]*C_phy[c+2*6];
                        Ke_base[r + c*6] += w * nu_dc * dot;
                    }
                }
            }
        }
        
        size_t off_R = e * 6 * numHarm;
        int idxR = 0;
        for (int h=0; h<numHarm; ++h) {
            for (int k=0; k<6; ++k) {
                R_rows[off_R + idxR] = dofs[k];
                R_cols[off_R + idxR] = (double)(h + 1); 
                R_val_R[off_R + idxR] = Re_harm_R[k + h*6] * s[k];
                R_val_I[off_R + idxR] = Re_harm_I[k + h*6] * s[k];
                idxR++;
            }
        }
        
        if (requestJ) {
            size_t off_K = e * 36 * numHarm;
            int idxK = 0;
            for (int h=0; h<numHarm; ++h) {
                double factor = Scalings[h];
                for (int c=0; c<6; ++c) {
                    for (int r=0; r<6; ++r) {
                        double val = Ke_base[r + c*6] * factor;
                        I_out[off_K + idxK] = dofs[r] + (h * numDofs); 
                        J_out[off_K + idxK] = dofs[c] + (h * numDofs);
                        V_out[off_K + idxK] = val * s[r] * s[c];
                        idxK++;
                    }
                }
            }
        }
    }
}