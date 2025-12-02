/**
 * assemble_hbfem_kernel_mex.cpp (Final Robust Version)
 * 功能: HBFEM 刚度矩阵与残差组装
 * * 修复日志:
 * 1. [Fix] 修复全局自由度推断逻辑，防止生成尺寸错误的稀疏矩阵导致 MATLAB "维度不匹配" 错误。
 * 2. [Safety] 强制所有索引计算使用 1-based double 类型，严格对齐 MATLAB sparse() 接口。
 * 3. [Robust] 增加输入维度校验，防止 SolR 与 Harmonics 配置不一致。
 */

#include "mex.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>
#include <omp.h>

// 包含用户提供的头文件
#include "MexElemUtils.hpp"
#include "MexMaterialUtils.hpp"

using namespace std;

// 复数运算辅助宏
#define C_MUL_R(ar, ai, br, bi) ((ar)*(br) - (ai)*(bi))
#define C_MUL_I(ar, ai, br, bi) ((ar)*(bi) + (ai)*(br))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // --- 1. 参数完整性检查 ---
    if (nrhs < 24) mexErrMsgIdAndTxt("HBFEM:Args", "Need 24 inputs.");

    // 基础数据指针
    const double* P = mxGetPr(prhs[0]);
    const double* T = mxGetPr(prhs[1]);
    const double* Dofs = mxGetPr(prhs[2]); // [6 x Nelems] (1-based indices from MATLAB)
    const double* Signs = mxGetPr(prhs[3]);
    const double* Tags = mxGetPr(prhs[4]);
    
    // DFT/IDFT 矩阵
    const double* Dmat_R = mxGetPr(prhs[7]); // [Nt x Nh]
    const double* Dmat_I = mxGetPr(prhs[8]);
    size_t numTime = mxGetM(prhs[7]);
    size_t numHarm = mxGetN(prhs[7]); // K (Harmonics Count)
    
    const double* Pmat_R = mxGetPr(prhs[9]); 
    const double* Pmat_I = mxGetPr(prhs[10]);
    
    // HBFEM 解向量
    const double* SolR = mxGetPr(prhs[5]);
    const double* SolI = mxGetPr(prhs[6]);
    size_t solRows = mxGetM(prhs[5]);
    size_t solCols = mxGetN(prhs[5]);

    // [CRITICAL logic] 推断空间自由度数 (N)
    // 逻辑: SolR 可以是 [N x K] 矩阵 或者 [NK x 1] 向量
    size_t numGlobalDofs = 0;
    if (solCols == numHarm) {
        numGlobalDofs = solRows; // Standard Matrix Form
    } else if (solCols == 1) {
        if (solRows % numHarm != 0) {
             mexErrMsgIdAndTxt("HBFEM:DimError", "SolR vector length must be multiple of NumHarmonics.");
        }
        numGlobalDofs = solRows / numHarm; // Flattened Vector Form
    } else {
        mexErrMsgIdAndTxt("HBFEM:DimMismatch", 
            "Solution dimension mismatch. Expected NxK or (NK)x1. "
            "Sol: %dx%d, Harmonics: %d", solRows, solCols, numHarm);
    }
    
    // 积分与基函数
    const double* Qw = mxGetPr(prhs[11]);
    size_t numQP = mxGetM(prhs[11]);
    const double* RefCurl = mxGetPr(prhs[12]); 
    
    // 材料评估器初始化
    MaterialEvaluator matEval(
        mxGetPr(prhs[13]), mxGetPr(prhs[14]), mxGetPr(prhs[15]), 
        mxGetPr(prhs[16]), mxGetPr(prhs[17]), mxGetPr(prhs[18]), 
        mxGetPr(prhs[19]), mxGetPr(prhs[20]), (int)mxGetM(prhs[13])-1
    );

    // 求解控制
    bool calcJ = (bool)mxGetScalar(prhs[22]);
    const double* Scalings = mxGetPr(prhs[23]); 

    size_t numElems = mxGetN(prhs[1]);
    
    // --- 2. 输出内存分配 ---
    size_t blockSize = 6 * numHarm;
    size_t estNNZ = numElems * blockSize * blockSize;
    size_t estRes = numElems * blockSize;

    // Jacobian Triplets
    plhs[0] = mxCreateDoubleMatrix(calcJ ? estNNZ : 0, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(calcJ ? estNNZ : 0, 1, mxREAL); // J
    plhs[2] = mxCreateDoubleMatrix(calcJ ? estNNZ : 0, 1, mxCOMPLEX); // V
    
    // Residual (COO)
    plhs[3] = mxCreateDoubleMatrix(estRes, 1, mxREAL); // R_rows
    plhs[4] = mxCreateDoubleMatrix(estRes, 1, mxREAL); // R_cols
    plhs[5] = mxCreateDoubleMatrix(estRes, 1, mxREAL); // R_val_R
    plhs[6] = mxCreateDoubleMatrix(estRes, 1, mxREAL); // R_val_I

    double* outI = mxGetPr(plhs[0]);
    double* outJ = mxGetPr(plhs[1]);
    double* outV_r = mxGetPr(plhs[2]);
    double* outV_i = mxGetPi(plhs[2]); // -R2017b compatible

    double* outR_r = mxGetPr(plhs[3]);
    double* outR_c = mxGetPr(plhs[4]);
    double* outR_vR = mxGetPr(plhs[5]);
    double* outR_vI = mxGetPr(plhs[6]);

    // 全局计数器 (OpenMP atomic update)
    size_t global_nnz = 0;
    size_t global_res = 0;

    // --- 3. 并行组装 ---
    #pragma omp parallel
    {
        // 线程局部缓冲区 (复用内存)
        vector<double> B_freq_r(3 * numHarm), B_freq_i(3 * numHarm);
        vector<double> B_time_x(numTime), B_time_y(numTime), B_time_z(numTime);
        vector<double> nu_t(numTime), dnu_t(numTime); 
        vector<double> M_iso_r(numHarm * numHarm), M_iso_i(numHarm * numHarm);
        
        vector<double> Ke_r, Ke_i;
        if (calcJ) { Ke_r.resize(blockSize * blockSize); Ke_i.resize(blockSize * blockSize); }
        vector<double> Re_r(blockSize), Re_i(blockSize);
        vector<double> A_r(6 * numHarm), A_i(6 * numHarm);

        #pragma omp for
        for (long e = 0; e < numElems; e++) {
            
            // A. Geometry
            double p_local[12];
            for (int n = 0; n < 4; n++) {
                long nid = (long)T[n + e*4] - 1; 
                p_local[n*3+0] = P[nid*3+0];
                p_local[n*3+1] = P[nid*3+1];
                p_local[n*3+2] = P[nid*3+2];
            }
            
            Mat3x3 J_mat, J_inv;
            double detJ = MexElemUtils::compute_jacobian_3d(p_local, J_mat);
            MexElemUtils::compute_inverse_3x3(J_mat, detJ, J_inv);
            double absDetJ = std::abs(detJ);
            double invDetJ = 1.0 / detJ;
            
            int tag = (int)Tags[e];
            double elem_signs[6];
            long elem_dofs[6]; 
            
            for(int i=0; i<6; i++) {
                elem_signs[i] = Signs[i + e*6];
                elem_dofs[i]  = (long)Dofs[i + e*6]; // 1-based Global DoF
            }

            // Extract Solution
            for(int i=0; i<6; i++) {
                long gdof = elem_dofs[i] - 1; // 0-based for C++ indexing
                
                // Safety: check index bounds
                if (gdof >= numGlobalDofs) {
                    // This implies DofHandler map is larger than Solution vector
                    // Just define as 0 to avoid crash
                     fill(A_r.begin()+i*numHarm, A_r.begin()+(i+1)*numHarm, 0.0);
                     fill(A_i.begin()+i*numHarm, A_i.begin()+(i+1)*numHarm, 0.0);
                     continue;
                }

                for(int h=0; h<numHarm; h++) {
                    size_t idx;
                    if (solCols == numHarm) {
                        idx = gdof + h * solRows; // Matrix: col-major
                    } else {
                        idx = gdof + h * numGlobalDofs; // Vector: stacked
                    }
                    A_r[i + h*6] = SolR[idx];
                    A_i[i + h*6] = SolI[idx];
                }
            }

            // Clear Buffers
            if (calcJ) { fill(Ke_r.begin(), Ke_r.end(), 0.0); fill(Ke_i.begin(), Ke_i.end(), 0.0); }
            fill(Re_r.begin(), Re_r.end(), 0.0); fill(Re_i.begin(), Re_i.end(), 0.0);

            // B. Quadrature Loop
            for(int q=0; q<numQP; q++) {
                double w = Qw[q] * absDetJ;
                
                // Piola
                double C_phy[6][3];
                size_t q_offset = q * 18; 
                for(int i=0; i<6; i++) {
                    double ref_c[3] = {RefCurl[i+q_offset], RefCurl[i+6+q_offset], RefCurl[i+12+q_offset]};
                    double temp[3];
                    J_mat.multVec(ref_c, temp); 
                    C_phy[i][0] = temp[0] * invDetJ;
                    C_phy[i][1] = temp[1] * invDetJ;
                    C_phy[i][2] = temp[2] * invDetJ;
                }
                
                // B_freq
                fill(B_freq_r.begin(), B_freq_r.end(), 0.0);
                fill(B_freq_i.begin(), B_freq_i.end(), 0.0);
                for(int h=0; h<numHarm; h++) {
                    for(int i=0; i<6; i++) {
                        double ar = A_r[i + h*6], ai = A_i[i + h*6];
                        for(int d=0; d<3; d++) {
                            B_freq_r[d + h*3] += C_phy[i][d] * ar;
                            B_freq_i[d + h*3] += C_phy[i][d] * ai;
                        }
                    }
                }
                
                // B_time (IDFT)
                for(int k=0; k<numTime; k++) {
                    double bx=0, by=0, bz=0;
                    for(int h=0; h<numHarm; h++) {
                        double Dr = Dmat_R[k + h*numTime], Di = Dmat_I[k + h*numTime];
                        double Brx = B_freq_r[0+h*3], Bix = B_freq_i[0+h*3];
                        bx += C_MUL_R(Dr, Di, Brx, Bix);
                        double Bry = B_freq_r[1+h*3], Biy = B_freq_i[1+h*3];
                        by += C_MUL_R(Dr, Di, Bry, Biy);
                        double Brz = B_freq_r[2+h*3], Biz = B_freq_i[2+h*3];
                        bz += C_MUL_R(Dr, Di, Brz, Biz);
                    }
                    B_time_x[k]=bx; B_time_y[k]=by; B_time_z[k]=bz;
                }
                
                // Material
                for(int k=0; k<numTime; k++) {
                    double b2 = B_time_x[k]*B_time_x[k] + B_time_y[k]*B_time_y[k] + B_time_z[k]*B_time_z[k];
                    double dnu_val = 0;
                    nu_t[k] = matEval.evaluate(tag, b2, dnu_val);
                    dnu_t[k] = 2.0 * dnu_val; 
                }
                
                // Residual (DFT)
                for(int h=0; h<numHarm; h++) {
                    double H_fr_x=0, H_fi_x=0, H_fr_y=0, H_fi_y=0, H_fr_z=0, H_fi_z=0;
                    for(int k=0; k<numTime; k++) {
                        double Pr=Pmat_R[h+k*numHarm], Pi=Pmat_I[h+k*numHarm];
                        double nuk=nu_t[k];
                        double Htx=nuk*B_time_x[k], Hty=nuk*B_time_y[k], Htz=nuk*B_time_z[k];
                        H_fr_x += Pr*Htx; H_fi_x += Pi*Htx;
                        H_fr_y += Pr*Hty; H_fi_y += Pi*Hty;
                        H_fr_z += Pr*Htz; H_fi_z += Pi*Htz;
                    }
                    double scale = Scalings[h];
                    for(int i=0; i<6; i++) {
                        double dot_r = C_phy[i][0]*H_fr_x + C_phy[i][1]*H_fr_y + C_phy[i][2]*H_fr_z;
                        double dot_i = C_phy[i][0]*H_fi_x + C_phy[i][1]*H_fi_y + C_phy[i][2]*H_fi_z;
                        double s_fac = elem_signs[i] * w * scale;
                        Re_r[i + h*6] += dot_r * s_fac;
                        Re_i[i + h*6] += dot_i * s_fac;
                    }
                }
                
                // Jacobian
                if (calcJ) {
                    fill(M_iso_r.begin(), M_iso_r.end(), 0.0); fill(M_iso_i.begin(), M_iso_i.end(), 0.0);
                    for(int mh=0; mh<numHarm; mh++) {
                        for(int nh=0; nh<numHarm; nh++) {
                            double val_r=0, val_i=0;
                            for(int k=0; k<numTime; k++) {
                                double Pr = Pmat_R[mh + k*numHarm], Pi = Pmat_I[mh + k*numHarm];
                                double Dr = Dmat_R[k + nh*numTime], Di = Dmat_I[k + nh*numTime];
                                double nuk = nu_t[k];
                                double tmp_r = C_MUL_R(Pr, Pi, Dr, Di);
                                double tmp_i = C_MUL_I(Pr, Pi, Dr, Di);
                                val_r += tmp_r * nuk; val_i += tmp_i * nuk;
                            }
                            M_iso_r[mh+nh*numHarm]=val_r; M_iso_i[mh+nh*numHarm]=val_i;
                        }
                    }
                    
                    for(int i=0; i<6; i++) {
                        for(int j=0; j<6; j++) {
                            double geom_dot = C_phy[i][0]*C_phy[j][0] + C_phy[i][1]*C_phy[j][1] + C_phy[i][2]*C_phy[j][2];
                            double s_fac = elem_signs[i] * elem_signs[j] * w;
                            
                            // Iso
                            for(int mh=0; mh<numHarm; mh++) {
                                double scale = Scalings[mh] * s_fac;
                                for(int nh=0; nh<numHarm; nh++) {
                                    int idx = (i + mh*6) + (j + nh*6)*blockSize;
                                    double mr = M_iso_r[mh + nh*numHarm], mi = M_iso_i[mh + nh*numHarm];
                                    Ke_r[idx] += mr * geom_dot * scale;
                                    Ke_i[idx] += mi * geom_dot * scale;
                                }
                            }
                            // Aniso
                            for(int mh=0; mh<numHarm; mh++) {
                                for(int nh=0; nh<numHarm; nh++) {
                                    double sum_r=0, sum_i=0;
                                    for(int k=0; k<numTime; k++) {
                                        double CiBk = C_phy[i][0]*B_time_x[k] + C_phy[i][1]*B_time_y[k] + C_phy[i][2]*B_time_z[k];
                                        double CjBk = C_phy[j][0]*B_time_x[k] + C_phy[j][1]*B_time_y[k] + C_phy[j][2]*B_time_z[k];
                                        double factor = dnu_t[k] * CiBk * CjBk;
                                        double Pr = Pmat_R[mh + k*numHarm], Pi = Pmat_I[mh + k*numHarm];
                                        double Dr = Dmat_R[k + nh*numTime], Di = Dmat_I[k + nh*numTime];
                                        double tmp_r = C_MUL_R(Pr, Pi, Dr, Di);
                                        double tmp_i = C_MUL_I(Pr, Pi, Dr, Di);
                                        sum_r += tmp_r * factor; sum_i += tmp_i * factor;
                                    }
                                    int idx = (i + mh*6) + (j + nh*6)*blockSize;
                                    double scale = Scalings[mh] * s_fac;
                                    Ke_r[idx] += sum_r * scale; Ke_i[idx] += sum_i * scale;
                                }
                            }
                        }
                    }
                }
            } // End Quad

            // C. Write Output (Safe Types)
            size_t r_offset;
            #pragma omp atomic capture
            { r_offset = global_res; global_res += blockSize; }
            for(int k=0; k<blockSize; k++) {
                int i=k%6, h=k/6;
                // R_rows: Global DoF (1-based)
                outR_r[r_offset+k] = (double)elem_dofs[i]; 
                // R_cols: Harmonic (1-based)
                outR_c[r_offset+k] = (double)(h+1); 
                outR_vR[r_offset+k] = Re_r[k]; outR_vI[r_offset+k] = Re_i[k];
            }

            if (calcJ) {
                size_t j_offset;
                size_t n_local_triplets = blockSize * blockSize;
                #pragma omp atomic capture
                { j_offset = global_nnz; global_nnz += n_local_triplets; }
                
                size_t count=0;
                for(int c=0; c<blockSize; c++) {
                    int j_b = c%6, h_c = c/6;
                    double gdof_j = (double)elem_dofs[j_b]; 
                    // Block Indexing: GlobalDoF_j + h_c * N
                    // Ensure type promotion to double before addition to avoid overflow (unlikely but safe)
                    double gcol = gdof_j + (double)h_c * (double)numGlobalDofs;
                    
                    for(int r=0; r<blockSize; r++) {
                        int i_b = r%6, h_r = r/6;
                        double gdof_i = (double)elem_dofs[i_b];
                        double grow = gdof_i + (double)h_r * (double)numGlobalDofs;
                        
                        outI[j_offset+count] = grow;
                        outJ[j_offset+count] = gcol;
                        outV_r[j_offset+count] = Ke_r[r+c*blockSize];
                        outV_i[j_offset+count] = Ke_i[r+c*blockSize];
                        count++;
                    }
                }
            }
        }
    }

    if (calcJ) {
        mxSetM(plhs[0], global_nnz); mxSetM(plhs[1], global_nnz); mxSetM(plhs[2], global_nnz);
    }
    mxSetM(plhs[3], global_res); mxSetM(plhs[4], global_res); mxSetM(plhs[5], global_res); mxSetM(plhs[6], global_res);
}