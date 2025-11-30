#include "mex.h"
#include <vector>
#include <cmath>
#include <omp.h>

// 简单的 3x3 矩阵乘法与行列式计算辅助函数
// 注意：在实际工程中可以使用 Eigen 库，但为了不依赖第三方库，这里手动展开

void compute_jacobian(const double* P, const double* T_col, int n_nodes, 
                      double* J_mat, double& detJ, double* invDetJ_J_trans) {
    // 获取四面体 4 个顶点的坐标
    // T_col 是当前单元的节点索引 (1-based from MATLAB)
    int n1 = (int)T_col[0] - 1;
    int n2 = (int)T_col[1] - 1;
    int n3 = (int)T_col[2] - 1;
    int n4 = (int)T_col[3] - 1;

    double x1 = P[3*n1], y1 = P[3*n1+1], z1 = P[3*n1+2];
    double x2 = P[3*n2], y2 = P[3*n2+1], z2 = P[3*n2+2];
    double x3 = P[3*n3], y3 = P[3*n3+1], z3 = P[3*n3+2];
    double x4 = P[3*n4], y4 = P[3*n4+1], z4 = P[3*n4+2];

    // d1 = x2 - x1, d2 = x3 - x1, d3 = x4 - x1
    double d1x = x2-x1, d1y = y2-y1, d1z = z2-z1;
    double d2x = x3-x1, d2y = y3-y1, d2z = z3-z1;
    double d3x = x4-x1, d3y = y4-y1, d3z = z4-z1;

    // Jacobian Matrix J = [d1, d2, d3]
    // detJ
    detJ = d1x*(d2y*d3z - d2z*d3y) - d1y*(d2x*d3z - d2z*d3x) + d1z*(d2x*d3y - d2y*d3x);
    double invDetJ = 1.0 / detJ;

    // Piola 变换核心: inv(J) * detJ (即伴随矩阵的转置)
    // 用于 H(curl) 变换: C_phy = (1/detJ) * J * C_ref
    // 这里我们直接计算变换矩阵 M = J * (1/detJ) ??? 
    // 等等，您的 MATLAB 代码是: C_phy = (C_ref * J_mat') * invDetJ
    // 这意味着 C_phy_row = C_ref_row * J^T * (1/detJ)
    
    // 我们在这里预计算 J * (1/detJ) 这一项并不划算，
    // 因为 C_ref 是多个积分点。不如直接存 J_mat。
    
    J_mat[0] = d1x; J_mat[1] = d2x; J_mat[2] = d3x;
    J_mat[3] = d1y; J_mat[4] = d2y; J_mat[5] = d3y;
    J_mat[6] = d1z; J_mat[7] = d2z; J_mat[8] = d3z;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // 输入检查
    if (nrhs < 7) mexErrMsgIdAndTxt("MagPackage:Assemble", "Need 7 inputs: P, T, Dofs, Signs, Nu, Qw, CurlRef");

    // 1. 解析输入数据 (使用指针直接访问，零拷贝)
    double *P = mxGetPr(prhs[0]); // [3 x Np]
    size_t n_nodes = mxGetN(prhs[0]);
    
    double *T = mxGetPr(prhs[1]); // [4 x Ne]
    size_t n_elems = mxGetN(prhs[1]);
    
    double *Dofs = mxGetPr(prhs[2]); // [6 x Ne]
    double *Signs = mxGetPr(prhs[3]); // [6 x Ne]
    double *Nu = mxGetPr(prhs[4]); // [Ne x 1] or Scalar
    bool nu_is_scalar = (mxGetM(prhs[4]) * mxGetN(prhs[4]) == 1);
    
    double *Qw = mxGetPr(prhs[5]); // [Nq x 1]
    size_t n_q = mxGetM(prhs[5]) * mxGetN(prhs[5]);
    
    double *CurlRef = mxGetPr(prhs[6]); // [6 x 3 x Nq] - 务必注意维度顺序
    // MATLAB 存储顺序为 (6, 3, Nq)，即前 6 个是第一个积分点第一个分量...

    // 2. 准备输出 (Triplets)
    // 每个单元产生 36 个非零元
    size_t n_triplets = n_elems * 36;
    plhs[0] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL); // I
    plhs[1] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL); // J
    plhs[2] = mxCreateDoubleMatrix(n_triplets, 1, mxREAL); // V
    
    double *I_out = mxGetPr(plhs[0]);
    double *J_out = mxGetPr(plhs[1]);
    double *V_out = mxGetPr(plhs[2]);

    // 3. 并行计算 (OpenMP)
    // 确保在 MATLAB 编译时开启 OpenMP (mex -v CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" ...)
    
    #pragma omp parallel for schedule(static)
    for (size_t e = 0; e < n_elems; ++e) {
        
        // --- 几何计算 ---
        double T_local[4];
        for(int k=0; k<4; k++) T_local[k] = T[k + e*4]; // 读取 T 的第 e 列
        
        double J_mat[9];
        double detJ;
        double dummy[9];
        compute_jacobian(P, T_local, n_nodes, J_mat, detJ, dummy);
        double absDetJ = std::abs(detJ);
        double invDetJ = 1.0 / detJ;
        
        double nu_e = nu_is_scalar ? Nu[0] : Nu[e];
        
        // --- 单元刚度矩阵 Ke (6x6) ---
        double Ke[36] = {0}; // Row-major or Col-major? We will flatten to Col-major later
        
        for (int q = 0; q < n_q; ++q) {
            double w = Qw[q] * absDetJ * nu_e;
            
            // 读取当前积分点的 CurlRef (6x3)
            // CurlRef index: (row, col, q) -> row + 6*col + 18*q
            double C_phy[18]; // 6x3
            
            for (int i = 0; i < 6; ++i) {
                // C_ref row i: [cx, cy, cz]
                double cr_x = CurlRef[i + 0*6 + q*18];
                double cr_y = CurlRef[i + 1*6 + q*18];
                double cr_z = CurlRef[i + 2*6 + q*18];
                
                // MATLAB logic: C_phy = (C_ref * J_mat') * invDetJ
                // Row_phy = Row_ref * J^T * invDetJ
                //         = [cr_x, cr_y, cr_z] * [d1x d1y d1z; d2x...; d3x...]^T * invDetJ
                //         = [cr_x, cr_y, cr_z] * [d1x d2x d3x; d1y...; ...] * invDetJ
                
                double cp_x = (cr_x * J_mat[0] + cr_y * J_mat[1] + cr_z * J_mat[2]) * invDetJ;
                double cp_y = (cr_x * J_mat[3] + cr_y * J_mat[4] + cr_z * J_mat[5]) * invDetJ;
                double cp_z = (cr_x * J_mat[6] + cr_y * J_mat[7] + cr_z * J_mat[8]) * invDetJ;
                
                C_phy[i + 0*6] = cp_x;
                C_phy[i + 1*6] = cp_y;
                C_phy[i + 2*6] = cp_z;
            }
            
            // 累加 Ke += w * C_phy * C_phy'
            // Ke(i, j) += w * dot(C_phy_i, C_phy_j)
            for (int c = 0; c < 6; ++c) {     // col
                for (int r = 0; r < 6; ++r) { // row
                    double dot_val = C_phy[r + 0*6] * C_phy[c + 0*6] + 
                                     C_phy[r + 1*6] * C_phy[c + 1*6] + 
                                     C_phy[r + 2*6] * C_phy[c + 2*6];
                    Ke[r + c*6] += w * dot_val; // Column-major storage in Ke array
                }
            }
        }
        
        // --- 符号修正与填充 ---
        // Signs: [6 x Ne]
        double s[6];
        for(int k=0; k<6; k++) s[k] = Signs[k + e*6];
        
        double dofs[6];
        for(int k=0; k<6; k++) dofs[k] = Dofs[k + e*6];
        
        size_t offset = e * 36;
        int idx = 0;
        
        // 按照 MATLAB sparse 的要求填充: 列优先
        for (int c = 0; c < 6; ++c) {
            for (int r = 0; r < 6; ++r) {
                // Apply signs: Ke_new(r,c) = Ke(r,c) * s(r) * s(c)
                double val = Ke[r + c*6] * s[r] * s[c];
                
                I_out[offset + idx] = dofs[r];
                J_out[offset + idx] = dofs[c];
                V_out[offset + idx] = val;
                idx++;
            }
        }
    }
}