/**
 * MexElemUtils.hpp
 * 通用有限元数学工具库 (v1.1 - Added Inverse)
 */

#ifndef MEX_ELEM_UTILS_HPP
#define MEX_ELEM_UTILS_HPP

#include "mex.h"
#include <cmath>
#include <algorithm>
#include <vector>

struct Mat3x3 {
    double data[3][3];
    
    // y = M * x
    void multVec(const double* x, double* y) const {
        for(int i=0; i<3; i++) {
            y[i] = data[i][0]*x[0] + data[i][1]*x[1] + data[i][2]*x[2];
        }
    }

    // y = M^T * x
    void multVecTrans(const double* x, double* y) const {
        for(int i=0; i<3; i++) {
            y[i] = data[0][i]*x[0] + data[1][i]*x[1] + data[2][i]*x[2];
        }
    }
};

class MexElemUtils {
public:
    // 计算 J = [x2-x1, ...] 并返回 detJ
    static inline double compute_jacobian_3d(const double* P_elem, Mat3x3& J) {
        double x1 = P_elem[0]; double y1 = P_elem[1]; double z1 = P_elem[2];
        double x2 = P_elem[3]; double y2 = P_elem[4]; double z2 = P_elem[5];
        double x3 = P_elem[6]; double y3 = P_elem[7]; double z3 = P_elem[8];
        double x4 = P_elem[9]; double y4 = P_elem[10]; double z4 = P_elem[11];

        J.data[0][0] = x2 - x1; J.data[1][0] = y2 - y1; J.data[2][0] = z2 - z1;
        J.data[0][1] = x3 - x1; J.data[1][1] = y3 - y1; J.data[2][1] = z3 - z1;
        J.data[0][2] = x4 - x1; J.data[1][2] = y4 - y1; J.data[2][2] = z4 - z1;

        double det = J.data[0][0] * (J.data[1][1]*J.data[2][2] - J.data[1][2]*J.data[2][1])
                   - J.data[0][1] * (J.data[1][0]*J.data[2][2] - J.data[1][2]*J.data[2][0])
                   + J.data[0][2] * (J.data[1][0]*J.data[2][1] - J.data[1][1]*J.data[2][0]);
        return det;
    }

    // 计算 J 的逆矩阵 J_inv
    static inline void compute_inverse_3x3(const Mat3x3& J, double det, Mat3x3& J_inv) {
        double invDet = 1.0 / det;
        // Row 1
        J_inv.data[0][0] = (J.data[1][1]*J.data[2][2] - J.data[1][2]*J.data[2][1]) * invDet;
        J_inv.data[0][1] = (J.data[0][2]*J.data[2][1] - J.data[0][1]*J.data[2][2]) * invDet;
        J_inv.data[0][2] = (J.data[0][1]*J.data[1][2] - J.data[0][2]*J.data[1][1]) * invDet;
        // Row 2
        J_inv.data[1][0] = (J.data[1][2]*J.data[2][0] - J.data[1][0]*J.data[2][2]) * invDet;
        J_inv.data[1][1] = (J.data[0][0]*J.data[2][2] - J.data[0][2]*J.data[2][0]) * invDet;
        J_inv.data[1][2] = (J.data[0][2]*J.data[1][0] - J.data[0][0]*J.data[1][2]) * invDet;
        // Row 3
        J_inv.data[2][0] = (J.data[1][0]*J.data[2][1] - J.data[1][1]*J.data[2][0]) * invDet;
        J_inv.data[2][1] = (J.data[0][1]*J.data[2][0] - J.data[0][0]*J.data[2][1]) * invDet;
        J_inv.data[2][2] = (J.data[0][0]*J.data[1][1] - J.data[0][1]*J.data[1][0]) * invDet;
    }
};

#endif