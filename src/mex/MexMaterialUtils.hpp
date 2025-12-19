#ifndef MEX_MATERIAL_UTILS_HPP
#define MEX_MATERIAL_UTILS_HPP

#include <cmath>
#include <vector>
#include <algorithm>

// 材料数据结构 (SoA Layout)
struct MexMaterialData {
    const double* linearNu;     
    const double* isNonlinear;  
    const double* maxBSq;       
    const double* breaks;       
    const double* coefs;        
    const double* splineStart;  
    const double* splineCount;  
    const double* splineCoeStart; // [新增] 专门用于系数寻址，修复偏移 Bug
};

/**
 * 计算磁阻率及其对 B^2 的导数
 * s = B^2
 */
inline void evaluate_nu_derivative(double s, int tag, const MexMaterialData& mat, 
                                   double& nu, double& dnu_ds) {
    // 材质索引通常为 Tag - 1
    int idx = tag - 1; 
    
    // 1. 线性材料处理
    if (mat.isNonlinear[idx] == 0.0) {
        nu = mat.linearNu[idx];
        dnu_ds = 0.0;
        return;
    }
    
    // 2. 获取定义域与外推检查
    double max_s = mat.maxBSq[idx];
    bool is_extrap = (s > max_s);
    double s_eval = is_extrap ? max_s : s;

    // 3. 定位区间 (Binary Search)
    int start_idx = (int)mat.splineStart[idx];
    int coe_start_idx = (int)mat.splineCoeStart[idx]; // [关键修复]
    int count = (int)mat.splineCount[idx]; 
    const double* b_ptr = &mat.breaks[start_idx];
    
    int k = 0;
    int low = 0, high = count - 1;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (s_eval >= b_ptr[mid]) {
            k = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    if (k >= count) k = count - 1;
    if (k < 0) k = 0;
    
    // 4. 多项式求值 (Cubic PCHIP)
    // [关键修复]: 使用 coe_start_idx 而非 start_idx * 4，解决多材质偏移累积问题
    double dx = s_eval - b_ptr[k];
    const double* c = &mat.coefs[coe_start_idx + k * 4]; 
    
    double nu_eval = ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3];
    double dnu_eval = (3.0 * c[0] * dx + 2.0 * c[1]) * dx + c[2];

    // 5. 应用外推与物理约束
    if (is_extrap) {
        nu = nu_eval + dnu_eval * (s - max_s);
        dnu_ds = dnu_eval; 
    } else {
        nu = nu_eval;
        dnu_ds = dnu_eval;
    }
    
    // 真空极限截断保护
    double nu0 = mat.linearNu[idx];
    if (nu > nu0) {
        nu = nu0;
        dnu_ds = 0.0;
    } else if (nu < 1e-7) {
        nu = nu0;
        dnu_ds = 0.0;
    }
}

#endif