#ifndef MEX_MATERIAL_UTILS_HPP
#define MEX_MATERIAL_UTILS_HPP

#include <cmath>
#include <vector>
#include <algorithm>

// 材料数据结构 (SoA Layout) - 保持不变
// 兼容 MATLAB 端 PCHIP 生成的 [c0, c1, c2, c3] 系数结构 (Stride=4)
struct MexMaterialData {
    const double* linearNu;     // [MaxTag+1] (对于非线性材料，此处存储 nu0)
    const double* isNonlinear;  // [MaxTag+1]
    const double* maxBSq;       // [MaxTag+1]
    
    // Spline Data (Flattened / 扁平化存储)
    const double* breaks;       
    const double* coefs;        
    
    // Indexing
    const double* splineStart;  
    const double* splineCount;  
};

// 类似于 MATLAB 的 ppval，同时计算值和一阶导数
// s = B^2 (磁通密度平方)
inline void evaluate_nu_derivative(double s, int tag, const MexMaterialData& mat, 
                                   double& nu, double& dnu_ds) {
    int idx = tag - 1; 
    
    // 1. 线性材料直接返回
    if (mat.isNonlinear[idx] == 0.0) {
        nu = mat.linearNu[idx];
        dnu_ds = 0.0;
        return;
    }
    
    // 获取定义域上限 maxBSq (样条定义的最后一个断点)
    double max_s = mat.maxBSq[idx];
    
    // 2. 准备外推逻辑 (Linear Extrapolation Setup)
    // 如果 s > max_s，我们需要基于 max_s 处的状态进行线性外推
    // PCHIP 保证了末端导数的连续性，因此线性外推是光滑的 (C1)
    bool is_extrap = (s > max_s);
    
    // 确定用于查表的点
    double s_eval = is_extrap ? max_s : s;

    // 3. 定位区间 (Binary Search)
    int start_idx = (int)mat.splineStart[idx];
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
    // 边界保护
    if (k >= count) k = count - 1;
    if (k < 0) k = 0;
    
    // 4. 多项式求值 (Cubic PCHIP Evaluation)
    // MATLAB PCHIP 系数顺序为 [c0, c1, c2, c3] 对应 c0*x^3 + ...
    
    double dx = s_eval - b_ptr[k];
    const double* c = &mat.coefs[(start_idx + k) * 4]; // Stride = 4
    
    // 值: nu = c0*x^3 + c1*x^2 + c2*x + c3
    double nu_eval = ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3];
    
    // 导数: dnu/ds = 3*c0*x^2 + 2*c1*x + c2
    double dnu_eval = (3.0 * c[0] * dx + 2.0 * c[1]) * dx + c[2];

    // 5. 应用外推 (Linear Extrapolation)
    if (is_extrap) {
        // 沿最后一点的切线方向延伸
        // 由于 MATLAB 端已做过 PCHIP 平滑，这里的 dnu_eval 是准确的切线斜率
        double delta = s - max_s;
        nu = nu_eval + dnu_eval * delta;
        dnu_ds = dnu_eval; 
    } else {
        // [正常插值]
        nu = nu_eval;
        dnu_ds = dnu_eval;
    }
    
    // ============================================================
    // [PHYSICS CONSTRAINT] 物理防御机制 (已启用)
    // 对于高阶 SDC 算法，防止数值超调至关重要
    // ============================================================
    
    // 获取真空磁阻率 nu0
    double nu0 = mat.linearNu[idx];
    
    // 1. 真空极限截断 (Vacuum Limit)
    // 物理上 nu <= nu0。如果外推导致略微超过 nu0，强制截断。
    if (nu > nu0) {
        nu = nu0;
        dnu_ds = 0.0; // 达到真空极限，视为常数，导数为 0
    }
    // 2. 正值性约束 (Positivity / Stability)
    // 防止数值错误导致负值或接近 0 (这会导致刚度矩阵奇异)
    else if (nu < 1e-7) {
        nu = nu0; // 回退到 nu0 以保证求解器继续运行 (比报错 NaN 更好)
        dnu_ds = 0.0;
    }
}

#endif