#ifndef MEX_MATERIAL_UTILS_HPP
#define MEX_MATERIAL_UTILS_HPP

#include <cmath>
#include <vector>
#include <algorithm>

// 材料数据结构 (SoA Layout) - 保持不变
// 注意: 这里的 coefs 依然被假定为每段 4 个系数 (Cubic layout)
// 对应 MATLAB 端构造的 [0, 0, slope, intercept] 系数结构
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
    // 这是判断是否需要外推的边界
    double max_s = mat.maxBSq[idx];
    
    // 2. 准备外推逻辑 (Linear Extrapolation Setup)
    // 如果 s > max_s，我们需要基于 max_s 处的状态进行线性外推
    bool is_extrap = (s > max_s);
    
    // 确定用于查表的点: 如果在外推区，强制使用末端点 max_s
    // 这样可以获取到最后一段（通常是真空切线修正段）的斜率和值
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
    // 边界保护: 确保 k 在有效区间 [0, count-1] 内
    if (k >= count) k = count - 1;
    if (k < 0) k = 0;
    
    // 4. 多项式求值 (Cubic - 保持数据结构兼容性)
    // [CRITICAL] 保持 stride = 4，不改变输入输出结构。
    // MATLAB 端已将线性系数补零为 [0, 0, slope, intercept]。
    // 公式会自动退化为线性计算，无需分支判断。
    
    double dx = s_eval - b_ptr[k];
    const double* c = &mat.coefs[(start_idx + k) * 4]; // Stride = 4 (保持不变)
    
    // 值: nu = c0*x^3 + c1*x^2 + c2*x + c3
    // 若 c0=c1=0，则 nu = c2*x + c3 (线性)
    double nu_eval = ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3];
    
    // 导数: dnu/ds = 3*c0*x^2 + 2*c1*x + c2
    // 若 c0=c1=0，则 dnu/ds = c2 (常数斜率)
    double dnu_eval = (3.0 * c[0] * dx + 2.0 * c[1]) * dx + c[2];

    // 5. 应用外推 (Linear Extrapolation)
    if (is_extrap) {
        // [线性外推]
        // 沿最后一点的切线方向延伸
        // nu(s) = nu(max_s) + slope * (s - max_s)
        double delta = s - max_s;
        nu = nu_eval + dnu_eval * delta;
        dnu_ds = dnu_eval; // 导数保持为末端斜率
    } else {
        // [正常插值]
        nu = nu_eval;
        dnu_ds = dnu_eval;
    }
    
    // ============================================================
    // [PHYSICS CONSTRAINT] 物理防御机制
    // ============================================================
    
    // 获取真空磁阻率 nu0 (存储在 linearNu 中)
    double nu0 = mat.linearNu[idx];
    
    // 1. 真空极限截断 (Vacuum Limit)
    // 物理上 mu_r >= 1 => nu <= nu0。
    // 线性外推可能会导致 nu 超过 nu0 (虽然修正后斜率指向 nu0，但数值误差可能微超)，需强制截断。
    // if (nu > nu0) {
    //     nu = nu0;
    //     dnu_ds = 0.0; // 达到极限后视为常数 (真空)，导数归零
    // }
    // // 2. 正值性约束 (Positivity / Stability)
    // // 防止数值错误导致负值
    // else if (nu < 1e-6) {
    //     nu = nu0; // 回退到 nu0 以避免矩阵奇异
    //     dnu_ds = 0.0;
    // }
}

#endif