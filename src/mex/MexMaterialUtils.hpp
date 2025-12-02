#ifndef MEX_MATERIAL_UTILS_HPP
#define MEX_MATERIAL_UTILS_HPP

#include <cmath>
#include <vector>
#include <algorithm>

// 材料数据结构 (SoA Layout)
struct MexMaterialData {
    const double* linearNu;     // [MaxTag+1] (Stores nu0 for Nonlinear materials)
    const double* isNonlinear;  // [MaxTag+1]
    const double* maxBSq;       // [MaxTag+1]
    
    // Spline Data (Flattened)
    const double* breaks;       
    const double* coefs;        
    
    // Indexing
    const double* splineStart;  
    const double* splineCount;  
};

// 类似于 MATLAB 的 ppval，同时计算值和一阶导数
// s = B^2
inline void evaluate_nu_derivative(double s, int tag, const MexMaterialData& mat, 
                                   double& nu, double& dnu_ds) {
    int idx = tag - 1; 
    
    // 1. 线性材料直接返回
    if (mat.isNonlinear[idx] == 0.0) {
        nu = mat.linearNu[idx];
        dnu_ds = 0.0;
        return;
    }

    // 2. 输入截断 (Clamping)
    // 防止超出样条范围导致的数值爆炸
    double max_s = mat.maxBSq[idx];
    bool is_clamped = false;
    if (s > max_s) {
        s = max_s;
        is_clamped = true;
    }
    
    // 3. 定位区间 (Binary Search)
    int start_idx = (int)mat.splineStart[idx];
    int count = (int)mat.splineCount[idx]; 
    const double* b_ptr = &mat.breaks[start_idx];
    
    int k = 0;
    int low = 0, high = count - 1;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (s >= b_ptr[mid]) {
            k = mid;
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    if (k >= count) k = count - 1;
    if (k < 0) k = 0;
    
    // 4. 多项式求值 (Cubic)
    double dx = s - b_ptr[k];
    const double* c = &mat.coefs[(start_idx + k) * 4];
    
    // nu = c1*x^3 + c2*x^2 + c3*x + c4
    nu = ((c[0] * dx + c[1]) * dx + c[2]) * dx + c[3];
    
    // dnu/ds = 3*c1*x^2 + 2*c2*x + c3
    dnu_ds = (3.0 * c[0] * dx + 2.0 * c[1]) * dx + c[2];

    // ============================================================
    // [CRITICAL FIX] 物理防御机制 (Robustness)
    // ============================================================
    
    // 获取真空磁阻率 nu0 (存储在 linearNu 中)
    double nu0 = mat.linearNu[idx];

    // 1. 深度饱和修正 (Deep Saturation Fix)
    // 如果 B 场极大，nu 不应超过 nu0 (空气)。
    // 样条外插或过冲可能导致 nu > nu0。
    if (nu > nu0) {
        nu = nu0;
        dnu_ds = 0.0; // 强制导数为 0，稳定 Newton 迭代
    }
    // 2. 非物理值修正 (Negative Nu Fix)
    // 样条震荡可能导致 nu < 0
    else if (nu < 0.0) {
        nu = nu0; // 回退到 nu0
        dnu_ds = 0.0;
    }
    // 3. 截断区的一致性修正 (Clamping Consistency)
    // 如果 B^2 超出了定义域，且我们截断了输入 s，
    // 意味着 nu 在此区域被视为常数 (nu(max_s))。
    // 因此，物理上 dnu/ds 应当为 0。
    // *注意*: 虽然 MATLAB 代码声称保留导数以平滑过渡，但在高度非线性下，
    // 不一致的 Jacobian (dnu!=0 但 nu 不变) 会导致死循环。
    // 这里我们选择更稳健的策略：进入平坦区后，导数归零。
    else if (is_clamped) {
        // 可选：为了极度稳健，建议将导数设为 0
        dnu_ds = 0.0; 
        
        // 鉴于 MATLAB 代码保留了导数，这里为了保持行为一致先不强制归零，
        // 除非上述 nu > nu0 触发 (通常深饱和都会触发 nu > nu0)。
    }
}

#endif