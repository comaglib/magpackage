/**
 * MexMaterialUtils.hpp (Fixed V7 - Failsafe)
 * 功能: 处理非线性材料属性
 * * 修复日志:
 * 1. [Safety] 增加 Nu 值物理完整性检查。
 * 如果计算出的 nu < 1e-9 (非物理/数据错误)，强制重置为真空磁阻率 (795774.7)。
 * 这能防止因材料数据全0导致的 "Line search stuck" 和 "Singular Matrix"。
 */

#ifndef MEX_MATERIAL_UTILS_HPP
#define MEX_MATERIAL_UTILS_HPP

#include <cmath>
#include <algorithm>
#include <vector>

// Vacuum Reluctivity = 1 / (4 * pi * 1e-7) approx 795774.715
#define NU_AIR 795774.7154594767

class MaterialEvaluator {
private:
    const double* m_lin;   
    const double* m_isNon; 
    const double* m_maxB;  
    const double* m_brk;   
    const double* m_coe;   
    const double* m_start; 
    const double* m_cnt;   
    const double* m_coe_start;
    
    int m_maxTag;

public:
    MaterialEvaluator(const double* lin, const double* isNon, const double* maxB,
                      const double* brk, const double* coe, const double* start, 
                      const double* cnt, const double* coe_start, int maxTag)
        : m_lin(lin), m_isNon(isNon), m_maxB(maxB), m_brk(brk), 
          m_coe(coe), m_start(start), m_cnt(cnt), m_coe_start(coe_start), m_maxTag(maxTag) 
    {}

    inline double evaluate(int tag, double b_sq, double& dnu_val) const {
        double nu_val = NU_AIR;
        dnu_val = 0.0;

        // 1. Boundary Check
        if (tag <= 0 || tag > m_maxTag) return nu_val;

        // 2. Linear Case
        if (m_isNon[tag] == 0.0) {
            nu_val = m_lin[tag];
            // Safety: Linear nu shouldn't be 0 either
            if (nu_val < 1e-9) nu_val = NU_AIR;
            return nu_val;
        }

        // 3. Nonlinear Case
        double x = b_sq;
        if (x < 0) x = 0;
        if (x > m_maxB[tag]) x = m_maxB[tag];

        int start_idx = (int)m_start[tag]; 
        int num_seg = (int)m_cnt[tag];
        
        const double* my_brks = &m_brk[start_idx];
        const double* it = std::upper_bound(my_brks, my_brks + num_seg + 1, x);
        int k = (int)(it - my_brks) - 1;
        
        if (k < 0) k = 0;
        if (k >= num_seg) k = num_seg - 1;

        double t = x - my_brks[k];
        
        // Use explicit offset passed from MATLAB
        long coe_base = (long)m_coe_start[tag];
        long coe_idx = coe_base + k * 4;
        
        double c1 = m_coe[coe_idx + 0];
        double c2 = m_coe[coe_idx + 1];
        double c3 = m_coe[coe_idx + 2];
        double c4 = m_coe[coe_idx + 3];

        // Evaluate Cubic
        nu_val = ((c1 * t + c2) * t + c3) * t + c4;
        dnu_val = (3.0 * c1 * t + 2.0 * c2) * t + c3;
        
        // [FAILSAFE CHECK]
        // 如果计算结果非物理 (例如数据全0导致 nu=0)，强制回退到空气
        if (nu_val < 1e-9) {
            nu_val = NU_AIR;
            dnu_val = 0.0; 
        }
        
        return nu_val;
    }
};

#endif