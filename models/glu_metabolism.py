#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
谷氨酸代谢模块 (Glutamate Metabolism Module)

功能:
- 接收上游T7聚合酶的活性作为输入。
- 模拟谷氨酸在细胞内的生产和向细胞外的分泌过程。
- 模型基于常微分方程(ODE)，描述了细胞内和细胞外谷氨酸浓度的动态变化。

如何使用:
- 实例化 `GluMetabolismModel`。
- 调用 `simulate` 方法，并提供T7聚合酶的活性，以运行模拟。
"""

import numpy as np
from scipy.integrate import odeint

class GluMetabolismModel:
    """
    谷氨酸生产和分泌的ODE模型，包含信号转导与关键酶表达动力学。
    
    该模型将T7聚合酶的活性与谷氨酸合成速率关联起来，并模拟
    谷氨酸通过分泌系统(如一个简单的转运蛋白)从细胞内运输到细胞外的过程。
    """
    def __init__(self, **params):
        """
        初始化模型参数。
        
        Parameters:
        -----------
        k_prod_max : float, optional
            最大谷氨酸生产速率 (mM/hr), 默认为 50.0。
        K_t7 : float, optional
            T7活性达到半最大生产速率时的值 (AU), 默认为 500.0。
        k_export_max : float, optional
            最大谷氨酸分泌速率 (mM/hr), 默认为 100.0。
        K_export : float, optional
            细胞内谷氨酸浓度达到半最大分泌速率时的值 (mM), 默认为 10.0。
        k_dilution : float, optional
            细胞生长或降解导致的稀释/降解速率 (1/hr), 默认为 0.1。
        V_intra_over_V_extra : float, optional
            细胞内总体积与细胞外总体积的比率, 默认为 0.01。
        k_syn_icd : float, optional
            Icd合成速率 (1/hr), 默认为 2.0。
        k_syn_gdhA : float, optional
            gdhA合成速率 (1/hr), 默认为 2.0。
        k_deg_icd : float, optional
            Icd降解速率 (1/hr), 默认为 0.2。
        k_deg_gdhA : float, optional
            gdhA降解速率 (1/hr), 默认为 0.2。
        Vmax_icd : float, optional
            Icd最大催化速率 (mM/hr), 默认为 100.0。
        K_icd : float, optional
            Icd底物常数 (mM), 默认为 5.0。
        Vmax_gdhA : float, optional
            gdhA最大催化速率 (mM/hr), 默认为 100.0。
        K_gdhA : float, optional
            gdhA底物常数 (mM), 默认为 5.0。
        n_hill : float, optional
            Hill系数，用于增强开关效应, 默认为 4.0。
        """
        self.k_prod_max = params.get('k_prod_max', 50.0)
        self.K_t7 = params.get('K_t7', 500.0)
        self.k_export_max = params.get('k_export_max', 100.0)
        self.K_export = params.get('K_export', 10.0)
        self.k_dilution = params.get('k_dilution', 0.1)
        self.V_ratio = params.get('V_intra_over_V_extra', 0.01) # V_intra / V_extra
        self.k_syn_icd = params.get('k_syn_icd', 2.0)      # Icd合成速率 (1/hr)
        self.k_syn_gdhA = params.get('k_syn_gdhA', 2.0)    # gdhA合成速率 (1/hr)
        self.k_deg_icd = params.get('k_deg_icd', 0.2)      # Icd降解速率 (1/hr)
        self.k_deg_gdhA = params.get('k_deg_gdhA', 0.2)    # gdhA降解速率 (1/hr)
        self.Vmax_icd = params.get('Vmax_icd', 100.0)      # Icd最大催化速率 (mM/hr)
        self.K_icd = params.get('K_icd', 5.0)              # Icd底物常数 (mM)
        self.Vmax_gdhA = params.get('Vmax_gdhA', 100.0)    # gdhA最大催化速率 (mM/hr)
        self.K_gdhA = params.get('K_gdhA', 5.0)            # gdhA底物常数 (mM)
        self.n_hill = params.get('n_hill', 4.0)            # Hill系数，用于增强开关效应

    def dydt(self, y, t, t7_activity):
        """
        定义谷氨酸代谢的常微分方程组。

        y: array
            状态变量 [Glu_intra, Glu_extra, Icd, gdhA]
            - Glu_intra: 细胞内谷氨酸浓度 (mM)
            - Glu_extra: 细胞外谷氨酸浓度 (mM)
            - Icd: Icd表达水平
            - gdhA: gdhA表达水平
        t: float
            时间
        t7_activity: float
            T7聚合酶的活性 (AU)
        """
        Glu_intra, Glu_extra, Icd, gdhA = y
        
        # 信号转导: T7驱动Icd/gdhA表达 (使用Hill函数实现开关效应)
        t7_signal = (t7_activity**self.n_hill) / (self.K_t7**self.n_hill + t7_activity**self.n_hill)
        dIcd_dt = self.k_syn_icd * t7_signal - self.k_deg_icd * Icd
        dgdhA_dt = self.k_syn_gdhA * t7_signal - self.k_deg_gdhA * gdhA
        
        # 代谢通路: 酶促反应速率简化为与酶浓度成正比
        v_prod = self.Vmax_gdhA * gdhA
        
        # 2. 谷氨酸分泌速率 (从细胞内到细胞外，Michaelis-Menten形式)
        v_export = self.k_export_max * Glu_intra / (self.K_export + Glu_intra)
        
        # 谷氨酸浓度变化
        dGlu_intra_dt = v_prod - v_export - self.k_dilution * Glu_intra
        dGlu_extra_dt = v_export * self.V_ratio - self.k_dilution * Glu_extra
        
        return [dGlu_intra_dt, dGlu_extra_dt, dIcd_dt, dgdhA_dt]

    def simulate(self, t7_activity, t_end=24.0, dt=0.1):
        """
        运行ODE模拟。

        Parameters:
        -----------
        t7_activity : float
            恒定的T7聚合酶活性 (AU)。
        t_end : float, optional
            模拟结束时间 (小时), 默认为 24.0。
        dt : float, optional
            模拟时间步长 (小时), 默认为 0.1。

        Returns:
        --------
        tuple: (t, solution)
            - t: 时间点数组
            - solution: 状态变量 [Glu_intra, Glu_extra] 在每个时间点的解
        """
        t = np.arange(0, t_end, dt)
        y0 = [0.0, 0.0, 0.0, 0.0]  # 初始: Glu_intra, Glu_extra, Icd, gdhA
        
        solution = odeint(self.dydt, y0, t, args=(t7_activity,))
        
        return t, solution

# --- 示例 ---
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import os

    # 解决中文字体显示问题
    try:
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
        plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号
    except Exception as e:
        print(f"设置中文字体失败: {e}")

    print("="*20 + " 谷氨酸代谢模块优化测试 " + "="*20)
    
    # --- 1. 定义优化后的参数 ---
    optimized_params = {
        'K_t7': 800.0,           # 降低激活阈值，使高T7(3000)能有效激活
        'k_syn_icd': 5.0,        # 提高合成速率，确保酶能积累
        'k_syn_gdhA': 5.0,       # 提高合成速率，确保酶能积累
        'k_deg_icd': 0.3,        # 适度降解速率
        'k_deg_gdhA': 0.3,       # 适度降解速率
        'k_dilution': 0.15,      # 适度稀释速率
        'k_export_max': 100.0,
        'K_export': 5.0,
        'Vmax_gdhA': 50.0,       # 调整最大生产速率，避免过度生产
        'n_hill': 3.0            # 适度的Hill系数
    }
    
    model = GluMetabolismModel(**optimized_params)
    print("模型已使用优化参数创建。")

    # --- 2. 定义高低T7活性 ---
    t7_activity_high = 3000  # 确保高于K_t7
    t7_activity_low = 50     # 确保远低于K_t7
    
    print(f"高T7活性: {t7_activity_high} AU, 低T7活性: {t7_activity_low} AU")

    # --- 3. 运行模拟 ---
    t_high, sol_high = model.simulate(t7_activity=t7_activity_high, t_end=48)
    t_low, sol_low = model.simulate(t7_activity=t7_activity_low, t_end=48)
    print("模拟完成。")

    # --- 4. 提取并分析结果 ---
    final_high = sol_high[:, 1][-1]
    final_low = sol_low[:, 1][-1]
    ratio = final_high / final_low if final_low > 1e-9 else float('inf')

    print("\n--- 结果分析 ---")
    print(f"最终细胞外谷氨酸浓度 (高T7): {final_high:.4f} mM")
    print(f"最终细胞外谷氨酸浓度 (低T7): {final_low:.4f} mM")
    print(f"比值 (高/低): {ratio:.2f}")
    
    # 检查酶表达水平
    icd_high_final = sol_high[-1, 2]
    gdhA_high_final = sol_high[-1, 3]
    icd_low_final = sol_low[-1, 2]
    gdhA_low_final = sol_low[-1, 3]
    
    print(f"\n--- 酶表达水平检查 ---")
    print(f"最终Icd水平 (高T7): {icd_high_final:.6f}")
    print(f"最终Icd水平 (低T7): {icd_low_final:.6f}")
    print(f"最终gdhA水平 (高T7): {gdhA_high_final:.6f}")
    print(f"最终gdhA水平 (低T7): {gdhA_low_final:.6f}")
    
    print("="*50)
