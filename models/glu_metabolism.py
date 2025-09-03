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
    谷氨酸生产和分泌的ODE模型。
    
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
        """
        self.k_prod_max = params.get('k_prod_max', 50.0)
        self.K_t7 = params.get('K_t7', 500.0)
        self.k_export_max = params.get('k_export_max', 100.0)
        self.K_export = params.get('K_export', 10.0)
        self.k_dilution = params.get('k_dilution', 0.1)
        self.V_ratio = params.get('V_intra_over_V_extra', 0.01) # V_intra / V_extra

    def dydt(self, y, t, t7_activity):
        """
        定义谷氨酸代谢的常微分方程组。

        y: array
            状态变量 [Glu_intra, Glu_extra]
            - Glu_intra: 细胞内谷氨酸浓度 (mM)
            - Glu_extra: 细胞外谷氨酸浓度 (mM)
        t: float
            时间
        t7_activity: float
            T7聚合酶的活性 (AU)
        """
        Glu_intra, Glu_extra = y
        
        # 1. 谷氨酸生产速率 (由T7活性驱动，Michaelis-Menten形式)
        v_prod = self.k_prod_max * t7_activity / (self.K_t7 + t7_activity)
        
        # 2. 谷氨酸分泌速率 (从细胞内到细胞外，Michaelis-Menten形式)
        v_export = self.k_export_max * Glu_intra / (self.K_export + Glu_intra)
        
        # 细胞内谷氨酸浓度变化
        dGlu_intra_dt = v_prod - v_export - self.k_dilution * Glu_intra
        
        # 细胞外谷氨酸浓度变化 (考虑体积比)
        dGlu_extra_dt = v_export * self.V_ratio - self.k_dilution * Glu_extra
        
        return [dGlu_intra_dt, dGlu_extra_dt]

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
        y0 = [0.0, 0.0]  # 初始条件：细胞内外谷氨酸浓度均为0
        
        solution = odeint(self.dydt, y0, t, args=(t7_activity,))
        
        return t, solution

# --- 示例 ---
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    print("="*20 + " 谷氨酸代谢模块测试 " + "="*20)
    
    model = GluMetabolismModel()
    
    # --- 条件1: 高T7活性 (AND门开启) ---
    t7_activity_high = 1000.0
    t_high, sol_high = model.simulate(t7_activity_high, t_end=48)
    
    # --- 条件2: 低T7活性 (AND门关闭) ---
    t7_activity_low = 50.0
    t_low, sol_low = model.simulate(t7_activity_low, t_end=48)
    
    # --- 绘图比较 ---
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # 高活性图
    ax1.plot(t_high, sol_high[:, 0], label='Intracellular Glu (High T7)', color='blue')
    ax1.plot(t_high, sol_high[:, 1], label='Extracellular Glu (High T7)', color='red', linestyle='--')
    ax1.set_title(f'Glutamate Dynamics (T7 Activity = {t7_activity_high} AU)')
    ax1.set_ylabel('Concentration (mM)')
    ax1.legend()
    ax1.grid(True)
    
    # 低活性图
    ax2.plot(t_low, sol_low[:, 0], label='Intracellular Glu (Low T7)', color='navy')
    ax2.plot(t_low, sol_low[:, 1], label='Extracellular Glu (Low T7)', color='indianred', linestyle='--')
    ax2.set_title(f'Glutamate Dynamics (T7 Activity = {t7_activity_low} AU)')
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Concentration (mM)')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.show()
    
    print("\n测试完成。图像显示了在高/低T7活性下，细胞内和细胞外谷氨酸浓度的变化。")
    print("="*54)
