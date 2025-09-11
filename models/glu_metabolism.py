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

import sys
import os
import numpy as np
from scipy.integrate import odeint

# 获取当前文件所在目录的父目录（即项目根目录）
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '..'))

# 将项目根目录添加到Python的模块搜索路径中
sys.path.insert(0, project_root)

from config_manager import ConfigManager # 导入ConfigManager

class GluMetabolismModel:
    def __init__(self, **params):
        # 如果没有传入参数，则从ConfigManager获取默认参数
        if not params:
            config_manager = ConfigManager()
            params = config_manager.get_params('glu_metabolism')

        self.k_prod_max = params.get('k_prod_max', 50.0)
        self.K_t7 = params.get('K_t7', 500.0)
        self.k_export_max = params.get('k_export_max', 100.0)
        self.K_export = params.get('K_export', 10.0)
        self.k_dilution = params.get('k_dilution', 0.1)
        self.V_ratio = params.get('V_intra_over_V_extra', 0.01) # 默认值，可在dydt中动态更新

        # Icd (Isocitrate dehydrogenase) parameters
        self.k_syn_icd = params.get('k_syn_icd', 2.0)
        self.k_deg_icd = params.get('k_deg_icd', 0.2)
        self.Vmax_icd = params.get('Vmax_icd', 100.0)
        self.K_icd = params.get('K_icd', 5.0)

        # gdhA (Glutamate dehydrogenase) parameters
        self.k_syn_gdhA = params.get('k_syn_gdhA', 2.0)
        self.k_deg_gdhA = params.get('k_deg_gdhA', 0.2)
        self.Vmax_gdhA = params.get('Vmax_gdhA', 100.0)
        self.K_gdhA = params.get('K_gdhA', 5.0)
        self.n_hill = params.get('n_hill', 4.0)

    def dydt(self, y, t, t7_activity):
        Glu_intra, Glu_extra, Icd, gdhA = y

        # T7 活性对谷氨酸生产的影响 (Hill function)
        t7_effect = (t7_activity**self.n_hill) / (self.K_t7**self.n_hill + t7_activity**self.n_hill)
        glu_prod_rate = self.k_prod_max * t7_effect

        # Icd and gdhA dynamics
        dIcd_dt = self.k_syn_icd * t7_effect - self.k_deg_icd * Icd
        dgdhA_dt = self.k_syn_gdhA * t7_effect - self.k_deg_gdhA * gdhA

        # Glutamate consumption by Icd and gdhA
        glu_consumption_icd = self.Vmax_icd * Glu_intra / (self.K_icd + Glu_intra) * Icd
        glu_consumption_gdhA = self.Vmax_gdhA * Glu_intra / (self.K_gdhA + Glu_intra) * gdhA

        # Glutamate export
        glu_export_rate = self.k_export_max * Glu_intra / (self.K_export + Glu_intra)

        # Intra-cellular glutamate dynamics
        dGlu_intra_dt = glu_prod_rate - glu_consumption_icd - glu_consumption_gdhA - glu_export_rate - self.k_dilution * Glu_intra

        # Extra-cellular glutamate dynamics (assuming V_ratio = V_intra / V_extra)
        dGlu_extra_dt = (glu_export_rate - self.k_dilution * Glu_extra) * self.V_ratio

        return [dGlu_intra_dt, dGlu_extra_dt, dIcd_dt, dgdhA_dt]

    def simulate(self, t_end=100.0, dt=0.1, t7_activity=100.0):
        t = np.arange(0, t_end, dt)
        # Initial conditions: [Glu_intra, Glu_extra, Icd, gdhA]
        y0 = [0.0, 0.0, 0.0, 0.0]
        solution = odeint(self.dydt, y0, t, args=(t7_activity,))
        return t, solution


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    print("Testing GluMetabolismModel...")

    # 使用默认参数初始化模型
    model = GluMetabolismModel()

    # 模拟不同T7活性下的谷氨酸代谢
    t_high_t7, sol_high_t7 = model.simulate(t7_activity=1000.0) # 高T7活性
    t_low_t7, sol_low_t7 = model.simulate(t7_activity=10.0)   # 低T7活性

    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.plot(t_high_t7, sol_high_t7[:, 0], label='Glu_intra (High T7)')
    plt.plot(t_high_t7, sol_high_t7[:, 1], label='Glu_extra (High T7)')
    plt.plot(t_low_t7, sol_low_t7[:, 0], label='Glu_intra (Low T7)', linestyle='--')
    plt.plot(t_low_t7, sol_low_t7[:, 1], label='Glu_extra (Low T7)', linestyle='--')
    plt.xlabel('Time (h)')
    plt.ylabel('Concentration (mM)')
    plt.title('Intra- and Extra-cellular Glutamate')
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(t_high_t7, sol_high_t7[:, 2], label='Icd (High T7)')
    plt.plot(t_high_t7, sol_high_t7[:, 3], label='gdhA (High T7)')
    plt.plot(t_low_t7, sol_low_t7[:, 2], label='Icd (Low T7)', linestyle='--')
    plt.plot(t_low_t7, sol_low_t7[:, 3], label='gdhA (Low T7)', linestyle='--')
    plt.xlabel('Time (h)')
    plt.ylabel('Expression Level')
    plt.title('Icd and gdhA Expression')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    print("GluMetabolismModel test finished.")
