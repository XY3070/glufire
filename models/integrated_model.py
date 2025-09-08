#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
整合治疗模型模块 (Integrated Therapy Model Module)

功能:
- 将AND门、谷氨酸代谢和肿瘤生长/死亡模型整合在一起。
- 模拟在特定环境条件下（温度和氧气），整个系统的治疗效果。
- 模型基于ODE，描述了肿瘤细胞、死亡细胞和细胞外谷氨酸浓度的动态变化。

如何使用:
- 实例化 `IntegratedTherapyModel`。
- 调用 `simulate` 方法，并提供环境条件字典，以运行模拟。
"""

import numpy as np
from scipy.integrate import odeint

# 导入上游模块
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))
from and_gate import SimpleANDGate
from glu_metabolism import GluMetabolismModel

class IntegratedTherapyModel:
    """
    整合的端到端治疗模型。
    
    该模型将环境信号（O2, Temp）转化为T7活性，再转化为谷氨酸生产，
    最终模拟谷氨酸如何影响肿瘤细胞的生长和死亡（铁死亡）。
    """
    def __init__(self, **params):
        """
        初始化整合模型，包括其子模块。
        """
        # 使用优化的谷氨酸代谢参数 - 确保T7活性能有效激活谷氨酸生产
        glu_params = {
            'K_t7': 800.0,          # 降低阈值，确保治疗条件下T7活性(~1200)能有效激活
            'k_syn_icd': 3.0,       # 提高合成速率
            'k_syn_gdhA': 3.0,      # 提高合成速率  
            'k_deg_icd': 0.3,       # 适度降解
            'k_deg_gdhA': 0.3,      # 适度降解
            'k_dilution': 0.15,     # 降低稀释速率
            'Vmax_gdhA': 100.0,     # 适中的最大生产速率
            'n_hill': 4.0,
            # V_ratio 在dydt中动态计算，这里无需设置
        }
        
        self.and_gate = SimpleANDGate()
        self.glu_metabolism = GluMetabolismModel(**glu_params)
        
        # 肿瘤生长和死亡参数
        self.r = params.get('r', 0.01)  # 大幅降低肿瘤生长速率，使其在治疗时间尺度内几乎不增长
        self.K_tumor = params.get('K_tumor', 1e9)  # 肿瘤承载能力 (细胞数)
        
        # 铁死亡参数 - 强效铁死亡以确保治疗效果显著
        self.k_ferroptosis_max = params.get('k_ferroptosis_max', 15.0) # 进一步提高最大铁死亡速率
        self.K_glu = params.get('K_glu', 30.0) # 设置谷氨酸铁死亡阈值为30mM
        self.n_glu = params.get('n_glu', 5.0) # 进一步增加Hill系数，增强开关效应

        # 工程细胞生长参数
        self.r_eng = params.get('r_eng', 0.2) # 工程细胞生长速率
        self.K_eng = params.get('K_eng', 5e8) # 工程细胞承载能力

        # 新增：用于计算谷氨酸浓度的体积参数
        # V_cell: 单个细胞体积 (L/cell), V_tumor_ext: 肿瘤间质液体积 (L)
        self.V_cell = params.get('V_cell', 2e-12) 
        self.V_tumor_ext = params.get('V_tumor_ext', 0.01)

    def get_t7_activity(self, env_conditions):
        """从AND门模块获取T7活性。"""
        return self.and_gate.get_t7_activity(
            env_conditions['O2_percent'],
            env_conditions['Temp_C']
        )

    def dydt(self, y, t, env_conditions):
        """
        定义整合模型的ODE系统。
        
        y: array [N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA]
            - N_tumor: 存活的肿瘤细胞数
            - D_tumor: 死亡的肿瘤细胞数
            - N_eng: 存活的工程细胞数
            - Glu_intra: 工程细胞内的谷氨酸浓度 (mM)
            - Glu_extra: 细胞外谷氨酸浓度 (mM)
            - Icd: Icd表达水平
            - gdhA: gdhA表达水平
        """
        N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA = y
        
        # 1. 从环境条件计算T7活性
        t7_activity = self.get_t7_activity(env_conditions)
        
        # 2. 调用谷氨酸代谢模块的ODE
        y_glu = [Glu_intra, Glu_extra, Icd, gdhA]
        
        # 动态计算体积比 V_ratio = (N_eng * V_cell) / V_tumor_ext
        current_V_ratio = (N_eng * self.V_cell) / self.V_tumor_ext
        self.glu_metabolism.V_ratio = max(current_V_ratio, 1e-12) # 防止除零

        dGlu_intra_dt, dGlu_extra_dt, dIcd_dt, dgdhA_dt = self.glu_metabolism.dydt(y_glu, t, t7_activity)

        # 3. 肿瘤细胞生长和死亡
        growth_rate_tumor = self.r * N_tumor * (1 - (N_tumor + N_eng) / self.K_tumor)
        ferroptosis_rate = self.k_ferroptosis_max * (Glu_extra**self.n_glu) / (self.K_glu**self.n_glu + Glu_extra**self.n_glu)
        death_term_tumor = ferroptosis_rate * N_tumor
        
        # 4. 工程细胞生长 (依赖T7活性的开关效应)
        t7_hill_factor = (t7_activity**self.glu_metabolism.n_hill) / \
                         (self.glu_metabolism.K_t7**self.glu_metabolism.n_hill + t7_activity**self.glu_metabolism.n_hill)
        growth_rate_eng = self.r_eng * N_eng * (1 - (N_tumor + N_eng) / self.K_tumor) * t7_hill_factor
        
        # 状态变量变化率
        dN_tumor_dt = growth_rate_tumor - death_term_tumor
        dD_tumor_dt = death_term_tumor
        # 工程细胞也受稀释/程序性死亡影响
        dN_eng_dt = growth_rate_eng - self.glu_metabolism.k_dilution * N_eng
        
        # 数值稳定性保护
        if N_tumor < 1.0: N_tumor = 0
        if N_eng < 1.0: N_eng = 0
        dN_tumor_dt = dN_tumor_dt if N_tumor > 0 else 0
        dN_eng_dt = dN_eng_dt if N_eng > 0 else 0
        
        return [dN_tumor_dt, dD_tumor_dt, dN_eng_dt, dGlu_intra_dt, dGlu_extra_dt, dIcd_dt, dgdhA_dt]

    def simulate(self, env_conditions, t_end=100.0, dt=0.5):
        """
        运行整合治疗模拟。
        """
        t = np.arange(0, t_end, dt)
        # 初始条件: [N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA]
        # 降低初始肿瘤细胞数量，使死亡效应更明显
        y0 = [1e6, 0, 5e5, 0, 0, 0, 0]  # 肿瘤细胞1M，工程细胞0.5M
        
        solution = odeint(self.dydt, y0, t, args=(env_conditions,))
        
        return t, solution

# --- 示例 ---
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import os

    # 解决中文字体显示问题
    try:
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
    except Exception as e:
        print(f"设置中文字体失败: {e}")

    print("="*20 + " 整合治疗模型测试 " + "="*20)
    
    model = IntegratedTherapyModel()
    
    # --- 条件1: 治疗开启 (低氧, 高温) ---
    therapy_conditions = {'O2_percent': 1.0, 'Temp_C': 42.0}
    
    # --- 条件2: 对照组 (正常氧, 体温) ---
    control_conditions = {'O2_percent': 21.0, 'Temp_C': 37.0}

    # --- 运行模拟 ---
    print("\n正在运行 '治疗' 条件下的模拟...")
    t_therapy, sol_therapy = model.simulate(env_conditions=therapy_conditions, t_end=200)
    print("正在运行 '对照' 条件下的模拟...")
    t_control, sol_control = model.simulate(env_conditions=control_conditions, t_end=200)
    print("模拟完成。")

    # --- 提取分析数据 ---
    # y: [N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA]
    N_tumor_therapy, D_tumor_therapy, N_eng_therapy = sol_therapy[:, 0], sol_therapy[:, 1], sol_therapy[:, 2]
    Glu_intra_therapy, Glu_extra_therapy = sol_therapy[:, 3], sol_therapy[:, 4]
    Icd_therapy, gdhA_therapy = sol_therapy[:, 5], sol_therapy[:, 6]
    
    N_tumor_control, D_tumor_control, N_eng_control = sol_control[:, 0], sol_control[:, 1], sol_control[:, 2]
    Glu_intra_control, Glu_extra_control = sol_control[:, 3], sol_control[:, 4]
    Icd_control, gdhA_control = sol_control[:, 5], sol_control[:, 6]
    
    # 计算T7活性
    t7_therapy = model.get_t7_activity(therapy_conditions)
    t7_control = model.get_t7_activity(control_conditions)
    
    # 计算治疗效果指标
    tumor_reduction_ratio = N_tumor_therapy[-1] / N_tumor_control[-1]
    glu_ratio = Glu_extra_therapy[-1] / (Glu_extra_control[-1] + 1e-9)
    
    print(f"\n=== 治疗效果分析 ===")
    print(f"T7活性 - 治疗组: {t7_therapy:.1f} AU, 对照组: {t7_control:.1f} AU")
    print(f"谷氨酸浓度 - 治疗组: {Glu_extra_therapy[-1]:.3f} mM, 对照组: {Glu_extra_control[-1]:.3f} mM")
    print(f"谷氨酸比值 (治疗/对照): {glu_ratio:.3f}")
    print(f"肿瘤细胞数量 - 治疗组: {N_tumor_therapy[-1]:.1e}, 对照组: {N_tumor_control[-1]:.1e}")
    print(f"死亡细胞数量 - 治疗组: {D_tumor_therapy[-1]:.1e}, 对照组: {D_tumor_control[-1]:.1e}")
    print(f"肿瘤细胞存活率 (治疗/对照): {tumor_reduction_ratio:.3f}")
    
    # 计算铁死亡速率检查
    ferroptosis_rate_therapy = model.k_ferroptosis_max * (Glu_extra_therapy[-1]**model.n_glu) / (model.K_glu**model.n_glu + Glu_extra_therapy[-1]**model.n_glu)
    ferroptosis_rate_control = model.k_ferroptosis_max * (Glu_extra_control[-1]**model.n_glu) / (model.K_glu**model.n_glu + Glu_extra_control[-1]**model.n_glu)
    print(f"铁死亡速率 - 治疗组: {ferroptosis_rate_therapy:.6f} /hr, 对照组: {ferroptosis_rate_control:.6f} /hr")
    
    print("="*50)
