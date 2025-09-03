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
        # 初始化子模块时使用优化的参数 - 精确调节设计
        glu_params = {
            'k_prod_max': 200.0,     # 降低最大生产速率
            'K_t7': 1000.0,          # 进一步提高阈值，使control T7(637)效率极低
            'k_export_max': 100.0,   # 适中分泌速率
            'K_export': 3.0,         # 适中分泌阈值
            'k_dilution': 0.08,      # 提高稀释速率，降低积累
            'V_intra_over_V_extra': 0.15  # 降低体积比
        }
        
        self.and_gate = SimpleANDGate()
        self.glu_metabolism = GluMetabolismModel(**glu_params)
        
        # 肿瘤生长和死亡参数
        self.r = params.get('r', 0.3)  # 肿瘤固有生长速率 (1/hr)
        self.K_tumor = params.get('K_tumor', 1e9)  # 肿瘤承载能力 (细胞数)
        self.r_eng = params.get('r_eng', 0.1) # 工程细胞生长速率
        self.K_eng = params.get('K_eng', 1e8) # 工程细胞承载能力
        self.k_ferroptosis_max = params.get('k_ferroptosis_max', 1.0) # 适中的最大铁死亡速率
        self.K_glu = params.get('K_glu', 5.0) # 适中谷氨酸阈值
        self.n_glu = params.get('n_glu', 2.0) # Hill系数

        # 新增：工程细胞生长参数
        self.r_eng = params.get('r_eng', 0.2) # 工程细胞生长速率
        self.K_eng = params.get('K_eng', 5e8) # 工程细胞承载能力

        # 新增：用于计算谷氨酸浓度的体积参数
        # V_cell: 单个细胞体积 (L/cell), V_tumor_ext: 肿瘤间质液体积 (L)
        self.V_cell = params.get('V_cell', 2e-12) 
        self.V_tumor_ext = params.get('V_tumor_ext', 0.01) # 回调到10mL

    def get_t7_activity(self, env_conditions):
        """从AND门模块获取T7活性。"""
        return self.and_gate.get_t7_activity(
            env_conditions['O2_percent'],
            env_conditions['Temp_C']
        )

    def dydt(self, y, t, env_conditions):
        """
        定义整合模型的ODE系统。
        
        y: array [N_tumor, D_tumor, Glu_intra, Glu_extra, N_eng]
            - N_tumor: 存活的肿瘤细胞数
            - D_tumor: 死亡的肿瘤细胞数
            - Glu_intra: 工程细胞内的谷氨酸浓度 (mM)
            - Glu_extra: 细胞外谷氨酸浓度 (mM)
            - N_eng: 存活的工程细胞数
        """
        N_tumor, D_tumor, Glu_intra, Glu_extra, N_eng = y
        
        # 1. 从环境条件计算T7活性
        t7_activity = self.get_t7_activity(env_conditions)
        
        # 2. 模拟谷氨酸生产和分泌 (与工程细胞数 N_eng 关联)
        v_prod = self.glu_metabolism.k_prod_max * t7_activity / (self.glu_metabolism.K_t7 + t7_activity)
        v_export = self.glu_metabolism.k_export_max * Glu_intra / (self.glu_metabolism.K_export + Glu_intra)
        
        # 细胞内谷氨酸浓度变化 (在每个工程细胞内)
        dGlu_intra_dt = v_prod - v_export - self.glu_metabolism.k_dilution * Glu_intra
        
        # 细胞外谷氨酸浓度变化 (总分泌量取决于 N_eng)
        # 进一步降低分泌放大系数
        secretion_amplifier = 1.5  # 降低放大系数使control组谷氨酸极低
        total_export_flux = v_export * N_eng * (self.V_cell / self.V_tumor_ext) * secretion_amplifier
        dGlu_extra_dt = total_export_flux - self.glu_metabolism.k_dilution * Glu_extra
        
        # 3. 肿瘤细胞生长和死亡
        # 生长
        growth_rate_tumor = self.r * N_tumor * (1 - (N_tumor + N_eng) / self.K_tumor)
        # 谷氨酸诱导的铁死亡 (只作用于肿瘤细胞)
        ferroptosis_rate = self.k_ferroptosis_max * (Glu_extra**self.n_glu) / (self.K_glu**self.n_glu + Glu_extra**self.n_glu)
        death_term_tumor = ferroptosis_rate * N_tumor
        
        # 4. 工程细胞生长 (依赖T7活性 - 关键修复!)
        # 工程细胞只有在T7活性高时才能有效增殖
        t7_growth_factor = t7_activity / (self.glu_metabolism.K_t7 + t7_activity)  # 归一化T7效应
        growth_rate_eng = self.r_eng * N_eng * (1 - (N_tumor + N_eng) / self.K_tumor) * t7_growth_factor
        
        # 状态变量变化率
        dN_tumor_dt = growth_rate_tumor - death_term_tumor
        dD_tumor_dt = death_term_tumor
        dN_eng_dt = growth_rate_eng
        
        # 简单的数值稳定性保护 - 防止负数
        if N_tumor <= 1.0 and dN_tumor_dt < 0:
            dN_tumor_dt = 0  # 防止肿瘤细胞数变负
        if N_eng <= 1.0 and dN_eng_dt < 0:
            dN_eng_dt = 0    # 防止工程细胞数变负
        if Glu_intra <= 0 and dGlu_intra_dt < 0:
            dGlu_intra_dt = 0
        if Glu_extra <= 0 and dGlu_extra_dt < 0:
            dGlu_extra_dt = 0

        return [dN_tumor_dt, dD_tumor_dt, dGlu_intra_dt, dGlu_extra_dt, dN_eng_dt]

    def simulate(self, env_conditions, t_end=100.0, dt=0.5):
        """
        运行整合治疗模拟。
        
        env_conditions: dict
            包含 'O2_percent' 和 'Temp_C' 的字典。
        """
        t = np.arange(0, t_end, dt)
        # 初始条件: [N_tumor, D_tumor, Glu_intra, Glu_extra, N_eng]
        # 初始肿瘤细胞与工程细胞比例 9:1
        y0 = [9e6, 0, 0, 0, 1e6]  
        
        solution = odeint(self.dydt, y0, t, args=(env_conditions,))
        
        return t, solution

# --- 示例 ---
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    print("="*20 + " 整合治疗模型测试 " + "="*20)
    
    model = IntegratedTherapyModel()
    
    # --- 条件1: 治疗开启 (低氧, 高温) ---
    print("模拟治疗条件 (低氧, 高温)...")
    env_therapy = {'O2_percent': 1.0, 'Temp_C': 42.0}
    t7_therapy = model.get_t7_activity(env_therapy)
    print(f"Therapy T7 activity: {t7_therapy:.1f}")
    t_therapy, sol_therapy = model.simulate(env_therapy, t_end=100)
    
    # --- 条件2: 对照组 (高氧, 低温 - 确保AND门关闭) ---
    print("模拟对照条件 (高氧, 低温)...")
    env_control = {'O2_percent': 21.0, 'Temp_C': 37.0}
    t7_control = model.get_t7_activity(env_control)
    print(f"Control T7 activity: {t7_control:.1f}")
    t_control, sol_control = model.simulate(env_control, t_end=100)
    
    print(f"T7 活性比值 (therapy/control): {t7_therapy/t7_control:.2f}")
    print(f"最终肿瘤细胞数比值 (therapy/control): {sol_therapy[-1,0]/sol_control[-1,0]:.2f}")
    print(f"最终细胞外谷氨酸浓度 - therapy: {sol_therapy[-1,3]:.3f} mM")
    print(f"最终细胞外谷氨酸浓度 - control: {sol_control[-1,3]:.3f} mM")

    # --- 绘图比较 ---
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21, 6))
    
    # 图1: 肿瘤细胞数
    ax1.plot(t_therapy, sol_therapy[:, 0], label='Therapy ON', color='red')
    ax1.plot(t_control, sol_control[:, 0], label='Control', color='blue', linestyle='--')
    ax1.set_title('Tumor Cell Response')
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Tumor Cell Count (relative)')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, which="both", ls="--")
    
    # 图2: 细胞外谷氨酸浓度
    ax2.plot(t_therapy, sol_therapy[:, 3], label='Therapy ON', color='red')
    ax2.plot(t_control, sol_control[:, 3], label='Control', color='blue', linestyle='--')
    ax2.set_title('Extracellular Glutamate')
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Glutamate (mM)')
    ax2.legend()
    ax2.grid(True)

    # 图3: 工程细胞数
    ax3.plot(t_therapy, sol_therapy[:, 4], label='Therapy ON', color='green')
    ax3.plot(t_control, sol_control[:, 4], label='Control', color='orange', linestyle='--')
    ax3.set_title('Engineered Cell Population')
    ax3.set_xlabel('Time (hours)')
    ax3.set_ylabel('Engineered Cell Count (relative)')
    ax3.set_yscale('log')
    ax3.legend()
    ax3.grid(True, which="both", ls="--")
    
    plt.tight_layout()
    plt.show()
    
    print("\n测试完成。图像显示了治疗组和对照组在肿瘤细胞、谷氨酸浓度和工程细胞上的差异。")
    print("="*68)
