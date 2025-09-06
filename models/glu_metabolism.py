#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
谷氨酸代谢模�?(Glutamate Metabolism Module)

功能:
- 接收上游T7聚合酶的活性作为输入�?
- 模拟谷氨酸在细胞内的生产和向细胞外的分泌过程�?
- 模型基于常微分方�?ODE)，描述了细胞内和细胞外谷氨酸浓度的动态变化�?

如何使用:
- 实例�?`GluMetabolismModel`�?
- 调用 `simulate` 方法，并提供T7聚合酶的活性，以运行模拟�?
"""

import numpy as np
from scipy.integrate import odeint

class GluMetabolismModel:
    """
    谷氨酸生产和分泌的ODE模型，包含信号转导与关键酶表达动力学�?
    
    该模型将T7聚合酶的活性与谷氨酸合成速率关联起来，并模拟
    谷氨酸通过分泌系统(如一个简单的转运蛋白)从细胞内运输到细胞外的过程�?
    """
    def __init__(self, **params):
        """
        初始化模型参数�?
        
        Parameters:
        -----------
        k_prod_max : float, optional
            最大谷氨酸生产速率 (mM/hr), 默认�?50.0�?
        K_t7 : float, optional
            T7活性达到半最大生产速率时的�?(AU), 默认�?500.0�?
        k_export_max : float, optional
            最大谷氨酸分泌速率 (mM/hr), 默认�?100.0�?
        K_export : float, optional
            细胞内谷氨酸浓度达到半最大分泌速率时的�?(mM), 默认�?10.0�?
        k_dilution : float, optional
            细胞生长或降解导致的稀�?降解速率 (1/hr), 默认�?0.1�?
        V_intra_over_V_extra : float, optional
            细胞内总体积与细胞外总体积的比率, 默认�?0.01�?
        k_syn_icd : float, optional
            Icd合成速率 (1/hr), 默认�?2.0�?
        k_syn_gdhA : float, optional
            gdhA合成速率 (1/hr), 默认�?2.0�?
        k_deg_icd : float, optional
            Icd降解速率 (1/hr), 默认�?0.2�?
        k_deg_gdhA : float, optional
            gdhA降解速率 (1/hr), 默认�?0.2�?
        Vmax_icd : float, optional
            Icd最大催化速率 (mM/hr), 默认�?100.0�?
        K_icd : float, optional
            Icd底物常数 (mM), 默认�?5.0�?
        Vmax_gdhA : float, optional
            gdhA最大催化速率 (mM/hr), 默认�?100.0�?
        K_gdhA : float, optional
            gdhA底物常数 (mM), 默认�?5.0�?
        n_hill : float, optional
            Hill系数，用于增强开关效�? 默认�?4.0�?
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
        self.n_hill = params.get('n_hill', 4.0)            # Hill系数，用于增强开关效�?

    def dydt(self, y, t, t7_activity):
        """
        定义谷氨酸代谢的常微分方程组�?

        y: array
            状态变�?[Glu_intra, Glu_extra, Icd, gdhA]
            - Glu_intra: 细胞内谷氨酸浓度 (mM)
            - Glu_extra: 细胞外谷氨酸浓度 (mM)
            - Icd: Icd表达水平
            - gdhA: gdhA表达水平
        t: float
            时间
        t7_activity: float
            T7聚合酶的活�?(AU)
        """
        Glu_intra, Glu_extra, Icd, gdhA = y
        
        # 信号转导: T7驱动Icd/gdhA表达 (使用Hill函数实现开关效�?
        t7_signal = (t7_activity**self.n_hill) / (self.K_t7**self.n_hill + t7_activity**self.n_hill)
        dIcd_dt = self.k_syn_icd * t7_signal - self.k_deg_icd * Icd
        dgdhA_dt = self.k_syn_gdhA * t7_signal - self.k_deg_gdhA * gdhA
        
        # 代谢通路: 酶促反应速率简化为与酶浓度成正�?
        v_prod = self.Vmax_gdhA * gdhA
        
        # 2. 谷氨酸分泌速率 (从细胞内到细胞外，Michaelis-Menten形式)
        v_export = self.k_export_max * Glu_intra / (self.K_export + Glu_intra)
        
        # 谷氨酸浓度变�?
        dGlu_intra_dt = v_prod - v_export - self.k_dilution * Glu_intra
        dGlu_extra_dt = v_export * self.V_ratio - self.k_dilution * Glu_extra
        
        return [dGlu_intra_dt, dGlu_extra_dt, dIcd_dt, dgdhA_dt]

    def simulate(self, t7_activity, t_end=24.0, dt=0.1):
        """
        运行ODE模拟�?

        Parameters:
        -----------
        t7_activity : float
            恒定的T7聚合酶活�?(AU)�?
        t_end : float, optional
            模拟结束时间 (小时), 默认�?24.0�?
        dt : float, optional
            模拟时间步长 (小时), 默认�?0.1�?

        Returns:
        --------
        tuple: (t, solution)
            - t: 时间点数�?
            - solution: 状态变�?[Glu_intra, Glu_extra] 在每个时间点的解
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

    print("="*20 + " 谷氨酸代谢模块优化测�?" + "="*20)
    
    # --- 1. 定义优化后的参数 ---
    optimized_params = {
        'K_t7': 800.0,           # 降低激活阈值，使高T7(3000)能有效激�?
        'k_syn_icd': 5.0,        # 提高合成速率，确保酶能积�?
        'k_syn_gdhA': 5.0,       # 提高合成速率，确保酶能积�?
        'k_deg_icd': 0.3,        # 适度降解速率
        'k_deg_gdhA': 0.3,       # 适度降解速率
        'k_dilution': 0.15,      # 适度稀释速率
        'k_export_max': 100.0,
        'K_export': 5.0,
        'Vmax_gdhA': 50.0,       # 调整最大生产速率，避免过度生�?
        'n_hill': 3.0            # 适度的Hill系数
    }
    
    model = GluMetabolismModel(**optimized_params)
    print("模型已使用优化参数创建�?)

    # --- 2. 定义高低T7活�?---
    t7_activity_high = 3000  # 确保高于K_t7
    t7_activity_low = 50     # 确保远低于K_t7
    
    print(f"高T7活�? {t7_activity_high} AU, 低T7活�? {t7_activity_low} AU")

    # --- 3. 运行模拟 ---
    t_high, sol_high = model.simulate(t7_activity=t7_activity_high, t_end=48)
    t_low, sol_low = model.simulate(t7_activity=t7_activity_low, t_end=48)
    print("模拟完成�?)

    # --- 4. 提取并分析结�?---
    final_high = sol_high[:, 1][-1]
    final_low = sol_low[:, 1][-1]
    ratio = final_high / final_low if final_low > 1e-9 else float('inf')

    print("\n--- 结果分析 ---")
    print(f"最终细胞外谷氨酸浓�?(高T7): {final_high:.4f} mM")
    print(f"最终细胞外谷氨酸浓�?(低T7): {final_low:.4f} mM")
    print(f"比�?(�?�?: {ratio:.2f}")
    
    # 检查酶表达水平
    icd_high_final = sol_high[-1, 2]
    gdhA_high_final = sol_high[-1, 3]
    icd_low_final = sol_low[-1, 2]
    gdhA_low_final = sol_low[-1, 3]
    
    print(f"\n--- 酶表达水平检�?---")
    print(f"最终Icd水平 (高T7): {icd_high_final:.6f}")
    print(f"最终Icd水平 (低T7): {icd_low_final:.6f}")
    print(f"最终gdhA水平 (高T7): {gdhA_high_final:.6f}")
    print(f"最终gdhA水平 (低T7): {gdhA_low_final:.6f}")
    
    # 验证Hill函数
    t7_signal_high = (t7_activity_high**model.n_hill) / (model.K_t7**model.n_hill + t7_activity_high**model.n_hill)
    t7_signal_low = (t7_activity_low**model.n_hill) / (model.K_t7**model.n_hill + t7_activity_low**model.n_hill)
    print(f"Hill信号 (高T7): {t7_signal_high:.6f}")
    print(f"Hill信号 (低T7): {t7_signal_low:.6f}")

    # --- 5. 计算分析指标 ---
    # 提取详细数据
    glu_intra_high, glu_extra_high, icd_high, gdhA_high = sol_high.T
    glu_intra_low, glu_extra_low, icd_low, gdhA_low = sol_low.T
    
    # 计算生产速率和分泌速率
    v_prod_high = model.Vmax_gdhA * gdhA_high
    v_prod_low = model.Vmax_gdhA * gdhA_low
    
    v_export_high = model.k_export_max * glu_intra_high / (model.K_export + glu_intra_high)
    v_export_low = model.k_export_max * glu_intra_low / (model.K_export + glu_intra_low)
    
    # 细胞外谷氨酸对比（归一化）
    glu_ratio = glu_extra_high / (glu_extra_low + 1e-9)
    
    # --- 6. 生成综合分析�?---
    fig = plt.figure(figsize=(18, 12))
    
    # 子图1: 谷氨酸浓度动�?- 分别显示高低T7
    ax1 = plt.subplot(2, 3, 1)
    # 只显示高T7的数据，因为低T7几乎�?
    plt.plot(t_high, glu_intra_high, label='细胞�?, color='blue', linewidth=2.5)
    plt.plot(t_high, glu_extra_high, label='细胞�?, color='red', linewidth=2.5)
    plt.title('谷氨酸浓度动�?(高T7)', fontsize=14)
    plt.xlabel('时间 (小时)')
    plt.ylabel('浓度 (mM)')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    # 添加文本说明低T7情况
    plt.text(0.6*t_high[-1], 0.8*max(glu_extra_high), 
             f'低T7条件�?\n细胞�? {final_low:.4f} mM\n细胞�? {glu_intra_low[-1]:.4f} mM', 
             fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7))
    
    # 子图2: 酶活性动�?- 只显示高T7，用插图显示低T7
    ax2 = plt.subplot(2, 3, 2)
    # 主图：高T7条件
    line1 = plt.plot(t_high, icd_high, label='Icd', color='green', linewidth=2.5)
    line2 = plt.plot(t_high, gdhA_high, label='gdhA', color='magenta', linewidth=2.5)
    plt.title('酶活性动�?(高T7)', fontsize=14)
    plt.xlabel('时间 (小时)')
    plt.ylabel('表达水平')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    
    # 在主图中添加小的插图显示低T7
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    ax2_inset = inset_axes(ax2, width="35%", height="35%", loc='upper right')
    ax2_inset.plot(t_low, icd_low, color='darkgreen', linewidth=1.5, label='Icd(低T7)')
    ax2_inset.plot(t_low, gdhA_low, color='purple', linewidth=1.5, label='gdhA(低T7)')
    ax2_inset.set_title('低T7', fontsize=9)
    ax2_inset.tick_params(labelsize=8)
    ax2_inset.grid(True, alpha=0.3)
    
    # 子图3: 细胞外谷氨酸对比
    ax3 = plt.subplot(2, 3, 3)
    plt.plot(t_high, glu_ratio, color='red', linewidth=3)
    plt.title('细胞外谷氨酸对比', fontsize=14)
    plt.xlabel('时间 (小时)')
    plt.ylabel('高T7/低T7比�?)
    plt.grid(True, alpha=0.3)
    plt.text(0.7*t_high[-1], 0.8*max(glu_ratio), f'最�? {glu_ratio[-1]:.0f}�?, 
             fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # 子图4: 谷氨酸生产速率
    ax4 = plt.subplot(2, 3, 4)
    plt.plot(t_high, v_prod_high, label='高T7', color='orange', linewidth=2)
    plt.plot(t_low, v_prod_low, label='低T7', color='brown', linestyle='--', linewidth=2)
    plt.title('谷氨酸生产速率', fontsize=14)
    plt.xlabel('时间 (小时)')
    plt.ylabel('生产速率 (mM/hr)')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # 子图5: 谷氨酸分泌速率
    ax5 = plt.subplot(2, 3, 5)
    plt.plot(t_high, v_export_high, label='高T7', color='cyan', linewidth=2)
    plt.plot(t_low, v_export_low, label='低T7', color='teal', linestyle='--', linewidth=2)
    plt.title('谷氨酸分泌速率', fontsize=14)
    plt.xlabel('时间 (小时)')
    plt.ylabel('分泌速率 (mM/hr)')
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # 子图6: 最终状态对�?- 分开显示不同数量级的变量
    ax6 = plt.subplot(2, 3, 6)
    
    # 创建双Y轴来显示不同数量级的数据
    ax6_enzyme = ax6.twinx()
    
    # 谷氨酸浓度数�?(左Y�?
    glu_categories = ['细胞外Glu', '细胞内Glu']
    glu_high_values = [final_high, glu_intra_high[-1]]
    glu_low_values = [final_low, glu_intra_low[-1]]
    
    x_glu = np.arange(len(glu_categories))
    width = 0.35
    
    bars1 = ax6.bar(x_glu - width/2, glu_high_values, width, label='高T7', color='red', alpha=0.8)
    bars2 = ax6.bar(x_glu + width/2, glu_low_values, width, label='低T7', color='blue', alpha=0.8)
    
    # 酶水平数�?(右Y�?
    enzyme_categories = ['Icd水平', 'gdhA水平']
    enzyme_high_values = [icd_high[-1], gdhA_high[-1]]
    enzyme_low_values = [icd_low[-1], gdhA_low[-1]]
    
    x_enzyme = np.arange(len(enzyme_categories)) + 3  # 偏移位置避免重叠
    
    bars3 = ax6_enzyme.bar(x_enzyme - width/2, enzyme_high_values, width, label='高T7(�?', color='orange', alpha=0.8)
    bars4 = ax6_enzyme.bar(x_enzyme + width/2, enzyme_low_values, width, label='低T7(�?', color='green', alpha=0.8)
    
    # 设置标签和标�?
    all_categories = glu_categories + enzyme_categories
    all_x = list(x_glu) + list(x_enzyme)
    
    ax6.set_title('最终状态对�?, fontsize=14)
    ax6.set_xlabel('变量')
    ax6.set_ylabel('谷氨酸浓�?(mM)', color='red')
    ax6_enzyme.set_ylabel('酶表达水�?, color='orange')
    
    # 设置X轴标�?
    ax6.set_xticks(all_x)
    ax6.set_xticklabels(all_categories, rotation=45, ha='right')
    
    # 设置Y轴颜�?
    ax6.tick_params(axis='y', labelcolor='red')
    ax6_enzyme.tick_params(axis='y', labelcolor='orange')
    
    # 合并图例
    lines1, labels1 = ax6.get_legend_handles_labels()
    lines2, labels2 = ax6_enzyme.get_legend_handles_labels()
    ax6.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)
    
    ax6.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # 保存图像
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, '..', 'results', 'glutamate_comprehensive_analysis.png')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n综合分析图已保存�? {os.path.abspath(output_path)}")
    
    try:
        plt.show()
    except Exception as e:
        print(f"无法自动显示图像: {e}")

    print("="*50)
