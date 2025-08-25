
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
转录翻译到谷氨酸生产模块 (改良版)

功能:
- 从T7活性到谷氨酸生产的完整动力学模拟
- 面向对象的模型设计
- 多条件比较和敏感性分析
- 高质量可视化

作者: iGEM建模团队
日期: 2025年8月
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import odeint
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib支持中文
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 加载参数文件
def load_parameters():
    """加载模型参数"""
    try:
        with open("params/promoters.json", "r", encoding='utf-8') as f: 
            promoters = json.load(f)
        with open("params/splitT7.json", "r", encoding='utf-8') as f: 
            splitp = json.load(f)
        print("已加载所有参数文件")
        return promoters, splitp
    except FileNotFoundError as e:
        print(f"警告: 参数文件缺失 {e}")
        print("使用默认参数，建议先运行前面的脚本生成参数")
        promoters = {
            "pPept": {"beta": 1000, "K": 1.0, "n": 2.0, "leaky": 100},
            "pLR": {"beta": 1000, "K": 37.0, "n": 2.0, "leaky": 100}
        }
        splitp = {"alpha": 1.0, "Kd": 1.0, "leaky": 0.0}
        return promoters, splitp

promoters, splitp = load_parameters()

def hill_act(x, beta, K, n, leaky): 
    """激活型Hill函数"""
    return leaky + beta * (x**n) / (K**n + x**n)

def hill_rep(x, beta, K, n, leaky): 
    """抑制型Hill函数"""
    return leaky + beta / (1.0 + (x / K)**n)

def pPept_output(O2):
    """pPept启动子对氧气的响应"""
    p = promoters["pPept"]
    return hill_rep(O2, p["beta"], p["K"], p["n"], p["leaky"])

def pLR_output(T):
    """pL/pR启动子对温度的响应"""
    p = promoters["pLR"]
    return hill_act(T, p["beta"], p["K"], p["n"], p["leaky"])

def split_t7_activity(A, B):
    """分裂T7活性计算"""
    return splitp["leaky"] + splitp["alpha"] * (A*B) / (splitp["Kd"] + A*B + 1e-12)

class TxTlGluModel:
    """转录翻译谷氨酸生产模型"""
    
    def __init__(self, **params):
        """
        初始化模型参数
        
        Parameters:
        -----------
        k_tx : float
            转录速率常数 (mRNA / (h · T7_unit))
        k_tl : float  
            翻译速率常数 (enzyme / (h · mRNA))
        dm : float
            mRNA降解速率 (1/h)
        dE : float
            酶降解/稀释速率 (1/h)
        kcat_gdh : float
            GDH酶催化常数 (1/h)
        Km_gdh : float
            GDH米氏常数 (mM)
        S_akg : float
            alpha-酮戊二酸浓度 (mM)
        k_sec : float
            胞外分泌速率 (1/h)
        """
        self.k_tx = params.get('k_tx', 1.0)
        self.k_tl = params.get('k_tl', 1.0) 
        self.dm = params.get('dm', 1.0)
        self.dE = params.get('dE', 0.1)
        self.kcat_gdh = params.get('kcat_gdh', 10.0)
        self.Km_gdh = params.get('Km_gdh', 0.1)
        self.S_akg = params.get('S_akg', 1.0)
        self.k_sec = params.get('k_sec', 0.05)
    
    def dydt(self, y, t, T7_activity):
        """ODE系统定义"""
        m, E, Cext = y
        
        # 转录速率
        dm_dt = self.k_tx * T7_activity - self.dm * m
        
        # 翻译速率
        dE_dt = self.k_tl * m - self.dE * E
        
        # 酶催化反应速率（Michaelis-Menten动力学）
        v_gdh = self.kcat_gdh * E * (self.S_akg / (self.Km_gdh + self.S_akg))
        
        # 胞外谷氨酸积累
        dC_dt = self.k_sec * v_gdh
        
        return [dm_dt, dE_dt, dC_dt]
    
    def simulate(self, O2_percent=1.0, Temp_C=40.0, t_end_h=24.0, dt_min=1.0):
        """
        模拟TX-TL到谷氨酸生产过程
        
        Parameters:
        -----------
        O2_percent : float
            氧气浓度百分比
        Temp_C : float
            温度（摄氏度）
        t_end_h : float
            模拟时长（小时）
        dt_min : float
            时间步长（分钟）
        
        Returns:
        --------
        tuple: (时间数组, mRNA, 酶浓度, 胞外谷氨酸)
        """
        # 计算T7活性
        A = pPept_output(O2_percent)
        B = pLR_output(Temp_C)
        T7 = split_t7_activity(A, B)
        
        print(f"模拟条件: O2={O2_percent}%, T={Temp_C}°C")
        print(f"pPept输出: {A:.1f}, pL/pR输出: {B:.1f}")
        print(f"T7活性: {T7:.1f}")
        
        # 时间数组
        t = np.linspace(0, t_end_h, int(t_end_h * 60 / dt_min))
        
        # 初始条件
        y0 = [0.0, 0.0, 0.0]  # [mRNA, enzyme, extracellular_glu]
        
        # 求解ODE
        solution = odeint(self.dydt, y0, t, args=(T7,))
        
        m = solution[:, 0]
        E = solution[:, 1] 
        Cext = solution[:, 2]
        
        return t, m, E, Cext
    
    def plot_results(self, t, m, E, Cext, save_fig=True):
        """绘制模拟结果"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
        
        # mRNA动力学
        ax1.plot(t, m, 'b-', linewidth=2)
        ax1.set_xlabel('时间 (h)')
        ax1.set_ylabel('mRNA (AU)')
        ax1.set_title('mRNA动力学')
        ax1.grid(True, alpha=0.3)
        
        # 酶浓度动力学
        ax2.plot(t, E, 'g-', linewidth=2)
        ax2.set_xlabel('时间 (h)')
        ax2.set_ylabel('酶浓度 (AU)')
        ax2.set_title('酶浓度动力学')
        ax2.grid(True, alpha=0.3)
        
        # 胞外谷氨酸积累
        ax3.plot(t, Cext, 'r-', linewidth=2)
        ax3.set_xlabel('时间 (h)')
        ax3.set_ylabel('胞外谷氨酸 (AU)')
        ax3.set_title('胞外谷氨酸积累')
        ax3.grid(True, alpha=0.3)
        
        # 生产速率
        production_rate = np.gradient(Cext, t)
        ax4.plot(t, production_rate, 'm-', linewidth=2)
        ax4.set_xlabel('时间 (h)')
        ax4.set_ylabel('生产速率 (AU/h)')
        ax4.set_title('谷氨酸生产速率')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_fig:
            plt.savefig('tx_tl_glu_dynamics.png', dpi=300, bbox_inches='tight')
            print("保存动力学图像: tx_tl_glu_dynamics.png")
        
        plt.show()
        
        return fig

def compare_conditions():
    """比较不同条件下的谷氨酸生产"""
    model = TxTlGluModel()
    
    conditions = [
        {"O2": 0.5, "T": 42.0, "label": "低氧+高温 (AND激活)"},
        {"O2": 21.0, "T": 42.0, "label": "常氧+高温 (部分激活)"},
        {"O2": 0.5, "T": 25.0, "label": "低氧+常温 (部分激活)"},
        {"O2": 21.0, "T": 25.0, "label": "常氧+常温 (基础水平)"}
    ]
    
    plt.figure(figsize=(10, 6))
    
    for cond in conditions:
        t, m, E, Cext = model.simulate(O2_percent=cond["O2"], Temp_C=cond["T"])
        plt.plot(t, Cext, linewidth=2, label=cond["label"])
    
    plt.xlabel('时间 (h)')
    plt.ylabel('胞外谷氨酸 (AU)')
    plt.title('不同条件下的谷氨酸生产比较')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('glu_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("保存比较图像: glu_comparison.png")

if __name__ == "__main__":
    print("=== TX-TL到谷氨酸生产模拟 ===")
    
    # 创建模型实例
    model = TxTlGluModel()
    
    # 默认模拟（AND门激活条件）
    print("\n默认模拟：低氧(0.5%) + 高温(42°C)")
    t, m, E, C = model.simulate(O2_percent=0.5, Temp_C=42.0)
    
    print(f"最终胞外谷氨酸浓度: {C[-1]:.2f} AU")
    print(f"最大酶浓度: {np.max(E):.2f} AU") 
    print(f"达到最大酶浓度时间: {t[np.argmax(E)]:.1f} h")
    
    # 绘制结果
    model.plot_results(t, m, E, C)
    
    # 比较不同条件
    print("\n=== 条件比较 ===")
    compare_conditions()
    
    # 参数敏感性分析
    print("\n=== 参数敏感性分析 ===")
    base_params = {'k_tx': 1.0, 'k_tl': 1.0, 'dm': 1.0, 'dE': 0.1, 
                   'kcat_gdh': 10.0, 'Km_gdh': 0.1, 'S_akg': 1.0, 'k_sec': 0.05}
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()
    
    sensitive_params = ['k_tx', 'k_tl', 'kcat_gdh', 'k_sec']
    factors = [0.5, 1.0, 2.0]
    
    for i, param in enumerate(sensitive_params):
        for factor in factors:
            params = base_params.copy()
            params[param] = base_params[param] * factor
            
            model_temp = TxTlGluModel(**params)
            t, _, _, C = model_temp.simulate(O2_percent=0.5, Temp_C=42.0)
            
            axes[i].plot(t, C, label=f'{param}×{factor}', linewidth=2)
        
        axes[i].set_xlabel('时间 (h)')
        axes[i].set_ylabel('胞外谷氨酸 (AU)')
        axes[i].set_title(f'{param} 敏感性')
        axes[i].legend()
        axes[i].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('parameter_sensitivity.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("保存敏感性分析图像: parameter_sensitivity.png")
