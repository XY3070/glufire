
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分裂T7 RNA聚合酶AND门建模模块

功能:
- 整合两个环境输入，建模分裂T7 RNA聚合酶活性
- 使用优化算法拟合参数
- 3D响应面可视化
- AND逻辑门验证

作者: iGEM建模团队
日期: 2025年8月
"""

import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib支持中文
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def split_t7_activity(A, B, alpha, Kd, leaky=0.0):
    """
    计算分裂T7 RNA聚合酶的有效活性
    
    Parameters:
    -----------
    A, B : float or array
        来自两个启动子的输出信号
    alpha : float
        最大活性参数
    Kd : float
        半饱和常数
    leaky : float
        基础泄漏活性
    
    Returns:
    --------
    float or array: T7 RNA聚合酶活性
    """
    product = A * B
    return leaky + alpha * product / (Kd + product + 1e-12)

# 加载启动子参数
try:
    with open("params/promoters.json", "r", encoding='utf-8') as f:
        promoters = json.load(f)
    print("已加载启动子参数")
except FileNotFoundError:
    print("警告: params/promoters.json 不存在，请先运行 01_promoter_fit.py")
    # 使用默认参数
    promoters = {
        "pPept": {"beta": 1000, "K": 1.0, "n": 2.0, "leaky": 100},
        "pLR": {"beta": 1000, "K": 37.0, "n": 2.0, "leaky": 100}
    }

# 从环境输入到启动子输出的映射函数
def hill_act(x, beta, K, n, leaky): 
    """激活型Hill函数"""
    return leaky + beta * (x**n) / (K**n + x**n)

def hill_rep(x, beta, K, n, leaky): 
    """抑制型Hill函数"""
    return leaky + beta / (1.0 + (x / K)**n)

def pPept_output(O2_percent):
    """pPept启动子对氧气浓度的响应"""
    p = promoters["pPept"]
    return hill_rep(O2_percent, p["beta"], p["K"], p["n"], p["leaky"])

def pLR_output(Temp_C):
    """pL/pR启动子对温度的响应"""
    p = promoters["pLR"]
    return hill_act(Temp_C, p["beta"], p["K"], p["n"], p["leaky"])

def optimize_split_T7_params(A, B, Y_observed):
    """
    优化分裂T7参数
    
    Parameters:
    -----------
    A, B : array
        输入信号A和B
    Y_observed : array
        观测到的T7活性
    
    Returns:
    --------
    dict: 优化后的参数
    """
    from scipy.optimize import minimize
    
    def objective(params):
        alpha, Kd, leaky = params
        if alpha <= 0 or Kd <= 0 or leaky < 0:
            return 1e6  # 惩罚非法参数
        Y_pred = split_t7_activity(A, B, alpha, Kd, leaky)
        mse = np.mean((Y_pred - Y_observed)**2)
        return mse
    
    # 初始猜测
    initial_guess = [np.max(Y_observed), np.mean(A*B), np.min(Y_observed)]
    
    # 参数边界
    bounds = [(1e-6, 10*np.max(Y_observed)),  # alpha
              (1e-6, 10*np.max(A*B)),         # Kd  
              (0, np.max(Y_observed))]         # leaky
    
    result = minimize(objective, initial_guess, bounds=bounds, method='L-BFGS-B')
    
    if result.success:
        alpha, Kd, leaky = result.x
        Y_pred = split_t7_activity(A, B, alpha, Kd, leaky)
        r2 = r2_score(Y_observed, Y_pred)
        
        return {
            "alpha": float(alpha),
            "Kd": float(Kd), 
            "leaky": float(leaky),
            "mse": float(result.fun),
            "r_squared": float(r2)
        }
    else:
        print("优化失败，使用默认参数")
        return {"alpha": 1.0, "Kd": 1.0, "leaky": 0.0, "mse": np.inf, "r_squared": 0.0}

def plot_split_T7_response(scan_data, params):
    """绘制分裂T7响应曲线"""
    if scan_data is None:
        return
        
    A = scan_data["A_u"].values
    B = scan_data["B_u"].values  
    Y_obs = scan_data["T7_reporter_MEFL"].values
    Y_pred = split_t7_activity(A, B, params["alpha"], params["Kd"], params["leaky"])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 预测vs观测
    ax1.scatter(Y_obs, Y_pred, alpha=0.6, s=50)
    min_val, max_val = min(np.min(Y_obs), np.min(Y_pred)), max(np.max(Y_obs), np.max(Y_pred))
    ax1.plot([min_val, max_val], [min_val, max_val], 'r--', label='理想拟合')
    ax1.set_xlabel('观测T7活性')
    ax1.set_ylabel('预测T7活性')
    ax1.set_title(f'预测vs观测 (R2 = {params["r_squared"]:.3f})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 3D响应面（投影）
    product = A * B
    sort_idx = np.argsort(product)
    ax2.plot(product[sort_idx], Y_obs[sort_idx], 'bo', markersize=4, label='观测数据')
    ax2.plot(product[sort_idx], Y_pred[sort_idx], 'r-', linewidth=2, label='模型预测')
    ax2.set_xlabel('A × B (输入信号乘积)')
    ax2.set_ylabel('T7活性')
    ax2.set_title('分裂T7响应曲线')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('splitT7_fit.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("保存分裂T7拟合图像: splitT7_fit.png")

# 主程序
if __name__ == "__main__":
    print("=== 分裂T7 AND门模型 ===")
    
    # 尝试加载分裂T7扫描数据
    scan_data = None
    try:
        scan_data = pd.read_csv("data/splitT7_scan.csv")
        print(f"加载分裂T7扫描数据: {len(scan_data)} 个数据点")
        
        A = scan_data["A_u"].values
        B = scan_data["B_u"].values
        Y = scan_data["T7_reporter_MEFL"].values
        
        print("优化分裂T7参数...")
        params = optimize_split_T7_params(A, B, Y)
        
        print(f"拟合结果 - R2: {params['r_squared']:.3f}, MSE: {params['mse']:.1f}")
        print(f"参数: alpha={params['alpha']:.3f}, Kd={params['Kd']:.3f}, leaky={params['leaky']:.3f}")
        
    except FileNotFoundError:
        print("未找到splitT7_scan.csv，使用默认参数")
        params = {"alpha": 1.0, "Kd": 1.0, "leaky": 0.0}
    
    # 保存参数
    params_dir = Path("params")
    params_dir.mkdir(exist_ok=True)
    
    params_file = params_dir / "splitT7.json"
    with open(params_file, "w", encoding='utf-8') as f:
        json.dump(params, f, indent=2, ensure_ascii=False)
    
    print(f"参数已保存至: {params_file}")
    
    # 绘制拟合结果
    if scan_data is not None:
        plot_split_T7_response(scan_data, params)
    
    # 演示AND门逻辑
    print("\n=== AND门逻辑演示 ===")
    print("条件: 低氧(0.5%) + 高温(42°C)")
    A_demo = pPept_output(0.5)  # 低氧条件
    B_demo = pLR_output(42.0)   # 高温条件
    T7_demo = split_t7_activity(A_demo, B_demo, params["alpha"], params["Kd"], params["leaky"])
    print(f"pPept输出: {A_demo:.1f}")
    print(f"pL/pR输出: {B_demo:.1f}")
    print(f"T7活性: {T7_demo:.1f}")
    
    print("\n条件: 常氧(21%) + 高温(42°C)")
    A_demo2 = pPept_output(21.0)  # 常氧条件
    T7_demo2 = split_t7_activity(A_demo2, B_demo, params["alpha"], params["Kd"], params["leaky"])
    print(f"pPept输出: {A_demo2:.1f}")
    print(f"T7活性: {T7_demo2:.1f} (应显著降低)")
