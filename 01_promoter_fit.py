#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
启动子传递函数拟合模块 (改良版)

功能:
- 使用Hill方程拟合启动子响应曲线
- 支持激活型和抑制型启动子
- 参数约束优化和拟合质量评估
- 高质量可视化和数据验证

数据说明:
- reporter_MEFL: Molecule    ax1.set_xlabel(    ax2.set_xlabel('温度 (°C)')
    ax2.set_ylabel('报告基因表达 (MEFL)')
    ax2.set_title(f'pL/pR启动子 (R2 = {pLR_pars["_r_squared"]:.3f})')浓度 (%)')
    ax1.set_ylabel('报告基因表达 (MEFL)')
    ax1.set_title(f'pPept启动子 (R2 = {pPept_pars["_r_squared"]:.3f})')f Equivalent Fluorescein per Liter
  这是标准化的荧光强度单位，表示启动子驱动的报告基因表达水平
  数值越高表示启动子活性越强

作者: iGEM建模团队
日期: 2025年8月
"""

import json
import os
import pandas as pd
import numpy as np
from lmfit import Model
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib支持中文
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# --- Hill/Logistic forms ---
def hill_act(x, beta, K, n, leaky):
    """
    激活型Hill函数 - 随着输入增加，输出增加
    
    Parameters:
    -----------
    x : array_like
        输入变量（如温度、浓度等）
    beta : float
        最大激活幅度 (MEFL单位)
    K : float
        半饱和常数 - 达到50%最大响应时的输入值
    n : float
        Hill系数 - 控制响应曲线的陡峭程度
    leaky : float
        基础表达水平 (MEFL单位)
    
    Returns:
    --------
    array_like : 启动子输出 (MEFL单位)
    """
    return leaky + beta * (x**n) / (K**n + x**n)

def hill_rep(x, beta, K, n, leaky):
    """
    抑制型Hill函数 - 随着输入增加，输出减少
    
    Parameters:
    -----------
    x : array_like
        输入变量（如氧气浓度等）
    beta : float
        最大抑制前的表达水平 (MEFL单位)
    K : float
        半抑制常数 - 50%抑制时的输入值
    n : float
        Hill系数 - 控制抑制曲线的陡峭程度
    leaky : float
        最小表达水平 (MEFL单位)
    
    Returns:
    --------
    array_like : 启动子输出 (MEFL单位)
    """
    return leaky + beta / (1.0 + (x / K)**n)

def fit_curve(df, xcol, ycol, mode="auto"):
    """
    拟合Hill函数到数据
    
    Parameters:
    -----------
    df : DataFrame
        包含数据的DataFrame
    xcol : str
        输入变量列名
    ycol : str
        输出变量列名
    mode : str
        拟合模式: "act"(激活), "rep"(抑制), "auto"(自动检测)
    
    Returns:
    --------
    tuple: (参数字典, 拟合结果对象)
    """
    x = df[xcol].values.astype(float)
    y = df[ycol].values.astype(float)
    
    # 数据验证
    if len(x) < 4:
        raise ValueError(f"数据点太少 ({len(x)} < 4), 无法进行可靠拟合")
    
    if mode == "act":
        f = Model(hill_act)
        params = f.make_params(
            beta=max(y)-min(y), 
            K=np.median(x)+1e-6, 
            n=2.0, 
            leaky=min(y)
        )
        # 设置参数约束
        params['beta'].min = 0
        params['K'].min = 1e-6
        params['n'].min = 0.1
        params['n'].max = 10
        params['leaky'].min = 0
        
    elif mode == "rep":
        f = Model(hill_rep)
        params = f.make_params(
            beta=max(y), 
            K=np.median(x)+1e-6, 
            n=2.0, 
            leaky=min(y)
        )
        # 设置参数约束
        params['beta'].min = 0
        params['K'].min = 1e-6
        params['n'].min = 0.1
        params['n'].max = 10
        params['leaky'].min = 0
        
    else:
        # 自动检测模式：基于相关性选择模型
        corr = np.corrcoef(x, y)[0,1]
        mode = "rep" if corr < 0 else "act"
        return fit_curve(df, xcol, ycol, mode)
    
    try:
        out = f.fit(y, x=x, params=params)
        pars = {k: float(v) for k, v in out.best_values.items()}
        pars["_mode"] = mode
        pars["_aic"] = float(out.aic)
        pars["_r_squared"] = 1 - out.residual.var() / np.var(y)
        return pars, out
    except Exception as e:
        print(f"拟合失败: {e}")
        return None, None

def plot_fits(pPept_df, pLR_df, pPept_mean, pLR_mean, pPept_fit, pLR_fit, pPept_pars, pLR_pars):
    """绘制拟合结果"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # pPept vs O2 plot
    ax1.scatter(pPept_df['O2_percent'], pPept_df['reporter_MEFL'], 
                alpha=0.6, c='lightblue', s=30, label='原始数据')
    ax1.plot(pPept_mean['O2_percent'], pPept_mean['reporter_MEFL'], 
             'bo-', markersize=6, label='平均值')
    
    x_fine = np.linspace(pPept_mean['O2_percent'].min(), 
                        pPept_mean['O2_percent'].max(), 100)
    y_pred = hill_rep(x_fine, pPept_pars['beta'], pPept_pars['K'], 
                     pPept_pars['n'], pPept_pars['leaky'])
    ax1.plot(x_fine, y_pred, 'r-', linewidth=2, label='Hill拟合')
    ax1.set_xlabel('氧气浓度 (%)')
    ax1.set_ylabel('报告基因 (MEFL)')
    ax1.set_title(f'pPept启动子 (R2 = {pPept_pars["_r_squared"]:.3f})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # pL/pR vs Temperature plot
    ax2.scatter(pLR_df['Temp_C'], pLR_df['reporter_MEFL'], 
                alpha=0.6, c='lightcoral', s=30, label='原始数据')
    ax2.plot(pLR_mean['Temp_C'], pLR_mean['reporter_MEFL'], 
             'ro-', markersize=6, label='平均值')
    
    x_fine = np.linspace(pLR_mean['Temp_C'].min(), 
                        pLR_mean['Temp_C'].max(), 100)
    y_pred = hill_act(x_fine, pLR_pars['beta'], pLR_pars['K'], 
                     pLR_pars['n'], pLR_pars['leaky'])
    ax2.plot(x_fine, y_pred, 'b-', linewidth=2, label='Hill拟合')
    ax2.set_xlabel('温度 (°C)')
    ax2.set_ylabel('报告基因 (MEFL)')
    ax2.set_title(f'pL/pR启动子 (R2 = {pLR_pars["_r_squared"]:.3f})')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('promoter_fits.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("保存拟合图像: promoter_fits.png")

# --- 主程序 ---
if __name__ == "__main__":
    print("=== 启动子传递函数拟合 ===")
    
    # 检查数据文件是否存在
    data_dir = Path("data")
    if not data_dir.exists():
        raise FileNotFoundError("数据目录 'data/' 不存在")
    
    pPept_file = data_dir / "pPept_O2_curve.csv"
    pLR_file = data_dir / "pLR_T_curve.csv"
    
    if not pPept_file.exists():
        raise FileNotFoundError(f"数据文件 {pPept_file} 不存在")
    if not pLR_file.exists():
        raise FileNotFoundError(f"数据文件 {pLR_file} 不存在")
    
    # 加载数据
    print("加载数据...")
    pPept_df = pd.read_csv(pPept_file)
    pLR_df = pd.read_csv(pLR_file)
    
    print(f"pPept数据: {len(pPept_df)} 行")
    print(f"pL/pR数据: {len(pLR_df)} 行")
    
    # 数据预处理：计算重复实验的平均值
    pPept_mean = pPept_df.groupby("O2_percent", as_index=False)["reporter_MEFL"].mean()
    pLR_mean = pLR_df.groupby("Temp_C", as_index=False)["reporter_MEFL"].mean()
    
    print("进行启动子响应曲线拟合...")
    
    # 拟合pPept启动子（氧气抑制型）
    print("拟合pPept启动子...")
    pPept_result = fit_curve(pPept_mean, "O2_percent", "reporter_MEFL", mode="rep")
    if pPept_result[0] is None:
        raise RuntimeError("pPept启动子拟合失败")
    pPept_pars, pPept_fit = pPept_result
    
    # 拟合pL/pR启动子（温度激活型）
    print("拟合pL/pR启动子...")
    pLR_result = fit_curve(pLR_mean, "Temp_C", "reporter_MEFL", mode="act")
    if pLR_result[0] is None:
        raise RuntimeError("pL/pR启动子拟合失败")
    pLR_pars, pLR_fit = pLR_result
    
    # 保存参数
    params = {"pPept": pPept_pars, "pLR": pLR_pars}
    params_dir = Path("params")
    params_dir.mkdir(exist_ok=True)
    
    params_file = params_dir / "promoters.json"
    with open(params_file, "w", encoding='utf-8') as f:
        json.dump(params, f, indent=2, ensure_ascii=False)
    
    print(f"参数已保存至: {params_file}")
    
    # 输出拟合结果摘要
    print("\n=== 拟合结果摘要 ===")
    print(f"pPept启动子 - R2: {pPept_pars['_r_squared']:.3f}, AIC: {pPept_pars['_aic']:.1f}")
    print(f"  beta={pPept_pars['beta']:.1f}, K={pPept_pars['K']:.2f}, n={pPept_pars['n']:.2f}")
    print(f"pL/pR启动子 - R2: {pLR_pars['_r_squared']:.3f}, AIC: {pLR_pars['_aic']:.1f}")
    print(f"  beta={pLR_pars['beta']:.1f}, K={pLR_pars['K']:.2f}, n={pLR_pars['n']:.2f}")
    
    # 绘制拟合图像
    plot_fits(pPept_df, pLR_df, pPept_mean, pLR_mean, 
             pPept_fit, pLR_fit, pPept_pars, pLR_pars)
