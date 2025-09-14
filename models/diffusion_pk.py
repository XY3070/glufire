#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
扩散与药代动力学模块 (Diffusion and Pharmacokinetics Module)
"""

import numpy as np
from scipy.integrate import odeint
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
from scipy.signal import savgol_filter  # ← 用于平滑导数
import sys
import os

# 将项目根目录添加到 Python 模块搜索路径
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from config_manager import ConfigManager


class MultiCompartmentPK:
    """多室药代动力学(PK)模型。"""
    def __init__(self, config_manager: ConfigManager = None, **params):
        if config_manager:
            pk_params = config_manager.get_params('diffusion_pk').get('MultiCompartmentPK', {})
        else:
            pk_params = {}
        pk_params.update(params)

        # 体积 (L)
        self.V_blood = pk_params.get('V_blood', 5.0)
        self.V_liver = pk_params.get('V_liver', 1.5)
        self.V_tumor = pk_params.get('V_tumor', 0.5)
        self.V_other = pk_params.get('V_other', 60.0)
        # 流速 (L/hr)
        self.q_blood_liver = pk_params.get('q_blood_liver', 90.0)
        self.q_blood_tumor = pk_params.get('q_blood_tumor', 10.0)
        self.q_blood_other = pk_params.get('q_blood_other', 200.0)
        # 清除率 (1/hr)
        self.k_elim_liver = pk_params.get('k_elim_liver', 0.5)
        self.k_uptake_tumor = pk_params.get('k_uptake_tumor', 0.2)

    def dydt(self, y, t, infusion_rate):
        C_blood, C_liver, C_tumor, C_other = y
        dC_blood_dt = (infusion_rate / self.V_blood +
                       (self.q_blood_liver * C_liver - self.q_blood_liver * C_blood) / self.V_blood +
                       (self.q_blood_tumor * C_tumor - self.q_blood_tumor * C_blood) / self.V_blood +
                       (self.q_blood_other * C_other - self.q_blood_other * C_blood) / self.V_blood)
        dC_liver_dt = ((self.q_blood_liver * C_blood - self.q_blood_liver * C_liver) / self.V_liver -
                       self.k_elim_liver * C_liver)
        dC_tumor_dt = ((self.q_blood_tumor * C_blood - self.q_blood_tumor * C_tumor) / self.V_tumor -
                       self.k_uptake_tumor * C_tumor)
        dC_other_dt = ((self.q_blood_other * C_blood - self.q_blood_other * C_other) / self.V_other)
        return [dC_blood_dt, dC_liver_dt, dC_tumor_dt, dC_other_dt]

    def simulate(self, infusion_rate_func, t_end=24.0, dt=0.1):
        t = np.arange(0, t_end, dt)
        y0 = [0, 0, 0, 0]
        infusion_rates = [infusion_rate_func(ti) for ti in t]
        solution = np.zeros((len(t), len(y0)))
        solution[0] = y0
        for i in range(1, len(t)):
            y_prev = solution[i-1]
            dydt_val = self.dydt(y_prev, t[i-1], infusion_rates[i-1])
            solution[i] = y_prev + np.array(dydt_val) * dt
        return t, solution


class TumorDiffusion:
    """肿瘤微环境中的反应-扩散模型 (1D)。"""
    def __init__(self, config_manager: ConfigManager = None, **params):
        if config_manager:
            diffusion_params = config_manager.get_params('diffusion_pk').get('TumorDiffusion', {})
        else:
            diffusion_params = {}
        diffusion_params.update(params)
        self.D = diffusion_params.get('D', 0.01)      # cm²/hr
        self.k_uptake = diffusion_params.get('k_uptake', 0.2)  # 1/hr
        self.L = diffusion_params.get('L', 1.0)       # cm
        self.Nx = diffusion_params.get('Nx', 50)
        self.dx = self.L / (self.Nx - 1)

    def simulate(self, C_boundary, t_end=10.0, dt=0.01):
        t = np.arange(0, t_end, dt)
        Nt = len(t)
        C = np.zeros((self.Nx, Nt))
        C[0, :] = C_boundary
        C[-1, :] = C_boundary

        A_diag = np.ones(self.Nx) * (1 + 2 * self.D * dt / self.dx**2 + self.k_uptake * dt)
        A_upper = np.zeros(self.Nx); A_upper[:-1] = -self.D * dt / self.dx**2
        A_lower = np.zeros(self.Nx); A_lower[1:]  = -self.D * dt / self.dx**2
        A = spdiags([A_lower, A_diag, A_upper], [-1, 0, 1], self.Nx, self.Nx, format='csc')

        for n in range(Nt - 1):
            b = C[:, n]
            b[0]  += self.D * dt / self.dx**2 * C_boundary
            b[-1] += self.D * dt / self.dx**2 * C_boundary
            C[:, n+1] = spsolve(A, b)
        return self.dx * np.arange(self.Nx), t, C


# === Adapter: Glu_extra(t) -> S_t(t) (μmol/h) -> 三室PK -> 风险评估 ===

def _finite_diff_flux_from_conc(
    t_h: np.ndarray,
    glu_extra_mM: np.ndarray,
    V_tumor_ext_L: float,
    f_leak: float = 0.05,   # 只有一小部分变化泄漏入系统；可调 0.01~0.1
) -> np.ndarray:
    """
    用平滑导数把 C_extra(t) 映射为“泄漏到系统循环”的等效通量 S_t(t)：
      1) Savitzky–Golay 平滑一阶导，抑制数值噪声尖峰；
      2) **允许负值**（回流/净清除），避免只进不出；
      3) 乘以漏系数 f_leak，表达局部→系统的弱耦合。
    单位换算：mM*L = mmol；×1000 → μmol。
    """
    t_h = np.asarray(t_h, dtype=float)
    c_mM = np.asarray(glu_extra_mM, dtype=float)
    dt = float(np.mean(np.diff(t_h)))
    win = max(5, int(round(0.5 / dt)) | 1)   # ~0.5 h 窗口，取奇数
    poly = 2
    dc_dt_mM_per_h = savgol_filter(c_mM, window_length=win, polyorder=poly,
                                   deriv=1, delta=dt, mode="interp")
    umol_per_h = 1000.0 * V_tumor_ext_L * dc_dt_mM_per_h * float(f_leak)
    return umol_per_h


def run_neurotox_from_glu_timeseries(
    t_h: np.ndarray,
    glu_extra_mM: np.ndarray,
    V_tumor_ext_L: float,
    pk_params=None,
    tox_thr=None,
    baseline_uM: float = 50.0,
    f_leak: float = 0.05,   # 对外暴露：敏感性分析时可调整
):
    """
    用 Glu_extra(t) 推导通量 S_t(t)，驱动三室 PK，并评估神经毒性。
    返回 dict: {'t_h','S_t_umol_per_h','Cb_uM','Ct_uM','Cn_uM','tox_report'}
    """
    from .pk_toxicity import (PKParams, ToxicityThresholds,
                              simulate_three_comp_pk, assess_neurotoxicity)

    pk_params = pk_params or PKParams()
    tox_thr   = tox_thr   or ToxicityThresholds()

    S_t = _finite_diff_flux_from_conc(
        t_h=t_h, glu_extra_mM=glu_extra_mM, V_tumor_ext_L=V_tumor_ext_L, f_leak=f_leak
    )

    Cb_uM, Ct_uM, Cn_uM = simulate_three_comp_pk(
        t_grid_h=t_h, S_t_umol_per_h=S_t, params=pk_params,
        Cb0_uM=baseline_uM, Ct0_uM=baseline_uM, Cn0_uM=baseline_uM, baseline_uM=baseline_uM
    )
    report = assess_neurotoxicity(Cb_uM, t_h, tox_thr)

    return {"t_h": t_h, "S_t_umol_per_h": S_t, "Cb_uM": Cb_uM, "Ct_uM": Ct_uM,
            "Cn_uM": Cn_uM, "tox_report": report}


# --- 示例 ---
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    print("--- 1. 测试多室PK模型 ---")
    pk_model = MultiCompartmentPK()

    def infusion(t):  # 前2小时恒定输注
        return 100.0 if t <= 2.0 else 0.0

    t_pk, sol_pk = pk_model.simulate(infusion, t_end=24)
    plt.figure(figsize=(10, 6))
    plt.plot(t_pk, sol_pk[:, 0], label='Blood')
    plt.plot(t_pk, sol_pk[:, 1], label='Liver')
    plt.plot(t_pk, sol_pk[:, 2], label='Tumor')
    plt.plot(t_pk, sol_pk[:, 3], label='Other Tissues')
    plt.title('Multi-Compartment PK Model')
    plt.xlabel('Time (hours)'); plt.ylabel('Concentration (mg/L)')
    plt.legend(); plt.grid(True); plt.show()
