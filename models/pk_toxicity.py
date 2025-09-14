# models/pk_toxicity.py
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp


# ---------------------- Data classes ---------------------- #

@dataclass
class PKParams:
    """
    三室PK参数（默认数值按小鼠量级）。
    单位：
      - 体积 V*: L
      - 交换/清除 k_*: h^-1
    """
    Vb: float = 0.002    # L, 血浆体积 ~2 mL
    Vt: float = 0.0005   # L, 肿瘤室体积 0.5 mL
    Vn: float = 0.02     # L, 其他组织等效体积 20 mL
    k_bt: float = 1.0    # h^-1, 血<->肿瘤 交换系数
    k_bn: float = 0.5    # h^-1, 血<->其他组织 交换系数
    k_b_clr: float = 0.5 # h^-1, 血浆清除
    k_t_clr: float = 0.2 # h^-1, 肿瘤清除
    k_n_clr: float = 0.2 # h^-1, 其他组织清除


@dataclass
class ToxicityThresholds:
    """神经毒性判定阈值（单位 μM）。"""
    caution_um: float = 100.0   # 提示阈值
    danger_um: float = 1000.0   # 危险阈值（1 mM）


# ---------------------- Helpers ---------------------- #

def _interp_series(t_query: float, t_src: np.ndarray, y_src: np.ndarray) -> float:
    """线性插值辅助：y(t_query) ← (t_src, y_src)。"""
    return float(np.interp(t_query, t_src, y_src))


# ---------------------- Core PK ---------------------- #

def simulate_three_comp_pk(
    t_grid_h: np.ndarray,
    S_t_umol_per_h: np.ndarray,
    params: PKParams = PKParams(),
    Cb0_uM: float = 50.0,   # 从基线起步，图更直观
    Ct0_uM: float = 50.0,
    Cn0_uM: float = 50.0,
    baseline_uM: float = 50.0,
    debug: bool = False,
):
    """
    三室PK（血浆 b、肿瘤 t、其他组织 n）。
    y = [Cb, Ct, Cn]（μM）；外源输入 S_t(t)（μmol/h）施加在肿瘤室。
    """
    t0, t1 = float(t_grid_h[0]), float(t_grid_h[-1])

    Vb, Vt, Vn = params.Vb, params.Vt, params.Vn
    kbt, kbn = params.k_bt, params.k_bn
    kbclr, ktclr, knclr = params.k_b_clr, params.k_t_clr, params.k_n_clr

    def rhs(t, y):
        Cb, Ct, Cn = y
        S_t = _interp_series(t, t_grid_h, S_t_umol_per_h)  # μmol/h

        if debug:
            if not hasattr(rhs, "_dbg"): rhs._dbg = 0
            if rhs._dbg < 3:
                print(f"[DBG rhs] t={t:.3f} h, S_t={S_t:.4g} μmol/h, Cb={Cb:.3f}, Ct={Ct:.3f}, Cn={Cn:.3f}")
                rhs._dbg += 1

        # 血浆
        dCb = (
            - kbt * (Cb - Ct)             # 血<->肿瘤
            - kbn * (Cb - Cn)             # 血<->其他
            - kbclr * (Cb - baseline_uM)  # 清除到基线
        )
        # 肿瘤
        dCt = (
            + (Vb / Vt) * kbt * (Cb - Ct) # 量守恒对应的反向交换
            - ktclr * (Ct - baseline_uM)  # 清除
            + S_t / Vt                    # μmol/h -> μM/h
        )
        # 其他组织
        dCn = (
            + (Vb / Vn) * kbn * (Cb - Cn)
            - knclr * (Cn - baseline_uM)
        )
        return [dCb, dCt, dCn]

    sol = solve_ivp(
        rhs, (t0, t1), [Cb0_uM, Ct0_uM, Cn0_uM],
        t_eval=t_grid_h, method="LSODA", rtol=1e-6, atol=1e-9
    )
    if not sol.success:
        raise RuntimeError(f"PK ODE failed: {sol.message}")

    Cb, Ct, Cn = sol.y
    return Cb, Ct, Cn


# ---------------------- Toxicity assessment ---------------------- #

def assess_neurotoxicity(
    Cb_uM: np.ndarray,
    t_h: np.ndarray,
    thr: ToxicityThresholds = ToxicityThresholds()
):
    """
    基于血浆曲线评估神经毒性风险。
    返回字典含：最大值、超过各阈值的累计时长、布尔标记等。
    """
    Cb_uM = np.asarray(Cb_uM, dtype=float)
    t_h   = np.asarray(t_h,   dtype=float)

    Cb_max = float(np.max(Cb_uM))
    above_caution = Cb_uM >= thr.caution_um
    above_danger  = Cb_uM >= thr.danger_um

    # 用时长积分近似累计超阈时间（单位：小时）
    t_above_caution_h = float(np.trapz(above_caution.astype(float), t_h))
    t_above_danger_h  = float(np.trapz(above_danger.astype(float),  t_h))

    return {
        "Cb_max_uM": Cb_max,
        "time_above_caution_h": t_above_caution_h,
        "time_above_danger_h":  t_above_danger_h,
        "flag_caution": bool(np.any(above_caution)),
        "flag_danger":  bool(np.any(above_danger)),
        "caution_threshold_uM": thr.caution_um,
        "danger_threshold_uM":  thr.danger_um,
    }
