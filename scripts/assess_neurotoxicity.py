# --- path bootstrap ---
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
# ----------------------

# scripts/assess_neurotoxicity.py
import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

# === 只为定位，不改模型逻辑 ===
import models.pk_toxicity as _pkt
print("[USING PK]", _pkt.__file__)  # 跑起来会打印真正加载的 pk_toxicity.py 路径
# ============================

try:
    from models.pk_toxicity import (
        PKParams, ToxicityThresholds,
        simulate_three_comp_pk, assess_neurotoxicity
    )
except Exception as e:
    print("[ERROR] failed to import models.pk_toxicity:", e, flush=True)
    raise


def trapezoid_flux(t, t_on, t_plateau, t_off, ramp_hours, peak_umol_per_h):
    """
    2 段线性渐变 + 平台的梯形通量（单位 μmol/h）：
    - 从 t_on 开始，用 ramp_hours 线性上升到 peak；
    - 到 t_plateau 进入平台，保持到 t_off；
    - 之后用 ramp_hours 线性下降到 0。
    """
    # 上升段
    if t_on <= t < t_on + ramp_hours:
        return peak_umol_per_h * (t - t_on) / ramp_hours
    # 平台段
    if t_on + ramp_hours <= t < t_off:
        if t < t_plateau:
            return peak_umol_per_h
        if t_plateau <= t < t_off:
            return peak_umol_per_h
    # 下降段
    if t_off <= t < t_off + ramp_hours:
        return peak_umol_per_h * (1 - (t - t_off) / ramp_hours)
    # 其他时间为 0
    return 0.0


def main():
    print("[INFO] preparing timeline & physiologic-ish flux...", flush=True)
    # ==== 配置区：可改 ====
    sim_hours = 48.0
    dt_h = 0.1
    baseline_uM = 50.0             # 血/其他室的基线
    intensity = "init-30mM"        # 标记一下这次情形
    peak_umol_per_h = 0.0          # 不额外加外泌通量，单看30 mM初始条件的外溢效应
    # 通量时序（保留以便之后要用S_t时再开）
    t_on = 8.0
    ramp_hours = 2.0
    t_plateau = t_on + ramp_hours
    t_off = 34.0
    # =====================

    # 时间网格
    t_h = np.linspace(0.0, sim_hours, int(sim_hours / dt_h) + 1)

    # 构造理想化梯形通量（作用在肿瘤室，单位 μmol/h）
    S_t = np.array([
        trapezoid_flux(t, t_on=t_on, t_plateau=t_plateau, t_off=t_off,
                       ramp_hours=ramp_hours, peak_umol_per_h=peak_umol_per_h)
        for t in t_h
    ], dtype=float)

    print(
        "[CHK] S_t stats -> min={:.3g} μmol/h, max={:.3g} μmol/h, nonzero={} / {}".format(
            float(S_t.min()), float(S_t.max()), int((S_t != 0).sum()), S_t.size
        ),
        flush=True
    )

    print("[INFO] running PK integration...", flush=True)
    # 关键：把血↔肿瘤耦合调小，避免“瞬时倾倒”；肿瘤清除也稍小一些以维持高浓度
    params = PKParams(
        k_bt=1e-3,     # 血<->肿瘤交换很弱（每小时 0.1% 级），避免瞬间洗出30 mM
        k_t_clr=0.05,  # 肿瘤清除放缓（从30 mM往50 μM回落的速度别太快）
        # 其余参数用默认：Vb=2 mL, Vt=0.5 mL, k_bn=0.5, k_b_clr=0.5, k_n_clr=0.2 ...
    )

    # 三个室的初值：血和其他维持50 μM，肿瘤从30 mM开始
    Cb, Ct, Cn = simulate_three_comp_pk(
        t_grid_h=t_h,
        S_t_umol_per_h=S_t * 0.0,      # 此次不叠加通量；若想再加外泌，改回 S_t
        params=params,
        Cb0_uM=baseline_uM,
        Ct0_uM=30000.0,                # 30 mM = 30000 μM
        Cn0_uM=baseline_uM,
        baseline_uM=baseline_uM,
        # debug=True,
    )

    # === control 组（基线起步 + 无外泌）===
    Cb_ctrl, Ct_ctrl, Cn_ctrl = simulate_three_comp_pk(
        t_grid_h=t_h,
        S_t_umol_per_h=S_t * 0.0,      # 对照为 0
        params=params,
        Cb0_uM=baseline_uM,
        Ct0_uM=baseline_uM,            # 肿瘤室也从 50 µM 起步
        Cn0_uM=baseline_uM,
        baseline_uM=baseline_uM,
    )

    # 数值上确认真实抬升（别只看图）
    dCb = Cb - baseline_uM
    dCt = Ct - baseline_uM
    print(f"[DBG] ΔCb_max = {float(np.max(dCb)):.6f} μM, ΔCt_max = {float(np.max(dCt)):.6f} μM")

    print("[INFO] assessing risk...", flush=True)
    thr = ToxicityThresholds()
    report = assess_neurotoxicity(Cb, t_h, thr)
    report_ctrl = assess_neurotoxicity(Cb_ctrl, t_h, thr)

    print("\n=== Neurotoxicity Risk Report (physio-ish) ===", flush=True)
    for k, v in report.items():
        print(f"{k}: {v}", flush=True)

    print("\n=== Neurotoxicity Risk Report (control) ===", flush=True)
    for k, v in report_ctrl.items():
        print(f"{k}: {v}", flush=True)

    # 输出目录（带时间戳与强度，避免覆盖）
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = Path("results") / f"neurotoxicity_{timestamp}_{intensity}"
    outdir.mkdir(parents=True, exist_ok=True)
    out_png = outdir / "plasma_glu_neurotoxicity.png"

    # ===== 画图：主轴画曲线+阈值；右轴只画分泌窗口并与左轴同步 =====
    print(f"[INFO] saving figure -> {out_png}", flush=True)
    fig, ax = plt.subplots(figsize=(9, 5))

    # 主轴：血浆
    ax.plot(t_h, Cb, label="Plasma Glu (therapy)")
    ax.plot(t_h, Cb_ctrl, "--", alpha=0.9, label="Plasma Glu (control)")

    # 根据治疗/对照的总体范围设置y轴
    span = float(max(Cb.max(), Cb_ctrl.max()) - min(Cb.min(), Cb_ctrl.min()))
    if span < 5.0:
        pad = max(1.0, 0.15 * max(1.0, span))
        ylo = float(min(Cb.min(), Cb_ctrl.min()) - pad)
        yhi = float(max(Cb.max(), Cb_ctrl.max()) + pad)
    else:
        lo = min(Cb.min(), Cb_ctrl.min())
        hi = max(Cb.max(), Cb_ctrl.max())
        ylo, yhi = min(45, lo - 2), max(110, hi + 2)

    # 确保能看到 caution=100 μM（必要时抬高上限）
    yhi = max(yhi, thr.caution_um * 1.05, 110)  # danger 若要可见可改为 thr.danger_um*1.05
    ax.set_ylim(ylo, yhi)

    # 基线 & 轴样式
    ax.axhline(baseline_uM, ls=":", lw=1.2, alpha=0.8, label=f"Baseline {baseline_uM:.0f} μM")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Concentration (μM)")
    ax.grid(True, alpha=0.25)

    # —— 阈值直接画在主轴上（与曲线同一坐标系，避免错位）——
    ax.axhline(thr.caution_um, ls="--", lw=1.2, alpha=0.6,
               label=f"Caution {thr.caution_um:.0f} μM")
    ax.axhline(thr.danger_um,  ls="--", lw=1.2, alpha=0.6,
               label="Danger 1.0 mM")

    # 右轴仅用于分泌窗口阴影；与左轴范围同步，避免视觉错位
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks([])

    def _sync_ylim(_):
        ax2.set_ylim(ax.get_ylim())
    ax.callbacks.connect('ylim_changed', _sync_ylim)

    ax2.axvspan(t_on, t_off, alpha=0.08, label="Tumor secretion window")

    # 合并图例
    handles, labels = [], []
    for a in (ax, ax2):
        h, l = a.get_legend_handles_labels()
        handles += h; labels += l
    seen, H, L = set(), [], []
    for h, l in zip(handles, labels):
        if l in seen:
            continue
        seen.add(l); H.append(h); L.append(l)
    ax.legend(H, L, loc="upper right", framealpha=0.9)

    ax.set_title(
        f"Plasma Glutamate vs Neurotoxicity Thresholds\n"
        f"(Intensity={intensity}, peak={peak_umol_per_h} μmol/h)"
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print("[INFO] done.", flush=True)


if __name__ == "__main__":
    print("Running assess_neurotoxicity.py ...", flush=True)
    try:
        main()
    except Exception as e:
        print("[FATAL] script crashed:", e, flush=True)
        raise
    print("Done.", flush=True)
