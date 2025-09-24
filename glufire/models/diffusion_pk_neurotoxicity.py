
from dataclasses import dataclass, asdict
from pathlib import Path
from datetime import datetime
from typing import Dict, Tuple

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ---------------------- Data classes ---------------------- #

@dataclass
class PKParams:
    """
    Three-compartment PK parameters (mouse-scale defaults).

    Units:
      - Volumes V*: L
      - Exchange/clearance k_*: h^-1
    """
    Vb: float = 0.002    # L, plasma volume ~ 2 mL
    Vt: float = 0.0005   # L, tumor interstitium ~ 0.5 mL
    Vn: float = 0.02     # L, lumped normal tissues ~ 20 mL
    k_bt: float = 1.0    # h^-1, blood <-> tumor exchange
    k_bn: float = 0.5    # h^-1, blood <-> normal tissues exchange
    k_b_clr: float = 0.5 # h^-1, plasma clearance (to baseline)
    k_t_clr: float = 0.2 # h^-1, tumor local clearance
    k_n_clr: float = 0.2 # h^-1, normal-tissue clearance


@dataclass
class ToxicityThresholds:
    """Neurotoxicity decision thresholds (µM)."""
    caution_um: float = 100.0   # caution threshold
    danger_um: float = 1000.0   # danger threshold (1 mM)


# ---------------------- Helpers ---------------------- #

def _interp_series(t_query: float, t_src: np.ndarray, y_src: np.ndarray) -> float:
    """Linear interpolation helper: y(t_query) ← (t_src, y_src)."""
    return float(np.interp(t_query, t_src, y_src))


def trapezoid_flux(t: float,
                   t_on: float, t_plateau: float, t_off: float,
                   ramp_hours: float, peak_umol_per_h: float) -> float:
    """
    Trapezoid secretion profile S_t(t) in µmol/h:
    - ramp up from t_on to t_on + ramp_hours
    - plateau until t_off
    - ramp down from t_off to t_off + ramp_hours
    Otherwise zero.
    """
    if t_on <= t < t_on + ramp_hours:
        return peak_umol_per_h * (t - t_on) / ramp_hours
    if t_on + ramp_hours <= t < t_off:
        return peak_umol_per_h
    if t_off <= t < t_off + ramp_hours:
        return peak_umol_per_h * (1 - (t - t_off) / ramp_hours)
    return 0.0


# ---------------------- Core PK ---------------------- #

def simulate_three_comp_pk(
    t_grid_h: np.ndarray,
    S_t_umol_per_h: np.ndarray,
    params: PKParams = PKParams(),
    Cb0_uM: float = 50.0,   # start from baseline for clarity on plots
    Ct0_uM: float = 50.0,
    Cn0_uM: float = 50.0,
    baseline_uM: float = 50.0,
    debug: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Three-compartment PK (blood/plasma b, tumor t, normal tissues n).
    State y = [Cb, Ct, Cn] in µM. External input S_t(t) (µmol/h) acts on the tumor compartment.

    Returns
    -------
    (Cb, Ct, Cn): tuple of arrays, each aligned to t_grid_h (µM).
    """
    t0, t1 = float(t_grid_h[0]), float(t_grid_h[-1])

    Vb, Vt, Vn = params.Vb, params.Vt, params.Vn
    kbt, kbn = params.k_bt, params.k_bn
    kbclr, ktclr, knclr = params.k_b_clr, params.k_t_clr, params.k_n_clr

    def rhs(t, y):
        Cb, Ct, Cn = y
        S_t = _interp_series(t, t_grid_h, S_t_umol_per_h)  # µmol/h

        if debug:
            if not hasattr(rhs, "_dbg"):
                rhs._dbg = 0
            if rhs._dbg < 3:
                print(f"[DBG rhs] t={t:.3f} h, S_t={S_t:.4g} µmol/h, "
                      f"Cb={Cb:.3f}, Ct={Ct:.3f}, Cn={Cn:.3f}")
                rhs._dbg += 1

        # Plasma
        dCb = (
            - kbt * (Cb - Ct)             # blood <-> tumor exchange
            - kbn * (Cb - Cn)             # blood <-> normal tissues
            - kbclr * (Cb - baseline_uM)  # clearance to baseline
        )
        # Tumor
        dCt = (
            + (Vb / Vt) * kbt * (Cb - Ct) # mass conservation exchange
            - ktclr * (Ct - baseline_uM)  # local clearance
            + S_t / Vt                    # µmol/h → µM/h
        )
        # Normal tissues
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
) -> Dict[str, float]:
    """
    Assess neurotoxicity risk from the plasma curve.

    Returns a dict with max value, cumulative time above thresholds, and flags.
    """
    Cb_uM = np.asarray(Cb_uM, dtype=float)
    t_h   = np.asarray(t_h,   dtype=float)

    Cb_max = float(np.max(Cb_uM))
    above_caution = Cb_uM >= thr.caution_um
    above_danger  = Cb_uM >= thr.danger_um

    # Approximate cumulative time above thresholds (hours) by trapezoid integration of boolean masks
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


# ---------------------- CLI runner & plotting ---------------------- #

def _build_cli():
    import argparse
    ap = argparse.ArgumentParser(
        description="Standalone neurotoxicity PK simulator (three-compartment)"
    )
    ap.add_argument("--hours", type=float, default=48.0, help="Total simulation horizon (h)")
    ap.add_argument("--dt", type=float, default=0.1, help="Time step for output grid (h)")
    ap.add_argument("--baseline", type=float, default=50.0, help="Baseline concentration for all compartments (µM)")
    ap.add_argument("--tumor-init-mM", type=float, default=30.0, help="Initial tumor glutamate (mM), worst-case")
    # Exchange/clearance overrides (optional)
    ap.add_argument("--k-bt", type=float, default=1e-3, help="Blood<->tumor exchange (h^-1), default 1e-3")
    ap.add_argument("--k-b-clr", type=float, default=0.5, help="Plasma clearance (h^-1)")
    ap.add_argument("--k-t-clr", type=float, default=0.05, help="Tumor clearance (h^-1)")
    # Optional secretion window (disabled by default when peak=0)
    ap.add_argument("--secretion-peak", type=float, default=0.0, help="Tumor secretion peak (µmol/h). 0 disables.")
    ap.add_argument("--t-on", type=float, default=8.0, help="Secretion start (h)")
    ap.add_argument("--t-off", type=float, default=34.0, help="Secretion end (h)")
    ap.add_argument("--ramp", type=float, default=2.0, help="Ramp duration for on/off (h)")
    ap.add_argument("--outdir", type=str, default="results", help="Output root directory")
    ap.add_argument("--title", type=str, default="Plasma Glutamate vs Neurotoxicity Thresholds",
                    help="Figure title")
    return ap


def _plot_and_save(t_h: np.ndarray,
                   Cb: np.ndarray, Cb_ctrl: np.ndarray,
                   baseline_uM: float,
                   thresholds: ToxicityThresholds,
                   t_on: float, t_off: float,
                   out_png: Path,
                   title: str = "Plasma Glutamate vs Neurotoxicity Thresholds") -> None:
    """Make the neurotoxicity figure and save it."""
    fig, ax = plt.subplots(figsize=(9, 5))

    # Main axis: plasma curves
    ax.plot(t_h, Cb, label="Plasma Glu (treatment)")
    ax.plot(t_h, Cb_ctrl, "--", alpha=0.9, label="Plasma Glu (control)")

    # Choose y-limits to show baseline and thresholds
    span = float(max(Cb.max(), Cb_ctrl.max()) - min(Cb.min(), Cb_ctrl.min()))
    if span < 5.0:
        pad = max(1.0, 0.15 * max(1.0, span))
        ylo = float(min(Cb.min(), Cb_ctrl.min()) - pad)
        yhi = float(max(Cb.max(), Cb_ctrl.max()) + pad)
    else:
        lo = min(Cb.min(), Cb_ctrl.min())
        hi = max(Cb.max(), Cb_ctrl.max())
        ylo, yhi = min(45, lo - 2), max(110, hi + 2)

    # Ensure the caution threshold is visible
    yhi = max(yhi, thresholds.caution_um * 1.05, 110.0)
    ax.set_ylim(ylo, yhi)

    # Baseline and thresholds
    ax.axhline(baseline_uM, ls=":", lw=1.2, alpha=0.8, label=f"Baseline {baseline_uM:.0f} µM")
    ax.axhline(thresholds.caution_um, ls="--", lw=1.2, alpha=0.6, label=f"Caution {thresholds.caution_um:.0f} µM")
    ax.axhline(thresholds.danger_um,  ls="--", lw=1.2, alpha=0.6, label="Danger 1.0 mM")

    # Right axis only for the secretion window shading; keep y synced for visual alignment
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks([])

    def _sync_ylim(_):
        ax2.set_ylim(ax.get_ylim())
    ax.callbacks.connect('ylim_changed', _sync_ylim)

    ax2.axvspan(t_on, t_off, alpha=0.08, label="Tumor secretion window")

    # Merge legends
    handles, labels = [], []
    for a in (ax, ax2):
        h, l = a.get_legend_handles_labels()
        handles += h
        labels  += l
    seen, H, L = set(), [], []
    for h, l in zip(handles, labels):
        if l in seen:
            continue
        seen.add(l); H.append(h); L.append(l)
    ax.legend(H, L, loc="upper right", framealpha=0.9)

    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Concentration (µM)")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def main():
    """CLI entry: run the PK simulation, assess thresholds, plot and save."""
    ap = _build_cli()
    args = ap.parse_args()

    sim_hours = float(args.hours)
    dt_h = float(args.dt)
    baseline_uM = float(args.baseline)
    Ct0_uM = float(args.tumor_init_mM) * 1000.0  # convert mM → µM

    # Time grid
    t_h = np.linspace(0.0, sim_hours, int(sim_hours / dt_h) + 1)

    # Optional trapezoid secretion profile (disabled if peak=0)
    if float(args.secretion_peak) > 0.0:
        S_t = np.array([
            trapezoid_flux(t,
                           t_on=float(args.t_on),
                           t_plateau=float(args.t_on) + float(args.ramp),
                           t_off=float(args.t_off),
                           ramp_hours=float(args.ramp),
                           peak_umol_per_h=float(args.secretion_peak))
            for t in t_h
        ], dtype=float)
    else:
        S_t = np.zeros_like(t_h)

    # PK parameters (use your tuned values: weak blood<->tumor coupling + slower tumor clearance)
    params = PKParams(
        k_bt=float(args.k_bt),     # default 1e-3 to avoid instant washout of 30 mM
        k_t_clr=float(args.k_t_clr),
        k_b_clr=float(args.k_b_clr),
        # others keep defaults: Vb=2 mL, Vt=0.5 mL, k_bn=0.5, k_n_clr=0.2
    )

    # Treatment: plasma starts 50 µM; tumor starts at 30 mM worst-case; no extra secretion by default
    Cb, Ct, Cn = simulate_three_comp_pk(
        t_grid_h=t_h,
        S_t_umol_per_h=S_t,
        params=params,
        Cb0_uM=baseline_uM,
        Ct0_uM=Ct0_uM,
        Cn0_uM=baseline_uM,
        baseline_uM=baseline_uM,
        # debug=True,
    )

    # Control: all compartments start at baseline; no secretion
    Cb_ctrl, Ct_ctrl, Cn_ctrl = simulate_three_comp_pk(
        t_grid_h=t_h,
        S_t_umol_per_h=np.zeros_like(S_t),
        params=params,
        Cb0_uM=baseline_uM,
        Ct0_uM=baseline_uM,
        Cn0_uM=baseline_uM,
        baseline_uM=baseline_uM,
    )

    # Numeric check of elevation (not just visually)
    dCb = Cb - baseline_uM
    dCt = Ct - baseline_uM
    print(f"[DBG] ΔCb_max = {float(np.max(dCb)):.6f} µM, ΔCt_max = {float(np.max(dCt)):.6f} µM")

    # Risk assessment
    thr = ToxicityThresholds()
    report_treat = assess_neurotoxicity(Cb, t_h, thr)
    report_ctrl  = assess_neurotoxicity(Cb_ctrl, t_h, thr)

    print("\n=== Neurotoxicity Risk Report (treatment) ===")
    for k, v in report_treat.items():
        print(f"{k}: {v}")

    print("\n=== Neurotoxicity Risk Report (control) ===")
    for k, v in report_ctrl.items():
        print(f"{k}: {v}")

    # Output path (timestamped folder; no "The figure is saved to:" line)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = Path(args.outdir) / f"neurotoxicity_{stamp}"
    out_png = outdir / "plasma_glu_neurotoxicity.png"

    print(f"\n[INFO] Writing figure -> {out_png}")
    _plot_and_save(
        t_h=t_h,
        Cb=Cb, Cb_ctrl=Cb_ctrl,
        baseline_uM=baseline_uM,
        thresholds=thr,
        t_on=float(args.t_on), t_off=float(args.t_off),
        out_png=out_png,
        title=args.title,
    )

    print("[INFO] Done.")
    print("[INFO] Params:", asdict(params))


if __name__ == "__main__":
    main()
