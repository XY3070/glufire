# 04_diffusion_pde.py
# 1D 球对称扩散-反应（老版） + 新增：Robin（血管壁交换） + Robin+一室血浆PK
# 运行后自动做三套情景，并始终保存图与CSV。
# - legacy：复现旧逻辑（Dirichlet 边界），文件名与旧版一致，保证兼容
# - robin：Robin 边界（计算进血通量）
# - robin_pk：Robin + 一室血浆PK（给出血浆浓度-时间曲线）

import numpy as np
import matplotlib.pyplot as plt
import argparse

# ---------------------------
# 工具函数
# ---------------------------
def rho_bac_profile(r, R, core_frac=0.4, core_density=1.0, shell_density=0.2):
    """单位化细菌密度：核心区高、外壳低。"""
    return np.where(r <= core_frac * R, core_density, shell_density)

def average_rho(core_frac=0.4, core_density=1.0, shell_density=0.2):
    """球体体积加权平均密度（用于把整体浓度增长换算为源强）。"""
    f = np.clip(core_frac, 0.0, 1.0)
    return (f**3) * core_density + (1 - f**3) * shell_density

def calibrate_q_from_bulk_delta(delta_mM, duration_h, core_frac, core_density, shell_density):
    """
    用“整体平均浓度 duration_h 小时增长 delta_mM”来标定 q_glu（上限、保守）。
    忽略扩散/清除，只做量级锚定。
    """
    # TODO(需实验数据)：当你们有“上清/总体浓度随时间”的数据时，用这个来估 q_glu
    if duration_h <= 0:
        raise ValueError("duration_h must be > 0")
    S_avg = delta_mM / (duration_h * 3600.0)  # mM/s
    q = S_avg / max(average_rho(core_frac, core_density, shell_density), 1e-12)
    return q

# ---------------------------
# 数值核心：单次情景模拟
# ---------------------------
def simulate_spherical(
    mode,                     # "legacy" | "robin" | "robin_pk"
    D=5e-6,                   # TODO(需实验数据)：有效扩散系数 [mm^2/s]
    kdeg=0.0,                 # TODO(需实验数据)：一阶清除(摄取/降解) [1/s]
    R=0.5,                    # 球半径 [mm]
    Nr=200,
    T_s=6*3600,               # 总时间 [s]
    q_glu=1e-7,               # TODO(需实验数据)：单位密度源强 [mM/s]
    core_frac=0.4,
    core_density=1.0,
    shell_density=0.2,
    h_mm_per_s=1e-3,          # TODO(需实验数据)：Robin 质传系数 [mm/s]
    Vd_L=5.0,                 # TODO(需实验数据/物种)：分布容积 [L]
    t_half_s=1800.0,          # TODO(需实验数据)：血浆半衰期 [s]
    plasma_baseline_mM=0.0,   # TODO(可填)：血浆基线 [mM]，成人常见 0.05~0.1 mM
    sample_every=100
):
    """
    返回：r, C_end, (times, flux_mmol_s), (times, C_plasma) 视情景而定
    """
    assert mode in ("legacy", "robin", "robin_pk")
    dr = R / Nr
    dt = 0.2 * dr**2 / D    # 显式稳定步长（保守）
    Nt = int(T_s / dt)

    r = np.linspace(0.0, R, Nr + 1)
    C = np.zeros_like(r)
    C_new = C.copy()

    rho = rho_bac_profile(r, R, core_frac, core_density, shell_density)

    area_mm2 = 4.0 * np.pi * R**2
    C_plasma = float(plasma_baseline_mM)
    k_el = np.log(2.0) / t_half_s if t_half_s > 0 else 0.0

    times = []
    flux_mmol_s = []
    plasma_mM = []

    use_robin = (mode in ("robin", "robin_pk"))
    use_pk    = (mode == "robin_pk")

    # Biot 数（仅供打印）
    if use_robin:
        Bi = h_mm_per_s * R / D
        print(f"[{mode}] Biot = h*R/D = {Bi:.3g}")

    for n in range(Nt):
        # 内部点：球坐标 Laplacian + 源项 + 清除
        for i in range(1, Nr):
            rp = r[i] + 1e-12
            d2 = (C[i+1] - 2*C[i] + C[i-1]) / dr**2
            d1 = (C[i+1] - C[i-1]) / (2*dr)
            lap = d2 + 2.0/rp * d1
            S = q_glu * rho[i]                     # mM/s
            C_new[i] = C[i] + dt * (D*lap - kdeg*C[i] + S)

        # r=0：对称 Neumann
        C_new[0] = C_new[1]

        # r=R：边界
        if use_robin:
            # Robin: -D*(C_N - C_{N-1})/dr = h*(C_N - C_plasma)
            denom = h_mm_per_s + D/dr
            C_R = (h_mm_per_s * C_plasma - (D/dr) * C_new[-2]) / denom
            C_new[-1] = max(C_R, 0.0)
        else:
            # 旧版：Dirichlet sink
            C_new[-1] = 0.0

        # 防止数值小负
        C_new[:] = np.maximum(C_new, 0.0)
        C, C_new = C_new, C

        # Robin 情景：通量 & 血浆
        if use_robin:
            j = h_mm_per_s * (C[-1] - C_plasma)        # mM*mm/s
            Rin = j * area_mm2 * 1e-6                  # mmol/s  (1 mm^3 = 1e-6 L)

            if use_pk:
                dCp = (Rin / Vd_L) - k_el * C_plasma   # mM/s
                C_plasma = max(C_plasma + dt * dCp, 0.0)

            if n % sample_every == 0:
                times.append(n*dt)
                flux_mmol_s.append(Rin)
                if use_pk:
                    plasma_mM.append(C_plasma)

    out_flux = (np.array(times), np.array(flux_mmol_s)) if use_robin else None
    out_plasma = (np.array(times), np.array(plasma_mM)) if use_pk else None
    return r, C, out_flux, out_plasma

# ---------------------------
# 统一画图/保存（始终出图）
# ---------------------------
def save_profile(prefix, r, C):
    np.savetxt(f"{prefix}.csv", np.c_[r, C], delimiter=",",
               header="radius_mm,glu_mM", comments="")
    plt.figure()
    plt.plot(r, C)
    plt.xlabel("radius (mm)")
    plt.ylabel("Glu (mM)")
    plt.title(prefix.replace("_", " "))
    plt.tight_layout()
    plt.savefig(f"{prefix}.png", dpi=200)
    plt.close()

def save_flux(prefix, times, flux):
    np.savetxt(f"{prefix}.csv", np.c_[times, flux], delimiter=",",
               header="time_s,flux_mmol_per_s", comments="")
    plt.figure()
    plt.plot(times, flux)
    plt.xlabel("time (s)")
    plt.ylabel("Flux into blood (mmol/s)")
    plt.title(prefix.replace("_", " "))
    plt.tight_layout()
    plt.savefig(f"{prefix}.png", dpi=200)
    plt.close()

def save_plasma(prefix, times, Cp):
    np.savetxt(f"{prefix}.csv", np.c_[times, Cp], delimiter=",",
               header="time_s,plasma_glu_mM", comments="")
    plt.figure()
    plt.plot(times, Cp)
    plt.xlabel("time (s)")
    plt.ylabel("Plasma Glu (mM)")
    plt.title(prefix.replace("_", " "))
    plt.tight_layout()
    plt.savefig(f"{prefix}.png", dpi=200)
    plt.close()

# ---------------------------
# 主入口：一次运行做三套情景
# ---------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Spherical diffusion with legacy (Dirichlet), Robin exchange, and Robin+plasma PK."
    )
    # ——基础物理参数（有默认值，可被命令行覆盖）——
    ap.add_argument("--D", type=float, default=5e-6, help="有效扩散 D [mm^2/s]  # TODO(需实验数据)")
    ap.add_argument("--kdeg", type=float, default=0.0, help="一阶清除 kdeg [1/s]  # TODO(需实验数据)")
    ap.add_argument("--R", type=float, default=0.5, help="球半径 R [mm]")
    ap.add_argument("--Nr", type=int, default=200, help="径向网格数")
    ap.add_argument("--T-h", type=float, default=6.0, help="总时间 T [小时]")

    # ——源项/密度——
    ap.add_argument("--q-glu", type=float, default=1e-7, help="源强 q_glu [mM/s/密度单位]  # TODO(需实验数据)")
    ap.add_argument("--core-frac", type=float, default=0.4, help="核心半径占比 (0~1)")
    ap.add_argument("--core-density", type=float, default=1.0, help="核心密度系数")
    ap.add_argument("--shell-density", type=float, default=0.2, help="外壳密度系数")

    # ——可选：用总体浓度增长来标定 q_glu（上限保守估计）——
    ap.add_argument("--calibrate-total-mM", type=float, default=None,
                    help="若提供：用整体平均浓度增长 ΔC[mM] 在 duration 内标定 q_glu")
    ap.add_argument("--calibrate-duration-h", type=float, default=None,
                    help="与 --calibrate-total-mM 配套的时长 [小时]")

    # ——Robin/PK 参数（无开关，脚本自动跑 robin 与 robin_pk）——
    ap.add_argument("--h", type=float, default=1e-3, help="Robin 质传系数 h [mm/s]  # TODO(需实验数据)")
    ap.add_argument("--Vd-L", type=float, default=5.0, help="分布容积 Vd [L]  # TODO(需实验数据/物种)")
    ap.add_argument("--t-half-s", type=float, default=1800.0, help="血浆半衰期 t1/2 [s]  # TODO(需实验数据)")
    ap.add_argument("--plasma-baseline-mM", type=float, default=0.0, help="血浆基线浓度 [mM] (可设0.05~0.1)")

    args = ap.parse_args()

    # 时间统一
    T_s = args.T_h * 3600.0

    # 若给了总体浓度增长 → 标定 q_glu
    q_glu = args.q_glu
    if args.calibrate_total_mM is not None and args.calibrate_duration_h is not None:
        q_glu = calibrate_q_from_bulk_delta(
            delta_mM=args.calibrate_total_mM,
            duration_h=args.calibrate_duration_h,
            core_frac=args.core_frac,
            core_density=args.core_density,
            shell_density=args.shell_density
        )
        print(f"[calibrate] q_glu = {q_glu:.3e} mM/s per density unit "
              f"(from ΔC={args.calibrate_total_mM} mM over {args.calibrate_duration_h} h)")

    # ===== 情景1：legacy（旧版，保持原文件名） =====
    r, C_end, _, _ = simulate_spherical(
        mode="legacy",
        D=args.D, kdeg=args.kdeg, R=args.R, Nr=args.Nr, T_s=T_s,
        q_glu=q_glu, core_frac=args.core_frac, core_density=args.core_density, shell_density=args.shell_density
    )
    # 保留旧文件名，保证兼容
    save_profile("spheroid_glu_profile", r, C_end)
    print("Saved spheroid_glu_profile.csv/png (legacy Dirichlet)")

    # ===== 情景2：robin（血管壁交换） =====
    r, C_end, flux, _ = simulate_spherical(
        mode="robin",
        D=args.D, kdeg=args.kdeg, R=args.R, Nr=args.Nr, T_s=T_s,
        q_glu=q_glu, core_frac=args.core_frac, core_density=args.core_density, shell_density=args.shell_density,
        h_mm_per_s=args.h
    )
    save_profile("spheroid_glu_profile_robin", r, C_end)
    if flux is not None:
        t_f, f_f = flux
        save_flux("boundary_flux_into_blood", t_f, f_f)
    print("Saved spheroid_glu_profile_robin.csv/png and boundary_flux_into_blood.csv/png")

    # ===== 情景3：robin_pk（血管交换 + 一室血浆） =====
    r, C_end, flux, plasma = simulate_spherical(
        mode="robin_pk",
        D=args.D, kdeg=args.kdeg, R=args.R, Nr=args.Nr, T_s=T_s,
        q_glu=q_glu, core_frac=args.core_frac, core_density=args.core_density, shell_density=args.shell_density,
        h_mm_per_s=args.h, Vd_L=args.Vd_L, t_half_s=args.t_half_s, plasma_baseline_mM=args.plasma_baseline_mM
    )
    save_profile("spheroid_glu_profile_robin_pk", r, C_end)
    if flux is not None:
        t_f, f_f = flux
        save_flux("boundary_flux_into_blood_robin_pk", t_f, f_f)
    if plasma is not None:
        t_p, Cp = plasma
        save_plasma("plasma_concentration", t_p, Cp)
    print("Saved spheroid_glu_profile_robin_pk.csv/png and plasma_concentration.csv/png")

if __name__ == "__main__":
    main()
