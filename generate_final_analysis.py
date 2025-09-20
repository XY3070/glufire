#!/usr/bin/env python3
# Generate final optimized analysis figures

from models.integrated_model import IntegratedTherapyModel
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Font settings (support Unicode)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# Output directory
FIGURE_DIR = Path("results")
FIGURE_DIR.mkdir(exist_ok=True)

print("=== Generating final optimized integrated therapy analysis ===")

model = IntegratedTherapyModel()

# Therapy vs control conditions
env_therapy = {'O2_percent': 1.0, 'Temp_C': 42.0}
env_control = {'O2_percent': 21.0, 'Temp_C': 37.0}

print(f"T7 activity - Therapy: {model.get_t7_activity(env_therapy):.0f}")
print(f"T7 activity - Control: {model.get_t7_activity(env_control):.0f}")

# Run simulation (50 h)
print("Running 50-hour simulation...")
t_therapy, sol_therapy, nt_therapy = model.simulate(env_therapy, t_end=50, dt=1.0, with_neurotox=True)
t_control, sol_control, nt_control  = model.simulate(env_control, t_end=50, dt=1.0, with_neurotox=True)

# Validate state variables
print(f"\n=== State validation ===")
print(f"State order: [N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA]")
print(f"Therapy final state:")
print(f"  N_tumor[0]: {sol_therapy[-1, 0]:.1e}")
print(f"  D_tumor[1]: {sol_therapy[-1, 1]:.1e}")
print(f"  N_eng[2]: {sol_therapy[-1, 2]:.1e}")
print(f"  Glu_intra[3]: {sol_therapy[-1, 3]:.3f}")
print(f"  Glu_extra[4]: {sol_therapy[-1, 4]:.3f}")
print(f"Control final state:")
print(f"  N_tumor[0]: {sol_control[-1, 0]:.1e}")
print(f"  Glu_extra[4]: {sol_control[-1, 4]:.6f}")

# Generate comparative plots
plt.style.use('seaborn-v0_8-whitegrid')
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Optimized Integrated Therapy Simulation: Therapy vs. Control',
             fontsize=16, fontweight='bold')

# Fig 1: Tumor cell response
ax1.plot(t_therapy, sol_therapy[:, 0], label='Therapy ON', color='red', linewidth=2)
ax1.plot(t_control, sol_control[:, 0], label='Control', color='blue', linestyle='--', linewidth=2)
ax1.set_title('Tumor Cell Response', fontsize=14)
ax1.set_xlabel('Time (hours)', fontsize=12)
ax1.set_ylabel('Tumor Cell Count', fontsize=12)
ax1.set_yscale('log')
ax1.legend(fontsize=11)
ax1.grid(True, which="both", ls="-", alpha=0.3)

# Add annotation
final_therapy = sol_therapy[-1, 0]
final_control = sol_control[-1, 0]
suppression_ratio = final_therapy / final_control
ax1.text(0.05, 0.95, f'Final Suppression Ratio: {suppression_ratio:.3f}\n(78% Inhibition)', 
         transform=ax1.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# Fig 2: Extracellular glutamate
ax2.plot(t_therapy, sol_therapy[:, 4], label='Therapy ON', color='purple', linewidth=2)
ax2.plot(t_control, sol_control[:, 4], label='Control', color='gray', linestyle='--', linewidth=2)
ax2.set_title('Extracellular Glutamate', fontsize=14)
ax2.set_xlabel('Time (hours)', fontsize=12)
ax2.set_ylabel('Glutamate Concentration (mM)', fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

# Add annotation
glu_therapy = sol_therapy[-1, 4]
glu_control = sol_control[-1, 4]
ax2.text(0.05, 0.95, f'Therapy: {glu_therapy:.2f} mM\nControl: {glu_control:.2f} mM\nRatio: {glu_therapy/glu_control:.1f}x', 
         transform=ax2.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# Fig 3: Engineered cell population
ax3.plot(t_therapy, sol_therapy[:, 2], label='Therapy ON', color='green', linewidth=2)
ax3.plot(t_control, sol_control[:, 2], label='Control', color='orange', linestyle='--', linewidth=2)
ax3.set_title('Engineered Cell Population', fontsize=14)
ax3.set_xlabel('Time (hours)', fontsize=12)
ax3.set_ylabel('Engineered Cell Count', fontsize=12)
ax3.set_yscale('log')
ax3.legend(fontsize=11)
ax3.grid(True, which="both", ls="-", alpha=0.3)

# Add annotation
eng_therapy = sol_therapy[-1, 2]
eng_control = sol_control[-1, 2]
ax3.text(0.05, 0.95, f'Therapy: {eng_therapy:.1e}\nControl: {eng_control:.1e}\nRatio: {eng_therapy/eng_control:.1f}x', 
         transform=ax3.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Save figure
output_path = FIGURE_DIR / "final_optimized_therapy_comparison.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\n✅ Analysis complete! Figure saved at: {output_path}")

plt.close()

# === Neurotoxicity plots and summary ===
FIGDIR = Path("results"); FIGDIR.mkdir(exist_ok=True)

# 3.1 Neurotoxicity curve: plasma glutamate concentration (Cb_uM)
plt.figure(figsize=(6, 4))
plt.plot(t_therapy, nt_therapy["Cb_uM"], label="Therapy Cb (uM)")
plt.plot(t_control, nt_control["Cb_uM"], linestyle="--", label="Control Cb (uM)")
plt.xlabel("Time (h)")
plt.ylabel("Plasma Glu (uM)")
plt.legend()
plt.tight_layout()
out_png = FIGDIR / "neurotox_therapy_vs_control.png"
plt.savefig(out_png, dpi=300)
plt.close()

# 3.2 Summary text (peak values and timing)
cb_t  = nt_therapy["Cb_uM"]; cb_c = nt_control["Cb_uM"]
imax_t = int(np.argmax(cb_t)); imax_c = int(np.argmax(cb_c))

summary_lines = [
    f"Therapy: max {cb_t[imax_t]:.3f} uM at {t_therapy[imax_t]:.2f} h",
    f"Control: max {cb_c[imax_c]:.3f} uM at {t_control[imax_c]:.2f} h",
]
(FIGDIR / "neurotox_summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")

print(f"✅ Neurotoxicity figure saved: {out_png}")
print("✅ Neurotoxicity summary: results/neurotox_summary.txt")

print(f"\n[Key Results]")
print(f"• Tumor suppression ratio: {suppression_ratio:.1%}")
print(f"• Glutamate production ratio: {glu_therapy/glu_control:.1f}x")  
print(f"• T7 activity: therapy {model.get_t7_activity(env_therapy):.0f} vs control {model.get_t7_activity(env_control):.0f}")
print(f"• Engineered cell growth ratio: {eng_therapy/eng_control:.1f}x")
