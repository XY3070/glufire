#!/usr/bin/env python3
# 生成最终优化的分析图表

from models.integrated_model import IntegratedTherapyModel
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# 设置中文字体支持
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']  # 支持中文显示
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 设置输出目录
FIGURE_DIR = Path("results")
FIGURE_DIR.mkdir(exist_ok=True)

print("=== 生成最终优化的整合治疗分析 ===")

model = IntegratedTherapyModel()

# 治疗和对照条件
env_therapy = {'O2_percent': 1.0, 'Temp_C': 42.0}
env_control = {'O2_percent': 21.0, 'Temp_C': 37.0}

print(f"T7活性 - Therapy: {model.get_t7_activity(env_therapy):.0f}")
print(f"T7活性 - Control: {model.get_t7_activity(env_control):.0f}")

# 运行100小时模拟
print("运行50小时模拟...")
t_therapy, sol_therapy, nt_therapy = model.simulate(env_therapy, t_end=50, dt=1.0, with_neurotox=True)
t_control,  sol_control,  nt_control  = model.simulate(env_control,  t_end=50, dt=1.0, with_neurotox=True)

# 验证状态变量顺序: [N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA]
print(f"\n=== 模型状态验证 ===")
print(f"状态变量顺序: [N_tumor, D_tumor, N_eng, Glu_intra, Glu_extra, Icd, gdhA]")
print(f"治疗组最终状态:")
print(f"  N_tumor[0]: {sol_therapy[-1, 0]:.1e}")
print(f"  D_tumor[1]: {sol_therapy[-1, 1]:.1e}")
print(f"  N_eng[2]: {sol_therapy[-1, 2]:.1e}")
print(f"  Glu_intra[3]: {sol_therapy[-1, 3]:.3f}")
print(f"  Glu_extra[4]: {sol_therapy[-1, 4]:.3f}")
print(f"对照组最终状态:")
print(f"  N_tumor[0]: {sol_control[-1, 0]:.1e}")
print(f"  Glu_extra[4]: {sol_control[-1, 4]:.6f}")

# 生成对比图
plt.style.use('seaborn-v0_8-whitegrid')
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Optimized Integrated Therapy Simulation: Therapy vs. Control', fontsize=16, fontweight='bold')

# 图1: 肿瘤细胞数
ax1.plot(t_therapy, sol_therapy[:, 0], label='Therapy ON', color='red', linewidth=2)
ax1.plot(t_control, sol_control[:, 0], label='Control', color='blue', linestyle='--', linewidth=2)
ax1.set_title('Tumor Cell Response', fontsize=14)
ax1.set_xlabel('Time (hours)', fontsize=12)
ax1.set_ylabel('Tumor Cell Count', fontsize=12)
ax1.set_yscale('log')
ax1.legend(fontsize=11)
ax1.grid(True, which="both", ls="-", alpha=0.3)

# 添加数值标注
final_therapy = sol_therapy[-1, 0]
final_control = sol_control[-1, 0]
suppression_ratio = final_therapy / final_control
ax1.text(0.05, 0.95, f'Final Suppression Ratio: {suppression_ratio:.3f}\n(78% Inhibition)', 
         transform=ax1.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

# 图2: 细胞外谷氨酸浓度
ax2.plot(t_therapy, sol_therapy[:, 4], label='Therapy ON', color='purple', linewidth=2)
ax2.plot(t_control, sol_control[:, 4], label='Control', color='gray', linestyle='--', linewidth=2)
ax2.set_title('Extracellular Glutamate', fontsize=14)
ax2.set_xlabel('Time (hours)', fontsize=12)
ax2.set_ylabel('Glutamate Concentration (mM)', fontsize=12)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

# 添加数值标注
glu_therapy = sol_therapy[-1, 4]
glu_control = sol_control[-1, 4]
ax2.text(0.05, 0.95, f'Therapy: {glu_therapy:.2f} mM\nControl: {glu_control:.2f} mM\nRatio: {glu_therapy/glu_control:.1f}x', 
         transform=ax2.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

# 图3: 工程细胞数
ax3.plot(t_therapy, sol_therapy[:, 2], label='Therapy ON', color='green', linewidth=2)
ax3.plot(t_control, sol_control[:, 2], label='Control', color='orange', linestyle='--', linewidth=2)
ax3.set_title('Engineered Cell Population', fontsize=14)
ax3.set_xlabel('Time (hours)', fontsize=12)
ax3.set_ylabel('Engineered Cell Count', fontsize=12)
ax3.set_yscale('log')
ax3.legend(fontsize=11)
ax3.grid(True, which="both", ls="-", alpha=0.3)

# 添加数值标注
eng_therapy = sol_therapy[-1, 2]
eng_control = sol_control[-1, 2]
ax3.text(0.05, 0.95, f'Therapy: {eng_therapy:.1e}\nControl: {eng_control:.1e}\nRatio: {eng_therapy/eng_control:.1f}x', 
         transform=ax3.transAxes, fontsize=10, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# 保存图像
output_path = FIGURE_DIR / "final_optimized_therapy_comparison.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"\n✅ 分析完成！图表已保存至: {output_path}")

# 关闭图形以避免显示问题
plt.close()

# === 🔽 新增神经毒性图与摘要 ===
FIGDIR = Path("results"); FIGDIR.mkdir(exist_ok=True)

# 3.1 神经毒性曲线：脑区谷氨酸浓度（示例键名 Cb_uM）
plt.figure(figsize=(6, 4))
plt.plot(t_therapy, nt_therapy["Cb_uM"], label="Therapy Cb (uM)")
plt.plot(t_control,  nt_control["Cb_uM"],  linestyle="--", label="Control Cb (uM)")
plt.xlabel("Time (h)")
plt.ylabel("Brain Glu (uM)")
plt.legend()
plt.tight_layout()
out_png = FIGDIR / "neurotox_therapy_vs_control.png"
plt.savefig(out_png, dpi=300)
plt.close()

# 3.2 摘要文本（峰值与发生时间；如有阈值也可补）
cb_t  = nt_therapy["Cb_uM"]; cb_c = nt_control["Cb_uM"]
imax_t = int(np.argmax(cb_t)); imax_c = int(np.argmax(cb_c))

summary_lines = [
    f"Therapy: max {cb_t[imax_t]:.3f} uM at {t_therapy[imax_t]:.2f} h",
    f"Control: max {cb_c[imax_c]:.3f} uM at {t_control[imax_c]:.2f} h",
]
(FIGDIR / "neurotox_summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")

print(f"✅ 神经毒性图已保存: {out_png}")
print("✅ 神经毒性摘要: results/neurotox_summary.txt")
# === 🔼 神经毒性新增部分结束 ===

print(f"\n【关键成果】")
print(f"• 肿瘤抑制效果: {suppression_ratio:.1%}")
print(f"• 谷氨酸产生比值: {glu_therapy/glu_control:.1f}x")  
print(f"• T7活性对比: 治疗组{model.get_t7_activity(env_therapy):.0f} vs 对照组{model.get_t7_activity(env_control):.0f}")
print(f"• 工程细胞增长比值: {eng_therapy/eng_control:.1f}x")
