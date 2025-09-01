
# iGEM AND 门 -> 谷氨酸 -> 铁死亡建模工具包（改良版）

本工具包提供从启动子拟合、AND 门建模、TX-TL 到组织扩散与全身血浆评估的一体化建模流程。默认导出图表（PNG）与数据表（CSV）。

## 改良特性

- 数据验证：检查数据完整性与参数范围
- 拟合改进：约束优化，多指标评估（R2、AIC）
- 可视化增强：统一风格的高分辨率图表（默认导出）
- 异常处理：更清晰的错误信息与提示
- 模块化设计：脚本解耦，易于扩展
- 敏感性分析：关键参数一键扫描
- 批量运行：`run_all.py` 串联主要流程
- 新增 04 扩散-PDE（球对称）
  - 兼容旧版 Dirichlet 边界（保留原文件名与行为）
  - 新增 Robin 边界（血管壁交换）
  - 新增 一室血浆 PK（将边界通量耦合到血浆浓度）

## 系统要求

- Python >= 3.8
- 建议使用 conda 虚拟环境

## 快速开始

1) 安装依赖
```bash
pip install -r requirements.txt
```

2) 检查依赖（可选）
```bash
python run_all.py --check-deps
```

3) 运行完整流程
```bash
python run_all.py
```

> 所有脚本默认会导出图表（PNG）与数据表（CSV）。

## 目录结构

```
iGEM-AND-GLU-Toolkit/
├── run_all.py                   # 主运行脚本
├── requirements.txt             # 依赖
├── README.md                    # 说明文档（本文件）
├── 01_promoter_fit.py           # 启动子传递函数拟合
├── 02_splitT7_AND_model.py      # 分裂 T7 AND 门建模
├── 03_tx_tl_to_glu.py           # TX-TL -> 谷氨酸生产
├── 04_diffusion_pde.py          # 组织扩散 PDE（新增 Robin/PK）
├── 05_cobrapy_dfba_template.py  # 动态 FBA（模板）
├── data/
│   ├── pPept_O2_curve.csv
│   ├── pLR_T_curve.csv
│   └── splitT7_scan.csv
└── params/                      # 自动生成
    ├── promoters.json
    └── splitT7.json
```

## 模型详解

### 1) 启动子传递函数拟合（01_promoter_fit.py）
功能：使用 Hill 方程拟合启动子响应  
特性：激活/抑制自动选择，参数约束，R2/AIC，可视化，数据校验  
输出：`params/promoters.json`，`promoter_fits.png`

### 2) 分裂 T7 AND 门建模（02_splitT7_AND_model.py）
功能：两输入整合，建模 AND 门输出  
特性：scipy 优化，3D 响应面，可选敏感性分析  
输出：`params/splitT7.json`，`splitT7_fit.png`

### 3) TX-TL -> 谷氨酸生产（03_tx_tl_to_glu.py）
功能：从 T7 活性到谷氨酸分泌的 ODE  
特性：面向对象，odeint 求解，多条件比较，敏感性分析  
输出：`tx_tl_glu_dynamics.png`，`glu_comparison.png`，`parameter_sensitivity.png`

### 4) 组织扩散 PDE（04_diffusion_pde.py）

作用：在球对称几何中模拟谷氨酸在组织内的扩散-反应，并评估向血液侧的交换及血浆浓度上升的可能性。

一次运行自动输出三套情景（脚本始终出图）：
1. legacy（兼容旧版）：`spheroid_glu_profile.csv/.png`
2. robin：`spheroid_glu_profile_robin.csv/.png`，`boundary_flux_into_blood.csv/.png`
3. robin_pk：`spheroid_glu_profile_robin_pk.csv/.png`，`boundary_flux_into_blood_robin_pk.csv/.png`，`plasma_concentration.csv/.png`

兼容性：legacy 情景沿用旧文件名，保证不影响既有下游流程；Robin/PK 的结果为新增文件。

#### 04 方程（ASCII 版本，无数学扩展）

```
# In-tissue PDE (spherical symmetry)
dC/dt = D*(d2C/dr2 + 2/r * dC/dr) - kdeg*C + q_glu * rho_bac(r)

# Boundary conditions
Legacy (Dirichlet):  C(R) = 0
Robin (exchange):    -D * (dC/dr at r=R) = h * ( C(R) - C_plasma )

# Plasma one-compartment PK
dC_plasma/dt = Rin(t)/Vd - (ln(2)/t_half) * C_plasma

# Biot number (diagnostic)
Bi = h * R / D
```

#### 04 常用命令

```bash
# 基线（默认参数，三情景一起生成）
python 04_diffusion_pde.py

# 用“24 h 上清 37.5 mM”来标定源强（保守上限），并指定 Robin/PK 参数
python 04_diffusion_pde.py --calibrate-total-mM 37.5 --calibrate-duration-h 24 \
  --h 0.001 --Vd-L 5 --t-half-s 1800 --plasma-baseline-mM 0.08

# 简单敏感性（扫描 h）
for H in 0.0001 0.001 0.01; do python 04_diffusion_pde.py --h $H; done
```

#### 04 关键参数与默认值（表格为纯 Markdown）

| 参数 | 单位 | 默认 | 说明 |
|---|---:|---:|---|
| --D | mm^2/s | 5e-6 | 有效扩散系数（介质/组织；建议做 2~5x 灵敏度） |
| --kdeg | 1/s | 0.0 | 一阶摄取/降解；未知时置 0 做上限 |
| --R | mm | 0.5 | 球半径 |
| --Nr | - | 200 | 径向网格数 |
| --T-h | h | 6.0 | 总时长 |
| --q-glu | mM/s/密度 | 1e-7 | 源强；或用“标定”自动给出 |
| --core-frac | - | 0.4 | 核心半径占比（rho_bac(r) 形状） |
| --core-density | - | 1.0 | 核心密度系数 |
| --shell-density | - | 0.2 | 外壳密度系数 |
| --h | mm/s | 1e-3 | Robin 质传系数（交换效率） |
| --Vd-L | L | 5.0 | 分布容积（按物种/体重） |
| --t-half-s | s | 1800 | 血浆半衰期 |
| --plasma-baseline-mM | mM | 0.0 | 血浆基线（可设 0.05~0.1 mM） |
| --calibrate-total-mM / --calibrate-duration-h | - | - | 用“总体浓度增长”标定 q_glu（例：24 h -> 37.5 mM） |

> 说明：需要实验/文献支撑的参数在代码中均以注释 “TODO(需实验数据)” 标注；可通过命令行覆盖默认值。

## 使用示例

```bash
# 启动子拟合
python 01_promoter_fit.py

# AND 门
python 02_splitT7_AND_model.py

# TX-TL -> 谷氨酸
python 03_tx_tl_to_glu.py

# 组织扩散（自动三情景 + 始终出图）
python 04_diffusion_pde.py
```

代码中自定义（示例）：
```python
from tx_tl_to_glu import TxTlGluModel
model = TxTlGluModel(k_tx=2.0, kcat_gdh=15.0, k_sec=0.1)
t, m, E, C = model.simulate(O2_percent=0.5, Temp_C=42.0)
```

## 故障排除

1) 依赖缺失  
```
pip install -r requirements.txt
```

2) 数据文件缺失  
- 确认 data/ 中 CSV 是否齐全，列名正确

3) 拟合失败  
- 检查数据点数量与范围；去除异常点后重试

4) 图像显示问题（服务器环境）  
```
export MPLBACKEND=Agg
```

## 扩展

- 添加新启动子：在 `01_promoter_fit.py` 中加载新 CSV，加入拟合与可视化
- 修改 TX-TL 动力学参数：在 `TxTlGluModel` 构造函数中传入自定义参数
- PDE 的空间分布：修改 `rho_bac_profile(...)` 或其三个标量参数（核心/外壳）

## 数据格式示例

pPept_O2_curve.csv
```
O2_percent,reporter_MEFL,replicate_id
0.0,1000,1
0.5,800,1
1.0,600,1
...
```

splitT7_scan.csv
```
A_u,B_u,T7_reporter_MEFL
100,50,200
150,75,350
...
```

## License

Apache License 2.0。详见 LICENSE 文件。

## 贡献

欢迎提交 Issue 和 Pull Request。

## 支持

遇到问题请先查看：
1) 本 README 的“故障排除”和“04 参数表”  
2) 代码中的注释与文档字符串  
3) 项目仓库的 Issue 区

--- 

备注：本 README 采用纯 ASCII 与基础 Markdown，不依赖 LaTeX/MathJax。

