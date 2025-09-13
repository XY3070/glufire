# 脚本机制与逻辑文档

## 1. 主要脚本说明

### 1.1 main.py
- **功能**: 程序入口，调用其他模块
- **输入**: 无
- **输出**: 打印状态信息
- **参数**: 无

### 1.2 generate_final_analysis.py
- **功能**: 生成最终优化的分析图表和模拟结果，包括肿瘤细胞响应、细胞外谷氨酸、工程细胞数量以及神经毒性分析。
- **输入**: 
  - `env_therapy`: 治疗组环境条件 (O2_percent, Temp_C)
  - `env_control`: 对照组环境条件 (O2_percent, Temp_C)
  - `t_end`: 模拟时长 (默认50小时)
  - `dt`: 时间步长 (默认1.0小时)
- **输出**: 
  - `final_optimized_therapy_comparison.png`: 整合治疗模拟对比图 (肿瘤细胞、细胞外谷氨酸、工程细胞)
  - `neurotox_therapy_vs_control.png`: 神经毒性对比图 (脑区谷氨酸浓度)
  - `neurotox_summary.txt`: 神经毒性摘要 (峰值及发生时间)
- **参数**: 无需命令行参数，所有参数在脚本内部定义。

### 1.3 run_analysis.py
- **功能**: 运行并可视化AND门响应面、谷氨酸生产动力学和整合治疗模型的效果。
- **输入**: 无需命令行参数，所有参数在脚本内部定义。
- **输出**: 
  - `and_gate_response.png`: AND门响应面图
  - `glutamate_production.png`: 谷氨酸生产和分泌动力学图
  - `integrated_therapy_comparison.png`: 整合治疗模拟对比图
- **参数**: 无需命令行参数，所有参数在脚本内部定义。

### 1.4 assess_neurotoxicity.py
- **功能**: 模拟谷氨酸的药代动力学 (PK) 并评估神经毒性风险。
- **输入**: 无需命令行参数，所有参数在脚本内部定义。
- **输出**: 
  - `plasma_glu_neurotoxicity.png`: 血浆谷氨酸浓度与神经毒性阈值对比图。
  - 神经毒性风险报告 (打印到控制台)。
- **参数**:
  - `sim_hours`: 模拟总时长 (小时，默认 48.0)。
  - `dt_h`: 模拟时间步长 (小时，默认 0.1)。
  - `baseline_uM`: 血浆/组织/脑脊液的谷氨酸基线浓度 (μM，默认 50.0)。
  - `intensity`: 谷氨酸分泌强度 ('mild' 或 'high'，默认 'mild')。
  - `peak_umol_per_h`: 谷氨酸峰值通量 (μmol/h，根据 `intensity` 自动设置)。
  - `t_on`: 谷氨酸通量开始时间 (小时，默认 8.0)。
  - `ramp_hours`: 谷氨酸通量上升和下降的持续时间 (小时，默认 2.0)。
  - `t_plateau`: 谷氨酸通量达到平台的时间 (小时，根据 `t_on` 和 `ramp_hours` 自动计算)。
  - `t_off`: 谷氨酸通量开始下降的时间 (小时，默认 34.0)。

## 2. 模型模块说明

### 2.1 and_gate.py
- **模型**: 简单AND门模型 (已外部化参数)
- **公式**: Hill函数
  - **激活型**: \( \text{Output} = \text{leaky} + \text{beta} \times \frac{x^n}{K^n + x^n} \)
  - **抑制型**: \( \text{Output} = \text{leaky} + \text{beta} \times \frac{K^n}{K^n + x^n} \)
  其中 \( x \) 是输入信号 (O2% 或温度)。
- **输入**: O2% (氧气浓度), 温度 (摄氏度)
- **输出**: T7活性 (任意单位)
- **参数**:
  - **所有参数均通过 `ConfigManager` 进行管理和加载。** 默认值在 `config_manager.py` 中定义，并可在需要时通过 JSON 配置文件覆盖。
  - **关键参数示例 (完整列表请参考 `config_manager.py` 中的 `and_gate` 部分)**:
    - `promoter_params`:
      - `pPept` (低氧诱导启动子):
        - `type`: 'rep' (抑制型)
        - `beta`: 最大表达量
        - `K`: 半数激活/抑制浓度
        - `n`: Hill系数
        - `leaky`: 基础泄漏表达
      - `pLR` (高温激活启动子):
        - `type`: 'act' (激活型)
        - `beta`: 最大表达量
        - `K`: 半数激活/抑制温度
        - `n`: Hill系数
        - `leaky`: 基础泄漏表达
    - `splitT7_params`:
      - `alpha`: T7组装效率
      - `Kd`: T7组装的半饱和常数
      - `leaky`: T7基础活性
    - `k_assembly`: T7组装速率常数
    - `k_disassembly`: T7解体速率常数
    - `k_deg`: T7降解速率常数
  - **辅助方法**:
    - `quick_diagnose()`: 用于快速诊断在典型条件下的表达与T7输出，帮助定位参数问题。

### 2.2 glu_metabolism.py
- **模型**: 谷氨酸生产和分泌的ODE模型 (已外部化参数)
- **公式**: 常微分方程组 (dydt 方法)
  - \( \frac{dGlc_{ext}}{dt} = -V_{max\_glc} \times \frac{Glc_{ext}}{K_{m\_glc} + Glc_{ext}} \times X \)
  - \( \frac{dNH4_{ext}}{dt} = -V_{max\_glc} \times \frac{Glc_{ext}}{K_{m\_glc} + Glc_{ext}} \times X \times Y_{X/Glc} \times Y_{NH4/X} \)
  - \( \frac{dICIT}{dt} = (V_{max\_base\_ICD} \times fold_{ICD} \times \frac{Glc_{ext}}{K_{m\_ICD} + Glc_{ext}}) - (V_{max\_base\_GDH} \times fold_{GDH} \times \frac{ICIT}{K_{m\_AKG} + ICIT}) - (k_{PPP} \times ICIT) \)
  - \( \frac{dAKG}{dt} = (V_{max\_base\_GDH} \times fold_{GDH} \times \frac{ICIT}{K_{m\_AKG} + ICIT}) - (k_{sec\_base} \times AKG) \)
  - \( \frac{dGlu_{in}}{dt} = (k_{sec\_base} \times AKG) - (k_{export\_base} \times Glu_{in}) - (k_{maintenance} \times Glu_{in}) \)
  - \( \frac{dNADPH}{dt} = (y_{ICD\_NADPH} \times V_{ICD}) - (\lambda_{NADPH} \times (NADPH - NADPH_{set})) \)
  - \( \frac{dX}{dt} = \mu_{max} \times \frac{Glc_{ext}}{K_{m\_glc} + Glc_{ext}} \times X \)
  - \( \frac{dGlu_{ext}}{dt} = (k_{export\_base} \times Glu_{in}) - (extracellular\_clearance\_rate \times Glu_{ext}) \)
  - \( \frac{dfold_{ICD}}{dt} = \frac{1}{\tau_{enzyme}} \times (\frac{T7_{activity}^{n_{hill}}}{K_{T7}^{n_{hill}} + T7_{activity}^{n_{hill}}} \times fold_{ICD\_max} - fold_{ICD}) \)
  - \( \frac{dfold_{GDH}}{dt} = \frac{1}{\tau_{enzyme}} \times (\frac{T7_{activity}^{n_{hill}}}{K_{T7}^{n_{hill}} + T7_{activity}^{n_{hill}}} \times fold_{GDH\_max} - fold_{GDH}) \)
- **输入**: T7聚合酶活性 (AU), 初始状态变量
- **输出**: 细胞内谷氨酸浓度 (Glu_in), 细胞外谷氨酸浓度 (Glu_ext), 其他状态变量
- **参数**:
  - **所有参数均通过 `ConfigManager` 进行管理和加载。** 默认值在 `config_manager.py` 中定义，并可在需要时通过 JSON 配置文件覆盖。
  - **关键参数示例 (完整列表请参考 `config_manager.py` 中的 `glu_metabolism_en` 部分)**:
    - `V_max_glc`: 最大葡萄糖摄取速率
    - `K_m_glc`: 葡萄糖米氏常数
    - `f_TCA`: TCA循环通量比例
    - `V_max_base_ICD`: 异柠檬酸脱氢酶 (ICD) 基础最大活性
    - `K_m_ICD`: ICD米氏常数
    - `V_max_base_GDH`: 谷氨酸脱氢酶 (GDH) 基础最大活性
    - `K_m_AKG`: α-酮戊二酸米氏常数
    - `K_m_NH4`: 铵离子米氏常数
    - `K_m_NADPH`: NADPH米氏常数
    - `k_PPP`: 磷酸戊糖途径 (PPP) 速率常数
    - `y_ICD_NADPH`: ICD生产NADPH的产率
    - `lambda_NADPH`: NADPH稳态调节强度
    - `NADPH_set`: NADPH设定点
    - `k_sec_base`: 基础分泌速率常数
    - `k_maintenance`: 维持代谢速率常数
    - `mu_max`: 最大比生长速率
    - `K_T7`: T7活性半饱和常数
    - `n_hill`: Hill系数
    - `tau_enzyme`: 酶表达时间常数
    - `fold_ICD_max`: ICD最大过表达倍数
    - `fold_GDH_max`: GDH最大过表达倍数
    - `homeostasis_strength`: 细胞内谷氨酸稳态调节强度
    - `accum_threshold`: 谷氨酸积累阈值
    - `export_accum_suppression`: 积累抑制下的谷氨酸输出比例
    - `postshock_export_boost`: 热激后谷氨酸输出增强因子
    - `extracellular_clearance_rate`: 细胞外谷氨酸清除率
    - `export_decay_rate`: 谷氨酸输出衰减率
    - `Glu_target`: 谷氨酸目标浓度

### 2.3 diffusion_pk.py
- **模型**: 扩散与药代动力学模块 (已外部化参数)
  1. **MultiCompartmentPK**: 多室药代动力学模型，描述药物在血液、肝脏、肿瘤等室的分布。
  2. **TumorDiffusion**: 肿瘤微环境中的反应-扩散模型 (1D)，描述药物在肿瘤中的空间分布。
- **公式**:
  1. **MultiCompartmentPK**: 常微分方程组 (ODE)
     - \( \frac{dC_{blood}}{dt} = \frac{Infusion}{V_{blood}} + \sum_{i \neq blood} \frac{q_{i \to blood} C_i - q_{blood \to i} C_{blood}}{V_{blood}} \)
     - \( \frac{dC_{liver}}{dt} = \frac{q_{blood \to liver} C_{blood} - q_{liver \to blood} C_{liver}}{V_{liver}} - k_{elim,liver} C_{liver} \)
     - \( \frac{dC_{tumor}}{dt} = \frac{q_{blood \to tumor} C_{blood} - q_{tumor \to blood} C_{tumor}}{V_{tumor}} - k_{uptake,tumor} C_{tumor} \)
     - \( \frac{dC_{other}}{dt} = \frac{q_{blood \to other} C_{blood} - q_{other \to blood} C_{other}}{V_{other}} \)
  2. **TumorDiffusion**: 偏微分方程 (PDE)
     - \( \frac{\partial C}{\partial t} = D \frac{\partial^2 C}{\partial x^2} - k_{uptake} C \)
- **输入**:
  1. **MultiCompartmentPK**: 药物输注速率 (infusion_rate_func), 模拟时长 (t_end), 时间步长 (dt)。
  2. **TumorDiffusion**: 肿瘤边界药物浓度 (C_boundary), 模拟时长 (t_end), 时间步长 (dt)。
- **输出**:
  1. **MultiCompartmentPK**: 各室药物浓度随时间的变化。
  2. **TumorDiffusion**: 肿瘤内药物浓度随空间和时间的变化。
- **参数**:
  - **所有参数均通过 `ConfigManager` 进行管理和加载。** 默认值在 `config_manager.py` 中定义，并可在需要时通过 JSON 配置文件覆盖。
  - **关键参数示例 (完整列表请参考 `config_manager.py` 中的 `diffusion_pk` 部分)**:
    - **MultiCompartmentPK**:
      - `V_blood`, `V_liver`, `V_tumor`, `V_other`: 各室体积 (L)
      - `q_blood_liver`, `q_blood_tumor`, `q_blood_other`: 血液与其他室之间的流速 (L/hr)
      - `k_elim_liver`: 肝脏清除率 (1/hr)
      - `k_uptake_tumor`: 肿瘤摄取率 (1/hr)
    - **TumorDiffusion**:
      - `D`: 扩散系数 (cm²/hr)
      - `k_uptake`: 肿瘤细胞摄取率 (1/hr)
      - `L`: 肿瘤直径 (cm)
      - `Nx`: 空间网格点数

### 2.4 integrated_model.py
- **模型**: 整合的端到端治疗模型
- **功能**: 将AND门、谷氨酸代谢和肿瘤生长/死亡模型整合在一起，模拟在特定环境条件下（温度和氧气），整个系统的治疗效果。模型基于ODE，描述了肿瘤细胞、死亡细胞和细胞外谷氨酸浓度的动态变化。
- **公式**: ODE系统 (dydt 方法)
  - \( \frac{dN_{tumor}}{dt} = r \times N_{tumor} \times (1 - \frac{N_{tumor} + N_{eng}}{K_{tumor}}) - k_{ferroptosis\_max} \times \frac{Glu_{extra}^{n_{glu}}}{K_{glu}^{n_{glu}} + Glu_{extra}^{n_{glu}}} \times N_{tumor} \)
  - \( \frac{dD_{tumor}}{dt} = k_{ferroptosis\_max} \times \frac{Glu_{extra}^{n_{glu}}}{K_{glu}^{n_{glu}} + Glu_{extra}^{n_{glu}}} \times N_{tumor} \)
  - \( \frac{dN_{eng}}{dt} = r_{eng} \times N_{eng} \times (1 - \frac{N_{tumor} + N_{eng}}{K_{tumor}}) \times \frac{T7_{activity}^{n_{hill}}}{K_{t7}^{n_{hill}} + T7_{activity}^{n_{hill}}} - k_{dilution} \times N_{eng} \)
  - 谷氨酸代谢相关的ODE (dGlu_intra/dt, dGlu_extra/dt, dIcd/dt, dgdhA/dt) 详见 `glu_metabolism.py` 部分。
- **输入**: 环境条件 (O2_percent, Temp_C)
- **输出**: 存活肿瘤细胞数 (N_tumor), 死亡肿瘤细胞数 (D_tumor), 存活工程细胞数 (N_eng), 细胞内谷氨酸浓度 (Glu_intra), 细胞外谷氨酸浓度 (Glu_extra), Icd表达水平, gdhA表达水平
- **参数**:
  - **固定参数 (通常无需调整)**:
    - `r`: 肿瘤生长速率 (默认 0.01)
    - `K_tumor`: 肿瘤承载能力 (细胞数, 默认 1e9)
    - `r_eng`: 工程细胞生长速率 (默认 0.2)
    - `K_eng`: 工程细胞承载能力 (默认 5e8)
    - `V_cell`: 单个细胞体积 (L/cell, 默认 2e-12)
    - `V_tumor_ext`: 肿瘤间质液体积 (L, 默认 0.01)
    - `k_ferroptosis_max`: 最大铁死亡速率 (默认 15.0)。影响谷氨酸诱导铁死亡的强度。
    - `K_glu`: 谷氨酸铁死亡阈值 (默认 30.0 mM)。当细胞外谷氨酸浓度达到此值时，铁死亡速率达到一半。可根据细胞系对谷氨酸缺乏的敏感性进行调整。
    - `n_glu`: 铁死亡Hill系数 (默认 5.0)。影响铁死亡响应谷氨酸浓度的陡峭程度。值越大，响应越接近开关效应。
    - 内部包含 `and_gate.py` 和 `glu_metabolism.py` 的参数，具体参数详见各自模块说明.

### 2.5 pk_toxicity.py
- **功能**: 模拟三室药代动力学 (PK) 并评估神经毒性。
- **输入**:
  - `t_grid_h`: 时间点数组 (小时)
  - `S_t_umol_per_h`: 肿瘤分泌通量 (μmol/h)
  - `params`: PKParams 对象，包含PK参数
  - `Cb0_uM`, `Ct0_uM`, `Cn0_uM`: 各室初始谷氨酸浓度 (μM)
  - `thr`: ToxicityThresholds 对象，包含毒性阈值
- **输出**:
  - `Cb_uM`, `Ct_uM`, `Cn_uM`: 血浆、肿瘤、正常组织谷氨酸浓度随时间的变化 (μM)
  - `tox_report`: 神经毒性评估报告 (字典，包含峰值、超阈值时间等)
- **参数**:
  - **PKParams (药代动力学参数)**:
    - `Vb`: 血浆体积 (L, 默认 0.002)
    - `Vt`: 肿瘤体积 (L, 默认 0.0005)
    - `Vn`: 正常组织等效体积 (L, 默认 0.02)
    - `k_bt`: 血<->肿瘤交换速率 (h^-1, 默认 1.0)
    - `k_bn`: 血<->正常组织交换速率 (h^-1, 默认 0.5)
    - `k_b_clr`: 血浆清除率 (h^-1, 默认 0.5)
    - `k_t_clr`: 肿瘤清除率 (h^-1, 默认 0.2)
    - `k_n_clr`: 正常组织清除率 (h^-1, 默认 0.2)
  - **ToxicityThresholds (毒性阈值参数)**:
    - `caution_um`: 提示阈值 (μM, 默认 100.0)
    - `danger_um`: 危险阈值 (μM, 默认 1000.0)



## 4. 数据文件说明

- `promoters.json`: 启动子参数
- `splitT7.json`: T7参数
- 各种CSV数据文件