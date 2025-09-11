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
- **模型**: 简单AND门模型
- **公式**: Hill函数
  - **激活型**: \( \text{Output} = \text{leaky} + \text{beta} \times \frac{x^n}{K^n + x^n} \)
  - **抑制型**: \( \text{Output} = \text{leaky} + \text{beta} \times \frac{K^n}{K^n + x^n} \)
  其中 \( x \) 是输入信号 (O2% 或温度)。
- **输入**: O2% (氧气浓度), 温度 (摄氏度)
- **输出**: T7活性 (任意单位)
- **参数**:
  - **固定参数 (通常无需调整)**:
    - `promoter_params` (从 `params/promoters.json` 加载):
      - `pPept` (低氧诱导启动子):
        - `type`: 'rep' (抑制型)
        - `beta`: 最大表达量 (默认1200.0)
        - `K`: 半数激活/抑制浓度 (默认5.0 %O2)
        - `n`: Hill系数 (默认2.0)
        - `leaky`: 基础泄漏表达 (默认50.0)
      - `pLR` (高温激活启动子):
        - `type`: 'act' (激活型)
        - `beta`: 最大表达量 (默认1500.0)
        - `K`: 半数激活/抑制温度 (默认40.0 °C)
        - `n`: Hill系数 (默认8.0)
        - `leaky`: 基础泄漏表达 (默认80.0)
    - `splitT7_params` (从 `params/splitT7.json` 加载):
      - `alpha`: T7组装效率 (默认1.0)
      - `Kd`: T7组装的半饱和常数 (默认2.0e5)
      - `leaky`: T7基础活性 (默认200.0)
  - **辅助方法**:
    - `quick_diagnose()`: 用于快速诊断在典型条件下的表达与T7输出，帮助定位参数问题。

### 2.2 glu_metabolism.py
- **模型**: 谷氨酸生产和分泌的ODE模型
- **公式**: 常微分方程组 (dydt 方法)
  - \( \frac{dIcd}{dt} = k_{syn\_icd} \times \frac{T7_{activity}^{n_{hill}}}{K_{t7}^{n_{hill}} + T7_{activity}^{n_{hill}}} - k_{deg\_icd} \times Icd \)
  - \( \frac{dgdhA}{dt} = k_{syn\_gdhA} \times \frac{T7_{activity}^{n_{hill}}}{K_{t7}^{n_{hill}} + T7_{activity}^{n_{hill}}} - k_{deg\_gdhA} \times gdhA \)
  - \( v_{prod} = Vmax_{gdhA} \times gdhA \)
  - \( v_{export} = k_{export\_max} \times \frac{Glu_{intra}}{K_{export} + Glu_{intra}} \)
  - \( \frac{dGlu_{intra}}{dt} = v_{prod} - v_{export} - k_{dilution} \times Glu_{intra} \)
  - \( \frac{dGlu_{extra}}{dt} = v_{export} \times V_{ratio} - k_{dilution} \times Glu_{extra} \)
- **输入**: T7聚合酶活性 (AU)
- **输出**: 细胞内谷氨酸浓度 (Glu_intra), 细胞外谷氨酸浓度 (Glu_extra), Icd表达水平, gdhA表达水平
- **参数**:
  - **固定参数 (通常无需调整)**:
    - `k_prod_max`: 最大谷氨酸生产速率 (mM/hr, 默认 50.0)
    - `K_t7`: T7活性达到半最大生产速率时的值 (AU, 默认 500.0)
    - `k_export_max`: 最大谷氨酸分泌速率 (mM/hr, 默认 100.0)
    - `K_export`: 细胞内谷氨酸浓度达到半最大分泌速率时的值 (mM, 默认 10.0)
    - `k_dilution`: 稀释/降解速率 (1/hr, 默认 0.1)
    - `V_intra_over_V_extra` (`V_ratio`): 细胞内总体积与细胞外总体积的比率 (默认 0.01)
    - `k_syn_icd`: Icd合成速率 (1/hr, 默认 2.0)
    - `k_syn_gdhA`: gdhA合成速率 (1/hr, 默认 2.0)
    - `k_deg_icd`: Icd降解速率 (1/hr, 默认 0.2)
    - `k_deg_gdhA`: gdhA降解速率 (1/hr, 默认 0.2)
    - `Vmax_icd`: Icd最大催化速率 (mM/hr, 默认 100.0)
    - `K_icd`: Icd底物常数 (mM, 默认 5.0)
    - `Vmax_gdhA`: gdhA最大催化速率 (mM/hr, 默认 100.0)
    - `K_gdhA`: gdhA底物常数 (mM, 默认 5.0)
    - `n_hill`: Hill系数 (默认 4.0)

### 2.3 diffusion_pk.py
- **功能**: 提供两种模型来描述治疗药物（如谷氨酸）在体内的分布和扩散：
  1. MultiCompartmentPK: 一个多室药代动力学(PK)模型，用于模拟药物在全身各主要器官（血液、肝脏、肿瘤等）的分布。基于ODE。
  2. TumorDiffusion: 一个反应-扩散模型，用于模拟药物在肿瘤微环境(TME)中的空间分布。基于偏微分方程(PDE)。
- **输入**:
  - MultiCompartmentPK: 药物输注速率 (`infusion_rate_func`)
  - TumorDiffusion: 肿瘤边界的药物浓度 (`C_boundary`)
- **输出**:
  - MultiCompartmentPK: 各室药物浓度随时间的变化
  - TumorDiffusion: 肿瘤内药物浓度随空间和时间的变化
- **参数**:
  - **MultiCompartmentPK (多室药代动力学模型)**:
    - `V_blood`: 血液体积 (L, 默认 5.0)
    - `V_liver`: 肝脏体积 (L, 默认 1.5)
    - `V_tumor`: 肿瘤体积 (L, 默认 0.5)
    - `V_other`: 其他组织体积 (L, 默认 60.0)
    - `q_blood_liver`: 血液到肝脏的流速 (L/hr, 默认 90.0)
    - `q_blood_tumor`: 血液到肿瘤的流速 (L/hr, 默认 10.0)
    - `q_blood_other`: 血液到其他组织的流速 (L/hr, 默认 200.0)
    - `k_elim_liver`: 药物在肝脏的清除率 (1/hr, 默认 0.5)
    - `k_uptake_tumor`: 肿瘤摄取率 (1/hr, 默认 0.2)
  - **TumorDiffusion (肿瘤扩散模型)**:
    - `D`: 扩散系数 (cm²/hr, 默认 0.01)
    - `k_uptake`: 肿瘤细胞摄取率 (1/hr, 默认 0.2)
    - `L`: 肿瘤直径 (cm, 默认 1.0)
    - `Nx`: 空间网格点数 (默认 50)

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