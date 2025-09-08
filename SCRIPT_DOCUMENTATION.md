# 脚本机制与逻辑文档

## 1. 主要脚本说明

### 1.1 main.py
- **功能**: 程序入口，调用其他模块
- **输入**: 无
- **输出**: 打印状态信息
- **参数**: 无

### 1.2 generate_final_analysis.py
- **功能**: 生成最终分析图表和模拟结果
- **输入**: 模型参数、模拟条件
- **输出**: 分析图表(png)
- **参数**: 
  - `env_conditions`: 环境条件字典(O2_percent, Temp_C)
  - `t_end`: 模拟时长(默认100.0)
  - `dt`: 时间步长(默认0.5)

### 1.3 run_analysis.py
- **功能**: 运行和可视化分析结果
- **输入**: 模型参数
- **输出**: 可视化结果
- **参数**:
  - `env_conditions`: 环境条件
  - `initial_conditions`: 初始条件

## 2. 模型模块说明

### 2.1 and_gate.py
- **模型**: 简单AND门模型
- **公式**: Hill函数
  - **激活型**: \( 	ext{Output} = 	ext{leaky} + 	ext{beta} 	imes \frac{x^n}{K^n + x^n} \)
  - **抑制型**: \( 	ext{Output} = 	ext{leaky} + 	ext{beta} 	imes \frac{K^n}{K^n + x^n} \)
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

### 2.3 integrated_model.py
- **模型**: 整合的端到端治疗模型
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
    - 内部包含 `and_gate.py` 和 `glu_metabolism.py` 的参数，具体参数详见各自模块说明.
  - **用户可调参数 (可根据实验数据或优化目标调整)**:
    - `k_ferroptosis_max`: 最大铁死亡速率 (默认 15.0)。影响谷氨酸诱导铁死亡的强度。
    - `K_glu`: 谷氨酸铁死亡阈值 (默认 30.0 mM)。当细胞外谷氨酸浓度达到此值时，铁死亡速率达到一半。可根据细胞系对谷氨酸缺乏的敏感性进行调整。
    - `n_glu`: 铁死亡Hill系数 (默认 5.0)。影响铁死亡响应谷氨酸浓度的陡峭程度。值越大，响应越接近开关效应。

## 3. 用户可调参数 (可根据实验数据或优化目标调整)

- `k_ferroptosis_max`: 最大铁死亡速率 (默认 15.0)。影响谷氨酸诱导铁死亡的强度。
- `K_glu`: 谷氨酸铁死亡阈值 (默认 30.0 mM)。当细胞外谷氨酸浓度达到此值时，铁死亡速率达到一半。可根据细胞系对谷氨酸缺乏的敏感性进行调整。
- `n_glu`: 铁死亡Hill系数 (默认 5.0)。影响铁死亡响应谷氨酸浓度的陡峭程度。值越大，响应越接近开关效应。
- `env_conditions`: 环境条件 (O2%, 温度)。用于模拟不同环境对AND门和T7活性的影响。例如，可以通过调整O2%和温度来模拟肿瘤微环境或热疗条件。
- `initial_conditions`: 初始细胞数量。包括初始肿瘤细胞数、死亡细胞数和工程细胞数。用于设置模拟的起始状态，例如，可以模拟不同肿瘤负荷或不同工程细胞接种量的情况。

## 4. 数据文件说明

- `promoters.json`: 启动子参数
- `splitT7.json`: T7参数
- 各种CSV数据文件