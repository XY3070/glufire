# iGEM AND门 → 谷氨酸 → 铁死亡建模工具包 (改良版)

这个工具包为您提供完整的建模框架，用于分析环境响应的AND逻辑门控制的谷氨酸生产系统。

## 🌟 改良特性

- **增强的数据验证**: 自动检查数据质量和参数合理性
- **改进的拟合算法**: 使用约束优化和多种评估指标
- **可视化增强**: 生成高质量的分析图表
- **错误处理**: 完善的异常处理和用户友好的错误信息
- **模块化设计**: 面向对象的模型设计，易于扩展
- **敏感性分析**: 自动进行参数敏感性分析
- **批量运行**: 一键运行所有分析流程

## 📋 系统要求

- Python >= 3.8
- 推荐使用conda环境

## 🚀 快速开始

### 1. 安装依赖
```bash
pip install -r requirements.txt
```

### 2. 检查依赖
```bash
python run_all.py --check-deps
```

### 3. 运行完整分析
```bash
python run_all.py
```

## 📊 文件结构
```
📦 iGEM AND门建模工具包/
├── 🐍 run_all.py                 # 主运行脚本
├── 📋 requirements.txt           # 依赖包列表
├── 📖 README.md                  # 项目说明
├── 🧬 01_promoter_fit.py         # 启动子传递函数拟合
├── 🔬 02_splitT7_AND_model.py    # 分裂T7 AND门建模  
├── ⚗️  03_tx_tl_to_glu.py         # TX-TL谷氨酸生产模拟
├── 🌊 04_diffusion_pde.py        # 扩散PDE模型
├── 🔬 05_cobrapy_dfba_template.py # COBRApy动态FBA
├── 📁 data/                      # 实验数据
│   ├── pPept_O2_curve.csv       # pPept启动子氧气响应数据
│   ├── pLR_T_curve.csv          # pL/pR启动子温度响应数据
│   └── splitT7_scan.csv         # 分裂T7扫描数据
└── 📁 params/                    # 模型参数 (自动生成)
    ├── promoters.json           # 启动子参数
    └── splitT7.json            # 分裂T7参数
```

## 🔬 模型详解

### 1. 启动子传递函数拟合 (`01_promoter_fit.py`)

**功能**: 使用Hill方程拟合启动子响应曲线

**改良特性**:
- 自动模型选择（激活型/抑制型）
- 参数约束优化
- 拟合质量评估（R²、AIC）
- 高质量可视化
- 数据验证和错误处理

**输出**:
- `params/promoters.json`: 拟合参数
- `promoter_fits.png`: 拟合结果图

### 2. 分裂T7 AND门建模 (`02_splitT7_AND_model.py`)

**功能**: 整合两个环境输入，建模分裂T7 RNA聚合酶活性

**改良特性**:
- 基于scipy的优化算法
- 3D响应面可视化
- AND逻辑门验证
- 参数敏感性分析

**输出**:
- `params/splitT7.json`: AND门参数
- `splitT7_fit.png`: 拟合和响应曲线

### 3. TX-TL谷氨酸生产模拟 (`03_tx_tl_to_glu.py`)

**功能**: 从T7活性到谷氨酸生产的完整动力学模拟

**改良特性**:
- 面向对象的模型设计
- 使用scipy.integrate.odeint求解ODE
- 多条件比较分析
- 参数敏感性分析
- 动力学可视化

**输出**:
- `tx_tl_glu_dynamics.png`: 动力学曲线
- `glu_comparison.png`: 条件比较
- `parameter_sensitivity.png`: 敏感性分析

## 📈 使用示例

### 单独运行脚本
```bash
# 拟合启动子响应
python 01_promoter_fit.py

# 建模AND门
python 02_splitT7_AND_model.py  

# 模拟谷氨酸生产
python 03_tx_tl_to_glu.py
```

### 自定义参数
```python
from tx_tl_to_glu import TxTlGluModel

# 创建自定义模型
model = TxTlGluModel(
    k_tx=2.0,      # 提高转录速率
    kcat_gdh=15.0, # 提高酶活性
    k_sec=0.1      # 提高分泌效率
)

# 运行模拟
t, m, E, C = model.simulate(O2_percent=0.5, Temp_C=42.0)
```

## 🔧 故障排除

### 常见问题

1. **依赖包缺失**
   ```bash
   pip install -r requirements.txt
   ```

2. **数据文件缺失**
   - 确保`data/`目录包含所需CSV文件
   - 检查文件格式是否正确

3. **拟合失败**
   - 检查数据质量（是否有足够数据点）
   - 验证数据范围的合理性

4. **图像显示问题**
   ```bash
   # 在服务器环境中
   export MPLBACKEND=Agg
   ```

## 📚 扩展功能

### 添加新的启动子
在`01_promoter_fit.py`中添加新的数据加载和拟合代码：

```python
# 加载新启动子数据
new_promoter_df = pd.read_csv("data/new_promoter_curve.csv")
new_promoter_mean = new_promoter_df.groupby("input_var")["output_var"].mean()

# 拟合
new_params, _ = fit_curve(new_promoter_mean, "input_var", "output_var")
```

### 修改动力学参数
```python
# 在TxTlGluModel中修改默认参数
model = TxTlGluModel(
    k_tx=1.5,      # 转录速率
    k_tl=2.0,      # 翻译速率  
    kcat_gdh=20.0, # 酶催化常数
    # ... 其他参数
)
```

## 📄 数据格式

### pPept_O2_curve.csv
```csv
O2_percent,reporter_MEFL,replicate_id
0.0,1000,1
0.5,800,1
1.0,600,1
...
```

### splitT7_scan.csv  
```csv
A_u,B_u,T7_reporter_MEFL
100,50,200
150,75,350
...
```

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个工具包！

## 📞 支持

如有问题，请查看：
1. 本README的故障排除部分
2. 代码注释和文档字符串
3. 提交Issue到项目仓库

---
*这个工具包专为iGEM团队设计，用于合成生物学系统的建模和分析。*
