# glufire 项目封装开发文档

## 项目结构和功能梳理

**项目名称:** `glufire`

**项目目标:**
`glufire` 是一个用于模拟和分析合成生物学系统的工具包，专注于通过环境响应的 AND 逻辑门控制谷氨酸生产，从而诱导肿瘤细胞铁死亡。

**核心模块:**

1.  **AND-gate (glufire/models/and_gate.py)**
    *   **功能:** 模拟环境响应的 AND 逻辑门。根据低氧 (1%) 和高温 (42°C) 等环境条件，控制 T7 RNA 聚合酶的活性。
    *   **输入:** 氧气浓度、温度。
    *   **输出:** T7 聚合酶活性。
    *   **关键概念:** 基因表达调控、环境传感器。

2.  **Glu Metabolism (glufire/models/glu_metabolism_en.py)**
    *   **功能:** 模拟谷氨酸的生产和分泌过程。T7 聚合酶的活性将直接影响谷氨酸的生成速率。
    *   **输入:** T7 聚合酶活性。
    *   **输出:** 谷氨酸浓度。
    *   **关键概念:** 代谢通路、酶动力学。

3.  **Diffusion (glufire/models/diffusion_pk_neurotoxicity.py)**
    *   **功能:** 模拟谷氨酸在生物体内的扩散和药代动力学 (PK) 过程，并评估其神经毒性。
    *   **输入:** 谷氨酸浓度。
    *   **输出:** 不同组织（如脑部、血浆）的谷氨酸浓度，以及神经毒性评估。
    *   **关键概念:** 药代动力学、多隔室模型、神经毒性阈值。

**辅助文件和目录:**

*   **`main.py`**: 项目入口脚本，可能包含一些测试或集成运行逻辑。
*   **`run_analysis.py`**: 用于运行所有核心模型并生成分析图表的主分析脚本。
*   **`generate_final_analysis.py`**: 生成最终优化分析的脚本，包含 therapy vs control 对比和脑部谷氨酸浓度曲线。
*   **`config_manager.py`**: 可能用于管理项目配置和参数。
*   **`data/`**: 存放实验数据，如 `pLR_T_curve.csv` 等。
*   **`params/`**: 存放模型参数，如 `promoters.json` 和 `splitT7.json`。
*   **`results/`**: 存放分析结果和生成的图表。
*   **`requirements.txt` / `uv.lock` / `pyproject.toml`**: 项目依赖管理文件。

## 项目封装提纲

为了将 `glufire` 封装成一个可安装的 Python 包，并提供 CLI 接口，我们将采取以下步骤：

**1. 更新 `uv` 配置 (`pyproject.toml`)**

*   **目的:** 定义包的元数据、依赖项和 CLI 入口点。
*   **修改内容:**
    *   在 `[project]` 部分添加 `name = "glufire"`, `version`, `description`, `authors` 等信息。
    *   在 `[project.dependencies]` 中列出所有项目所需的第三方库。
    *   在 `[project.scripts]` 或 `[project.entry-points."console_scripts"]` 中定义 CLI 入口点。例如：
        ```toml
        [project.scripts]
        glufire = "glufire.cli:main"
        ```
        或者更具体的子命令：
        ```toml
        [project.entry-points."console_scripts"]
        glufire = "glufire.cli:main"
        glufire-and-gate = "glufire.cli:and_gate_command"
        glufire-glu-metabolism = "glufire.cli:glu_metabolism_command"
        glufire-diffusion = "glufire.cli:diffusion_command"
        ```
        考虑到你希望 `glufire and-gate <option>` 这样的形式，我们可能需要一个主 `glufire` 命令，然后通过 `argparse` 或 `click` 等库实现子命令。

**2. 组织项目结构**

*   **目的:** 遵循 Python 包的最佳实践，使模块易于导入和管理。
*   **修改内容:**
    *   创建一个顶层包目录 `glufire/`。
    *   将核心模块 (`and_gate.py`, `glu_metabolism_en.py`, `diffusion_pk_neurotoxicity.py`) 移动到 `glufire/models/` 目录下。
    *   创建 `glufire/__init__.py` 文件，使其成为一个 Python 包。
    *   创建 `glufire/cli.py` 文件，用于处理命令行接口逻辑。
    *   将 `run_analysis.py` 和 `generate_final_analysis.py` 中的核心逻辑提取到 `glufire/analysis/` 模块中，或者作为 `glufire.cli` 的子命令。

**3. 实现命令行接口 (CLI)**

*   **目的:** 允许用户通过命令行运行项目的不同功能模块。
*   **实现方式:**
    *   使用 `argparse` (Python 标准库) 或 `click` (第三方库，更强大) 来构建 CLI。
    *   在 `glufire/cli.py` 中，定义一个主函数 `main()`，它将解析命令行参数并调用相应的模块函数。
    *   为 `and-gate`, `glu-metabolism`, `diffusion` 等模块创建子命令。
    *   每个子命令将负责调用其对应模块中的功能，并处理输入参数和输出结果。

**4. 模块化和重构**

*   **目的:** 确保每个模块都是独立的、可测试的，并且易于集成。
*   **修改内容:**
    *   检查 `and_gate.py`, `glu_metabolism_en.py`, `diffusion_pk_neurotoxicity.py` 中的代码，确保它们可以作为独立的函数或类被调用，而不是直接执行脚本。
    *   将 `data/` 和 `params/` 中的数据和参数加载逻辑封装到各自模块中，或者通过 `config_manager` 统一管理。

**5. 文档更新**

*   **目的:** 提供清晰的使用说明和 API 参考。
*   **修改内容:**
    *   更新 `README.md`，包含安装说明、CLI 使用示例和项目概览。
    *   为每个模块和 CLI 命令添加 docstrings。

**6. 测试**

*   **目的:** 确保所有功能正常工作。
*   **修改内容:**
    *   为每个模块和 CLI 命令编写单元测试和集成测试。

## CLI 使用说明

为了能够直接在命令行中运行 `glufire` 命令（例如 `glufire --help`），您需要确保 Python 解释器能够找到 `glufire` 包。这通常通过设置 `PYTHONPATH` 环境变量来实现。

**设置 PYTHONPATH 环境变量 (Windows PowerShell):**

在 PowerShell 中，您可以使用以下命令将项目根目录添加到 `PYTHONPATH`：

```powershell
$env:PYTHONPATH = "D:\pythonProject\glufire"
```

**重要提示:**

*   上述命令仅对当前的 PowerShell 会话有效。如果您关闭终端或打开新的终端，需要重新执行此命令。
*   如果您希望 `PYTHONPATH` 永久生效，您需要将其添加到系统环境变量中。具体操作请参考 Windows 操作系统文档。
*   在某些情况下，即使设置了 `PYTHONPATH`，系统可能仍然无法直接识别 `glufire` 命令。此时，您可以使用 `python -m glufire.cli <command> [OPTIONS]` 的形式来执行 CLI 命令，例如：
    ```bash
    python -m glufire.cli --help
    python -m glufire.cli and-gate --help
    ```
    这种方式更加健壮，因为它直接告诉 Python 解释器去哪里找并执行模块，不依赖于系统 `PATH` 环境变量中 `glufire` 可执行文件的位置。

---


`glufire` 包安装完成后，您可以通过命令行访问其核心功能。主命令是 `glufire`，它提供了三个子命令：`and-gate`、`glu-metabolism` 和 `diffusion`。

### 1. `glufire and-gate`

用于模拟 AND 逻辑门的行为。

**用法:**
```bash
glufire and-gate [OPTIONS]
```

**选项:**
*   `--oxygen <OXYGEN>`: 氧气浓度 (例如: `0.01` 代表 1%)。
*   `--temperature <TEMPERATURE>`: 温度 (例如: `42` 摄氏度)。
*   `--plot`: 显示结果图表。
*   `--save-path <SAVE_PATH>`: 保存结果图表的路径。
*   `--help`: 显示帮助信息。

**示例:**
模拟氧气浓度为 0.01，温度为 42 时的 AND 逻辑门，并显示图表：
```bash
glufire and-gate --oxygen 0.01 --temperature 42 --plot
```

### 2. `glufire glu-metabolism`

用于模拟谷氨酸代谢过程。

**用法:**
```bash
glufire glu-metabolism [OPTIONS]
```

**选项:**
*   `--t7-activity <T7_ACTIVITY>`: T7 聚合酶活性。
*   `--duration <DURATION>`: 模拟时长 (例如: `24` 小时)。
*   `--plot`: 显示结果图表。
*   `--save-path <SAVE_PATH>`: 保存结果图表的路径。
*   `--help`: 显示帮助信息。

**示例:**
模拟 T7 活性为 0.8，时长为 12 小时的谷氨酸代谢，并保存图表：
```bash
glufire glu-metabolism --t7-activity 0.8 --duration 12 --save-path "results/glu_metabolism_simulation.png"
```

### 3. `glufire diffusion`

用于模拟谷氨酸在体内的扩散和评估神经毒性。

**用法:**
```bash
glufire diffusion [OPTIONS]
```

**选项:**
*   `--glu-concentration <GLU_CONCENTRATION>`: 初始谷氨酸浓度。
*   `--duration <DURATION>`: 模拟时长 (例如: `48` 小时)。
*   `--plot`: 显示结果图表。
*   `--save-path <SAVE_PATH>`: 保存结果图表的路径。
*   `--help`: 显示帮助信息。

**示例:**
模拟初始谷氨酸浓度为 100，时长为 24 小时的扩散过程，并显示图表：
```bash
glufire diffusion --glu-concentration 100 --duration 24 --plot
```