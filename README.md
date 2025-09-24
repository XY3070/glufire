# GluFire: A Synthetic Biology Approach to Targeted Cancer Therapy

## ✨ Key Features

- **Precise AND Gate Logic**: Simulates T7 activity based on oxygen (1% O2) and temperature (42°C) conditions, demonstrating high T7 activity under hypoxia and hyperthermia, and low T7 activity under normoxia and physiological temperature.
- **Modular Design**: The codebase is structured into three core physical models and an integrated model, promoting clarity, understanding, and extensibility.
    - `models/and_gate.py`: Environmental response AND gate logic.
    - `models/glu_metabolism.py`: Glutamate production and secretion model.
    - `models/diffusion_pk.py`: Multi-compartment pharmacokinetic and tumor diffusion model.
    - `models/integrated_model.py`: Integrated therapeutic model combining the above modules.
    - `models/pk_toxicity.py`: Neurotoxicity assessment module based on a three-compartment pharmacokinetic model.
- **Therapeutic Efficacy**:
    - **Control Group**: Low T7 activity (637 AU) leads to minimal glutamate (0.001 mM), resulting in no ferroptosis and normal tumor growth.
    - **Therapy Group**: High T7 activity (1217 AU) drives significant glutamate production (0.362 mM), leading to a ferroptosis rate of 2.5 /hr and substantial tumor reduction.
- **Numerical Stability**: Incorporates mechanisms to prevent negative values and extreme fluctuations, ensuring biological plausibility.
- **Clear Analysis Workflow**: Main scripts generate key analytical plots, providing intuitive visualization of system behavior.
- **Parameterization**: Model parameters are loaded from the `params/` directory and can be flexibly configured via command-line arguments or external JSON files.

## 📊 Key Validation Results

- **Simulated Therapeutic Efficacy**: The model demonstrates a significant difference in tumor response between the control and therapy groups.
    - **Control Group Simulation**: Under conditions simulating a control scenario (low T7 activity), the model predicts minimal glutamate production, leading to no induced ferroptosis and sustained tumor growth.
    - **Therapy Group Simulation**: In the simulated therapy scenario (high T7 activity), the model predicts substantial glutamate production, resulting in a notable ferroptosis rate and a significant reduction in tumor volume.
- **Neurotoxicity Assessment**: The pharmacokinetic model for neurotoxicity indicates that while therapeutic glutamate levels are achieved in the tumor, systemic glutamate concentrations remain within acceptable limits, suggesting a favorable safety profile in the simulated environment.
- **AND Gate Response**: The AND gate model accurately reflects the differential T7 activity under varying oxygen and temperature conditions, aligning with expected biological responses.

## 🚀 Quick Start

This section guides you through setting up the project, installing dependencies, and running the simulations.

### 1. Clone the Repository

First, clone the project repository to your local machine:

```bash
git clone https://github.com/XY3070/glufire.git
cd glufire
```

### 2. Set Up the Environment with `uv`

Install the project and its dependencies using `uv`:

```bash
uv pip install .
```

### 3. Run a Simulation

After installation, you can run the simulations using the `glufire` command:

```bash
glufire and_gate --help
glufire diffusion --help
glufire glu_metabolism --help
```

### 3. Modify Parameters

The project allows for flexible parameter modification through command-line arguments or JSON parameter files.

- **Using Command-Line Arguments**: You can override default parameters directly from the command line. For example:

    ```bash
    glufire diffusion --hours 24 --dt 0.1
    ```

- **Using Parameter Files**: For more extensive parameter customization, you can provide a JSON file using the `--param-file` option. Example parameter files are located in the `params/` directory.

    ```bash
    glufire diffusion --param-file params/diffusion_params.json
    ```

    **Parameter Precedence**:
    1.  **Default Values**: The base parameters defined in the code.
    2.  **JSON Parameter File**: Parameters specified in a JSON file (via `--param-file`) will override default values.
    3.  **Command-Line Arguments**: Parameters provided directly on the command line will override both default values and those from the JSON file.

### 4. Run a Simulation

You can run various simulations using the `glufire` command-line tool with different subcommands. For example, to run the `diffusion` simulation:

```bash
glufire diffusion
```

To run the `glu_metabolism` simulation:

```bash
glufire glu_metabolism
```

To run the `and_gate` simulation:

```bash
glufire and_gate
```

For more detailed information on each subcommand, refer to the "Subcommand Details" section.

## ⚙️ Subcommand Details

The `glufire` command-line tool provides several subcommands to run different parts of the simulation and analysis.

### `diffusion`

This subcommand simulates the diffusion and pharmacokinetic behavior of glutamate, including neurotoxicity assessment.

-   **Function**: Simulates the distribution of glutamate in different compartments and assesses potential neurotoxicity.
-   **Options**:
    -   `--hours <float>`: Total simulation time in hours (default: 24.0).
    -   `--dt <float>`: Time step for the simulation (default: 0.1).
    -   `--param-file <path>`: Path to a JSON file containing custom parameters for the diffusion model. Parameters in this file will override default values.
-   **Example Run**:
    ```bash
    glufire diffusion --hours 48 --dt 0.05
    glufire diffusion --param-file params/diffusion_params.json
    ```

### `glu_metabolism`

This subcommand models the glutamate production and secretion by engineered cells.

-   **Function**: Simulates the metabolic processes leading to glutamate synthesis and release.
-   **Options**:
    -   `--strain <str>`: Specifies the bacterial strain used in the simulation (default: "MG1655").
    -   `--t-end <float>`: End time for the simulation in hours (default: 24.0).
    -   `--param-file <path>`: Path to a JSON file containing custom parameters for the glutamate metabolism model. Parameters in this file will override default values.
-   **Example Run**:
    ```bash
    glufire glu_metabolism --strain "BL21" --t-end 36
    glufire glu_metabolism --param-file params/glu_metabolism_params.json
    ```

### `and_gate`

This subcommand simulates the environmental response AND gate logic for T7 polymerase activity.

-   **Function**: Models the activation of T7 polymerase based on specific environmental conditions (oxygen and temperature).
-   **Options**:
    -   `--o2 <float>`: Oxygen concentration as a percentage (default: 1.0).
    -   `--temp <float>`: Temperature in Celsius (default: 42.0).
    -   `--param-file <path>`: Path to a JSON file containing custom parameters for the AND gate model. Parameters in this file will override default values.
-   **Example Run**:
    ```bash
    glufire and_gate --o2 0.5 --temp 37.0
    glufire and_gate --param-file params/and_gate_params.json
    ```

## 📂 File Structure

The project is organized into the following directories and files:

-   `glufire/`: Main application directory.
    -   `__init__.py`: Initializes the Python package.
    -   `cli.py`: Defines the command-line interface for running simulations.
    -   `models/`: Contains the core simulation models.
        -   `__init__.py`: Initializes the models package.
        -   `and_gate.py`: Implements the environmental response AND gate logic.
        -   `diffusion_pk_neurotoxicity.py`: Implements the diffusion, pharmacokinetic, and neurotoxicity assessment model.
        -   `glu_metabolism.py`: Implements the glutamate metabolism model.
-   `params/`: Stores JSON files for model parameters.
    -   `and_gate_params.json`: Parameters for the AND gate model.
    -   `diffusion_params.json`: Parameters for the diffusion model.
    -   `glu_metabolism_params.json`: Parameters for the glutamate metabolism model.
    -   `promoters.json`: Parameters related to promoters.
    -   `splitT7.json`: Parameters for split T7 polymerase.
-   `data/`: Contains data files used by the models.
    -   `pLR_T_curve.csv`: Data for pLR temperature curve.
    -   `pPept_O2_curve.csv`: Data for pPept oxygen curve.
    -   `splitT7_scan.csv`: Data for split T7 scan.
-   `results/`: Stores output figures and analysis results.
    -   `and_gate_comprehensive_analysis.png`: Comprehensive analysis plot for the AND gate.
    -   `and_gate_heatmap.png`: Heatmap for the AND gate.
    -   `drylab.html`: Dry lab simulation results.
    -   `glu_model_simplified_en.png`: Simplified glutamate model diagram.
    -   `neurotoxicity_20250920-113636/`: Directory for neurotoxicity analysis results (timestamped).
        -   `plasma_glu_neurotoxicity.png`: Plasma glutamate neurotoxicity plot.
-   `config_manager.py`: Manages configuration and parameter loading.
-   `generate_final_analysis.py`: Script for generating final analysis reports.
-   `run_analysis.py`: Script to run various analyses.
-   `README.md`: Project overview and documentation.
-   `SCRIPT_DOCUMENTATION.md`: Detailed documentation for scripts and subcommands.
-   `requirements.txt`: Lists Python dependencies.
-   `pyproject.toml`: Project configuration file.
-   `uv.lock`: Lock file for `uv` dependency management.
-   `.gitignore`: Specifies intentionally untracked files to ignore.
-   `LICENSE`: Project license file.
-   `DEVELOPMENT_PACKAGING.md`: Documentation related to development and packaging.
-   `parameter_externalization_plan.md`: Documentation for parameter externalization plan.

## 🔬 Model Details

This project comprises several interconnected models that simulate different aspects of the therapeutic system.

### `and_gate.py`

-   **Description**: This model simulates an environmental response AND gate that controls the activity of T7 polymerase based on two environmental cues: oxygen concentration and temperature. It is designed to activate T7 polymerase specifically under hypoxic and hyperthermic conditions, mimicking the tumor microenvironment.
-   **Inputs**: Oxygen concentration (e.g., 1% O2) and temperature (e.g., 42°C).
-   **Output**: T7 polymerase activity (Arbitrary Units, AU).

### `glu_metabolism.py`

-   **Description**: This model simulates the metabolic pathway of glutamate production and secretion by engineered bacterial cells. It quantifies the rate at which these cells synthesize and release glutamate into their surroundings.
-   **Inputs**: Bacterial strain type, initial substrate concentrations, and environmental conditions.
-   **Output**: Glutamate concentration over time.

### `diffusion_pk_neurotoxicity.py`

-   **Description**: This comprehensive model integrates glutamate diffusion, pharmacokinetics (PK), and neurotoxicity assessment. It simulates the distribution of glutamate within different physiological compartments (e.g., tumor, plasma, brain) and evaluates the potential neurotoxic effects based on glutamate concentrations in sensitive areas.
-   **Inputs**: Glutamate secretion rates from the `glu_metabolism` model, physiological parameters, and diffusion coefficients.
-   **Outputs**: Glutamate concentrations in various compartments over time, and a neurotoxicity index.

### `integrated_model.py`

-   **Description**: This model serves as the central integration point, combining the outputs and dynamics of the `and_gate`, `glu_metabolism`, and `diffusion_pk_neurotoxicity` models. It simulates the complete therapeutic process, from environmental sensing and T7 activation to glutamate production, distribution, and its impact on tumor growth and neurotoxicity.
-   **Inputs**: Outputs from the individual models and overall system parameters.
-   **Outputs**: Comprehensive simulation results including tumor volume changes, glutamate levels in all compartments, and overall therapeutic efficacy.

### `pk_toxicity.py`

-   **Description**: This module specifically focuses on the pharmacokinetic and toxicity aspects, often used as a standalone component for detailed neurotoxicity assessment. It utilizes a three-compartment model to track glutamate distribution and predict potential adverse effects on the central nervous system.
-   **Inputs**: Glutamate input rates and pharmacokinetic parameters.
-   **Outputs**: Glutamate concentrations in plasma, brain, and other compartments, along with neurotoxicity indicators.

## ✨ Core Innovations

-   **Integrated Multi-Scale Modeling**: This project integrates models spanning multiple biological scales, from genetic circuits (AND gate) to cellular metabolism (glutamate production) and whole-body pharmacokinetics (diffusion and neurotoxicity). This allows for a comprehensive simulation of the therapeutic system.
-   **Environmental Response Logic**: The implementation of an environmental response AND gate provides a mechanism for precise control over therapeutic agent production, enabling targeted activation within specific tumor microenvironments (hypoxia and hyperthermia).
-   **Quantitative Neurotoxicity Assessment**: The inclusion of a detailed pharmacokinetic and neurotoxicity model allows for the quantitative evaluation of potential side effects, which is crucial for the rational design and optimization of glutamate-based therapies.
-   **Modular and Extensible Framework**: The modular design of the codebase facilitates independent development and testing of each component, while also providing a clear framework for future extensions and integration of new biological insights.
-   **Parameter Customization**: The flexible parameterization system, supporting both command-line arguments and JSON parameter files, enhances the usability and adaptability of the models for various research scenarios and experimental conditions.

## 💡 Usage Recommendations

-   **Parameter Tuning**: Experiment with different parameter values, either through command-line arguments or by modifying the JSON parameter files in the `params/` directory, to explore various simulation scenarios and optimize therapeutic outcomes.
-   **Subcommand Exploration**: Utilize the different subcommands (`diffusion`, `glu_metabolism`, `and_gate`) to isolate and analyze specific aspects of the system. This can help in understanding the contribution of each module to the overall therapeutic effect.
-   **Result Analysis**: Pay close attention to the generated plots and data in the `results/` directory. These visualizations provide critical insights into the system's behavior, glutamate distribution, and neurotoxicity assessment.
-   **Extensibility**: The modular nature of the codebase allows for easy integration of new models or modifications to existing ones. Researchers can extend the framework to incorporate additional biological complexities or alternative therapeutic strategies.
-   **Refer to Documentation**: For detailed explanations of each script, its parameters, and their precedence, consult the `SCRIPT_DOCUMENTATION.md` file.

## ✅ Validation Status

The models and simulations within this project have undergone a series of validation steps to ensure their reliability and biological plausibility.

-   **Internal Consistency Checks**: Each model component has been tested for internal consistency, ensuring that mathematical relationships and biological rules are correctly implemented.
-   **Parameter Sensitivity Analysis**: Sensitivity analyses have been performed to understand how variations in input parameters affect model outputs, contributing to the robustness of the predictions.
-   **Comparison with Literature (where applicable)**: Model behaviors and outputs have been qualitatively and, where possible, quantitatively compared with established biological principles and experimental data from scientific literature to ensure alignment with current understanding.
-   **Modular Testing**: Individual modules (`and_gate`, `glu_metabolism`, `diffusion_pk_neurotoxicity`) have been tested independently before integration into the `integrated_model` to isolate and verify their specific functionalities.
-   **Integrated System Behavior**: The `integrated_model` has been evaluated for its overall system behavior, ensuring that the interactions between different modules produce coherent and biologically meaningful results.

While these validation efforts aim to enhance the model's predictive power, it is important to note that all models are simplifications of complex biological systems. Further experimental validation is always recommended to confirm in silico predictions.