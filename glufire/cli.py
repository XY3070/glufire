import click
import numpy as np
import os
from pathlib import Path
import json # Import the json module

from glufire.models.and_gate import SimpleANDGate
from glufire.models.diffusion_pk_neurotoxicity import simulate_three_comp_pk, assess_neurotoxicity, PKParams, ToxicityThresholds, trapezoid_flux, _plot_and_save
from glufire.models.glu_metabolism import GluModel, create_heat_shock_protocol


@click.group()
def main():
    """
    glufire CLI for simulating biological models.
    """
    pass


@main.command()
@click.option('--o2', default=1.0, type=float, help='Oxygen level in percent (e.g., 1.0 for 1%).')
@click.option('--temp', default=42.0, type=float, help='Temperature in Celsius (e.g., 42.0).')
@click.option('--diagnose', is_flag=True, help='Run quick diagnosis and print T7 activity for predefined conditions.')
@click.option('--plot', is_flag=True, help='Generate and show comprehensive analysis plot.')
@click.option('--save-path', type=click.Path(), help='Path to save the plot (e.g., results/and_gate_analysis.png).')
@click.option('--param-file', type=click.Path(exists=True), help='Path to a JSON file containing AND-gate parameters.')
def and_gate(o2, temp, diagnose, plot, save_path, param_file):
    """
    Simulate the AND-gate logic model.
    """
    click.echo(f"Running AND-gate simulation with O2={o2}%, Temp={temp}°C")

    params_from_file = {}
    if param_file:
        with open(param_file, 'r') as f:
            params_from_file = json.load(f)

    # Default parameters
    default_params = {
        'o2': 1.0,
        'temp': 42.0,
    }

    # Merge parameters: default_params < params_from_file < command_line_args
    merged_params = default_params.copy()
    merged_params.update(params_from_file)

    cmd_line_args = {
        'o2': o2, 'temp': temp
    }
    for key, value in cmd_line_args.items():
        if value is not None:
            merged_params[key] = value

    # Assign merged parameters to local variables
    o2 = merged_params['o2']
    temp = merged_params['temp']

    gate = SimpleANDGate()

    if diagnose:
        gate.quick_diagnose()

    if plot or save_path:
        if save_path:
            os.makedirs(Path(save_path).parent, exist_ok=True)
        click.echo("Generating comprehensive analysis plot...")
        gate.create_comprehensive_analysis(save_path=save_path, show_plot=plot)
        click.echo("AND-gate simulation completed.")


@main.command()
@click.option('--strain', default='engineered', type=click.Choice(['engineered', 'wildtype']), help='Strain type for glutamate metabolism model.')
@click.option('--t-end', default=48.0, type=float, help='Total simulation time in hours.')
@click.option('--dt', default=0.1, type=float, help='Time step for simulation.')
@click.option('--shock-start', default=8.0, type=float, help='Heat shock start time in hours.')
@click.option('--shock-duration', default=4.0, type=float, help='Heat shock duration in hours.')
@click.option('--t7-low', default=50.0, type=float, help='Low T7 activity during non-shock periods.')
@click.option('--t7-high', default=3000.0, type=float, help='High T7 activity during heat shock periods.')
@click.option('--plot', is_flag=True, help='Generate and show plots for glutamate metabolism.')
@click.option('--save-path', type=click.Path(), help='Path to save the plot (e.g., results/glu_metabolism_analysis.png).')
@click.option('--param-file', type=click.Path(exists=True), help='Path to a JSON file containing glutamate metabolism parameters.')
def glu_metabolism(strain, t_end, dt, shock_start, shock_duration, t7_low, t7_high, plot, save_path, param_file):
    """
    Simulate the glutamate metabolism model.
    """
    click.echo(f"Running glutamate metabolism simulation for {strain} strain...")

    params_from_file = {}
    if param_file:
        with open(param_file, 'r') as f:
            params_from_file = json.load(f)

    # Default parameters
    default_params = {
        'strain': 'engineered',
        't_end': 48.0,
        'dt': 0.1,
        'shock_start': 8.0,
        'shock_duration': 4.0,
        't7_low': 50.0,
        't7_high': 3000.0,
    }

    # Merge parameters: default_params < params_from_file < command_line_args
    merged_params = default_params.copy()
    merged_params.update(params_from_file)

    cmd_line_args = {
        'strain': strain, 't_end': t_end, 'dt': dt, 'shock_start': shock_start,
        'shock_duration': shock_duration, 't7_low': t7_low, 't7_high': t7_high
    }
    for key, value in cmd_line_args.items():
        if value is not None:
            merged_params[key] = value

    # Assign merged parameters to local variables
    strain = merged_params['strain']
    t_end = merged_params['t_end']
    dt = merged_params['dt']
    shock_start = merged_params['shock_start']
    shock_duration = merged_params['shock_duration']
    t7_low = merged_params['t7_low']
    t7_high = merged_params['t7_high']

    model = GluModel(strain_type=strain)

    heat_shock_protocol = create_heat_shock_protocol(
        shock_start=shock_start,
        shock_duration=shock_duration,
        t7_low=t7_low,
        t7_high=t7_high
    )

    t, solution = model.simulate(t7_activity_func=heat_shock_protocol, t_end=t_end, dt=dt)
    metrics = model.analyze_performance(t, solution)

    click.echo("\n=== Performance Analysis Results ===")
    click.echo(f"  Intracellular Glu Peak: {metrics['max_intracellular_glu']:.2f} mM")
    click.echo(f"  Extracellular Glu Peak: {metrics['max_extracellular_glu']:.2f} mM")
    click.echo(f"  Final Intracellular Glu: {metrics['final_intracellular_glu']:.2f} mM")
    click.echo(f"  Intracellular Glu>=45mM: {'Success' if metrics['targets_met']['intracellular_peak_50mM'] else 'Failed'}")
    click.echo(f"  Extracellular Glu>=30mM: {'Success' if metrics['targets_met']['extracellular_peak_30mM'] else 'Failed'}")
    click.echo(f"  Final Recovery~20mM: {'Success' if metrics['targets_met']['final_recovery_20mM'] else 'Failed'}")

    if plot or save_path:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle(f'{strain.capitalize()} Strain Glutamate Metabolism Analysis', fontsize=16, fontweight='bold')

        # Intracellular glutamate
        axes[0].plot(t, solution[:, 4], color='#FF8200', linewidth=2, label='Intracellular Glu')
        axes[0].axhline(y=45, color='#A0522D', linestyle=':', alpha=0.7, label='Target>=45mM')
        axes[0].axhline(y=20, color='green', linestyle=':', alpha=0.7, label='Recovery~20mM')
        axes[0].axvline(x=shock_start, color='#A0522D', linestyle='--', alpha=0.8, linewidth=1.5, label='Heat Shock Start')
        axes[0].set_xlabel('Time (h)')
        axes[0].set_ylabel('Intracellular Glu (mM)')
        axes[0].set_title('Intracellular Glutamate Dynamics')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        # Extracellular glutamate
        axes[1].plot(t, solution[:, 7], color='#515B87', linewidth=2, label='Extracellular Glu')
        axes[1].axhline(y=30, color='#A0522D', linestyle=':', alpha=0.7, label='Target>=30mM')
        axes[1].axvline(x=shock_start, color='#A0522D', linestyle='--', alpha=0.8, linewidth=1.5, label='Heat Shock Start')
        axes[1].set_xlabel('Time (h)')
        axes[1].set_ylabel('Extracellular Glu (mM)')
        axes[1].set_title('Extracellular Glutamate Dynamics')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        # Key metabolites
        axes[2].plot(t, solution[:, 3], color='#FF8200', linewidth=2, label='AKG')
        axes[2].plot(t, solution[:, 5], color='#515B87', linewidth=2, label='NADPH')
        axes[2].set_xlabel('Time (h)')
        axes[2].set_ylabel('Concentration (mM)')
        axes[2].set_title('Key Metabolites')
        axes[2].legend()
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            os.makedirs(Path(save_path).parent, exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            click.echo(f"Plot saved to: {save_path}")
        if plot:
            plt.show()
        else:
            plt.close(fig)
    click.echo("Glutamate metabolism simulation completed.")


@main.command()
@click.option('--hours', default=None, type=float, help='Total simulation horizon (h).')
@click.option('--dt', default=None, type=float, help='Time step for output grid (h).')
@click.option('--baseline', default=None, type=float, help='Baseline concentration for all compartments (µM).')
@click.option('--tumor-init-mm', default=None, type=float, help='Initial tumor glutamate (mM), worst-case.')
@click.option('--k-bt', default=None, type=float, help='Blood<->tumor exchange (h^-1).')
@click.option('--k-b-clr', default=None, type=float, help='Plasma clearance (h^-1).')
@click.option('--k-t-clr', default=None, type=float, help='Tumor clearance (h^-1).')
@click.option('--secretion-peak', default=None, type=float, help='Tumor secretion peak (µmol/h). 0 disables.')
@click.option('--t-on', default=None, type=float, help='Secretion start (h).')
@click.option('--t-off', default=None, type=float, help='Secretion end (h).')
@click.option('--ramp', default=None, type=float, help='Ramp duration for on/off (h).')
@click.option('--plot', is_flag=True, help='Generate and show neurotoxicity plot.')
@click.option('--save-path', type=click.Path(), help='Path to save the plot (e.g., results/neurotoxicity_analysis.png).')
@click.option('--title', default=None, help='Figure title.')
@click.option('--param-file', type=click.Path(exists=True), help='Path to a JSON file containing diffusion parameters.')
def diffusion(hours, dt, baseline, tumor_init_mm, k_bt, k_b_clr, k_t_clr, secretion_peak, t_on, t_off, ramp, plot, save_path, title, param_file):
    """
    Simulate the diffusion and neurotoxicity model.
    """
    click.echo("Running diffusion and neurotoxicity simulation...")

    params_from_file = {}
    if param_file:
        with open(param_file, 'r') as f:
            params_from_file = json.load(f)

    # Default parameters
    default_params = {
        'hours': 48.0,
        'dt': 0.1,
        'baseline': 50.0,
        'tumor_init_mm': 30.0,
        'k_bt': 1e-3,
        'k_b_clr': 0.5,
        'k_t_clr': 0.05,
        'secretion_peak': 0.0,
        't_on': 8.0,
        't_off': 34.0,
        'ramp': 2.0,
        'title': 'Plasma Glutamate vs Neurotoxicity Thresholds'
    }

    # Merge parameters: default_params < params_from_file < command_line_args
    # Command-line arguments that are not None override file parameters, which override defaults
    merged_params = default_params.copy()
    merged_params.update(params_from_file)

    cmd_line_args = {
        'hours': hours, 'dt': dt, 'baseline': baseline, 'tumor_init_mm': tumor_init_mm,
        'k_bt': k_bt, 'k_b_clr': k_b_clr, 'k_t_clr': k_t_clr, 'secretion_peak': secretion_peak,
        't_on': t_on, 't_off': t_off, 'ramp': ramp, 'title': title
    }
    for key, value in cmd_line_args.items():
        if value is not None:
            merged_params[key] = value

    # Assign merged parameters to local variables
    hours = merged_params['hours']
    dt = merged_params['dt']
    baseline = merged_params['baseline']
    tumor_init_mm = merged_params['tumor_init_mm']
    k_bt = merged_params['k_bt']
    k_b_clr = merged_params['k_b_clr']
    k_t_clr = merged_params['k_t_clr']
    secretion_peak = merged_params['secretion_peak']
    t_on = merged_params['t_on']
    t_off = merged_params['t_off']
    ramp = merged_params['ramp']
    title = merged_params['title']

    sim_hours = hours
    baseline_uM = baseline
    Ct0_uM = tumor_init_mm * 1000.0  # convert mM -> µM

    t_h = np.linspace(0.0, sim_hours, int(sim_hours / dt) + 1)

    if secretion_peak > 0.0:
        S_t = np.array([
            trapezoid_flux(t,
                           t_on=t_on,
                           t_plateau=t_on + ramp,
                           t_off=t_off,
                           ramp_hours=ramp,
                           peak_umol_per_h=secretion_peak)
            for t in t_h
        ], dtype=float)
    else:
        S_t = np.zeros_like(t_h)

    params = PKParams(
        k_bt=k_bt,
        k_t_clr=k_t_clr,
        k_b_clr=k_b_clr,
    )

    Cb, Ct, Cn = simulate_three_comp_pk(
        t_grid_h=t_h,
        S_t_umol_per_h=S_t,
        params=params,
        Cb0_uM=baseline_uM,
        Ct0_uM=Ct0_uM,
        Cn0_uM=baseline_uM,
        baseline_uM=baseline_uM,
    )

    Cb_ctrl, Ct_ctrl, Cn_ctrl = simulate_three_comp_pk(
        t_grid_h=t_h,
        S_t_umol_per_h=np.zeros_like(S_t),
        params=params,
        Cb0_uM=baseline_uM,
        Ct0_uM=baseline_uM,
        Cn0_uM=baseline_uM,
        baseline_uM=baseline_uM,
    )

    thr = ToxicityThresholds()
    report_treat = assess_neurotoxicity(Cb, t_h, thr)

    click.echo("\n=== Neurotoxicity Risk Report (treatment) ===")
    for k, v in report_treat.items():
        click.echo(f"{k}: {v}")

    if plot or save_path:
        if save_path:
            os.makedirs(Path(save_path).parent, exist_ok=True)
        click.echo("Generating neurotoxicity plot...")
        _plot_and_save(
            t_h=t_h,
            Cb=Cb, Cb_ctrl=Cb_ctrl,
            baseline_uM=baseline_uM,
            thresholds=thr,
            t_on=t_on,
            t_off=t_off,
            out_png=Path(save_path) if save_path else Path("results/neurotoxicity_analysis.png"),
            title=title,
        )
        click.echo("Diffusion and neurotoxicity simulation completed.")


if __name__ == '__main__':
    main()