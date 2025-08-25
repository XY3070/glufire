
# Skeleton for COBRApy dFBA linking expression to flux
# Requires: pip install cobra
import json, numpy as np, pandas as pd
from cobra.io import load_json_model

# --- User: set model path ---
MODEL_JSON = "models/e_coli_core.json"  # put your model here

# Load model
model = load_json_model(MODEL_JSON)

# Map: icd/gdhA expression level -> flux bounds (toy mapping; replace with your calibration)
def apply_expression_to_bounds(model, icd_scale=1.0, gdh_scale=1.0):
    # Example reaction IDs in E. coli core: ICDHyr (isocitrate dehydrogenase), GLUDy (glutamate dehydrogenase NADP-dependent)
    rxn_map = {"ICDHyr": icd_scale, "GLUDy": gdh_scale}
    for rid, scale in rxn_map.items():
        if rid in model.reactions:
            rxn = model.reactions.get_by_id(rid)
            # widen upper bound proportional to expression (placeholder logic)
            rxn.upper_bound = 1000.0 * float(scale)
    return model

# Simple dFBA loop with fixed time-steps and substrate drains
def run_dfba(T=24.0, dt=0.5, icd_scale=2.0, gdh_scale=3.0):
    t = 0.0; records = []
    while t < T:
        apply_expression_to_bounds(model, icd_scale, gdh_scale)
        sol = model.optimize()
        mu = sol.objective_value
        q_glu = 0.0
        # if you add an exchange reaction EX_glu(e), read its flux:
        if "EX_glu__L_e" in model.reactions:
            q_glu = model.reactions.get_by_id("EX_glu__L_e").flux
        records.append({"time_h": t, "mu": mu, "q_glu": q_glu})
        t += dt
    df = pd.DataFrame(records)
    df.to_csv("dfba_timeseries.csv", index=False)
    print("Saved dfba_timeseries.csv")

if __name__ == "__main__":
    run_dfba()
