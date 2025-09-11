import json
import os

class ConfigManager:
    def __init__(self):
        self.config = {}
        self._load_default_params()

    def _load_default_params(self):
        # 谷氨酸代谢模型 (GluMetabolismModel) 参数
        self.config['glu_metabolism'] = {
            'k_prod_max': 10.0,
            'K_t7': 100.0,
            'k_export_max': 5.0,
            'K_export': 1.0,
            'k_dilution': 0.1,
            'V_intra_over_V_extra': 0.01,
            'k_syn_icd': 0.5,
            'k_deg_icd': 0.1,
            'Vmax_icd': 50.0,
            'K_icd': 0.1,
            'k_syn_gdhA': 0.5,
            'k_deg_gdhA': 0.1,
            'Vmax_gdhA': 50.0,
            'K_gdhA': 0.1,
            'n_hill': 2.0
        }
        self.config['glu_metabolism_en'] = {
            'V_max_glc': 10.0,
            'K_m_glc': 1.0,
            'f_TCA': 0.6,
            'V_max_base_ICD': 25.0,
            'K_m_ICD': 0.029,
            'V_max_base_GDH': 30.0,
            'K_m_AKG': 0.64,
            'K_m_NH4': 1.1,
            'K_m_NADPH': 0.04,
            'k_PPP': 0.4,
            'y_ICD_NADPH': 1.0,
            'lambda_NADPH': 2.0,
            'NADPH_set': 0.15,
            'k_sec_base': 0.8,
            'k_maintenance': 0.08,
            'mu_max': 0.5,
            'K_T7': 800.0,
            'n_hill': 3.0,
            'tau_enzyme': 0.05,
            'fold_ICD_max': 1000.0,
            'fold_GDH_max': 1500.0,
            'homeostasis_strength': 2.0,
            'accum_threshold': 55.0,
            'export_accum_suppression': 0.05,
            'postshock_export_boost': 10.0,
            'extracellular_clearance_rate': 0.5,
            'export_decay_rate': 0.8,
            'Glu_target': 20.0
        }

    def load_config_from_json(self, json_path):
        if os.path.exists(json_path):
            with open(json_path, 'r', encoding='utf-8') as f:
                user_config = json.load(f)
            self._deep_update(self.config, user_config)
        else:
            print(f"Warning: Config file not found at {json_path}. Using default parameters.")

    def _deep_update(self, target, source):
        for k, v in source.items():
            if isinstance(v, dict) and k in target and isinstance(target[k], dict):
                target[k] = self._deep_update(target[k], v)
            else:
                target[k] = v
        return target

    def get_params(self, module_name):
        return self.config.get(module_name, {})

# 示例用法
if __name__ == '__main__':
    config_manager = ConfigManager()
    # 可以选择从JSON文件加载配置
    # config_manager.load_config_from_json('path/to/your/config.json')
    
    glu_params = config_manager.get_params('glu_metabolism')
    print("Glu metabolism parameters:", glu_params)
