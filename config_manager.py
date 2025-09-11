import json
import os

class ConfigManager:
    def __init__(self):
        self.config = {}
        self._load_default_params()

    def _load_default_params(self):
        # 谷氨酸代谢模型 (GluMetabolismModel) 参数
        self.config['glu_metabolism'] = {
            'k_prod_max': 50.0,  # 最大谷氨酸生产速率 (mM/hr)
            'K_t7': 500.0,       # T7活性达到半最大生产速率时的值 (AU)
            'k_export_max': 100.0, # 最大谷氨酸分泌速率 (mM/hr)
            'K_export': 10.0,     # 细胞内谷氨酸浓度达到半最大分泌速率时的值 (mM)
            'k_dilution': 0.1,   # 细胞生长或降解导致的稀释/降解速率 (1/hr)
            'V_intra_over_V_extra': 0.01, # 细胞内总体积与细胞外总体积的比率
            'k_syn_icd': 2.0,    # Icd合成速率 (1/hr)
            'k_syn_gdhA': 2.0,   # gdhA合成速率 (1/hr)
            'k_deg_icd': 0.2,    # Icd降解速率 (1/hr)
            'k_deg_gdhA': 0.2,   # gdhA降解速率 (1/hr)
            'Vmax_icd': 100.0,    # Icd最大催化速率 (mM/hr)
            'K_icd': 5.0,        # Icd底物常数 (mM)
            'Vmax_gdhA': 100.0,   # gdhA最大催化速率 (mM/hr)
            'K_gdhA': 5.0,       # gdhA底物常数 (mM)
            'n_hill': 4.0        # Hill系数，用于增强开关效应
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
