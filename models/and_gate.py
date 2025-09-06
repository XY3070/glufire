#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AND门逻辑模块 (AND Gate Logic Module)

低氧 + 高温 => 高 T7 活性
高氧 或 低温 任一不满足 => 低 T7 活性

修复说明 (2025-09-03):
- 之前 promoters.json 中使用字段 "_mode": "rep" / "_mode": "act"，但代码仅检查 'type'，导致 pPept 被当成激活型 → 氧气高时反而提升表达，逻辑颠倒。
- 本次更新在参数加载阶段归一化: 如果存在 _mode 或 mode 则写入标准键 'type'。
- Hill 函数增加冗余健壮性判断 (支持 rep/repressor/inh 关键词)。
- 新增 quick_diagnose() 便于快速打印在典型条件下的表达与T7输出，帮助定位温度 K 不匹配问题 (例如 pLR K=44.7 而实验设 42°C 会导致温度臂未充分激活 → 可考虑下调 K 或提高设定温度)。
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import odeint

# ------------------------------------------------------------------
# 参数加载
# ------------------------------------------------------------------
def load_model_parameters(promoter_params_path="params/promoters.json", splitT7_params_path="params/splitT7.json"):
    base_path = Path(__file__).resolve().parents[1]
    # 启动子参数
    try:
        with open(base_path / promoter_params_path, 'r', encoding='utf-8') as f:
            promoter_params = json.load(f)
    except Exception:
        # 默认: pPept 为低氧诱导 (repressor of O2), pLR 为高温激活
        promoter_params = {
            "pPept": {"type": "rep", "beta": 1200.0, "K": 5.0, "n": 2.0, "leaky": 50.0},  # K=5% O2
            "pLR":   {"type": "act", "beta": 1500.0, "K": 40.0, "n": 8.0,  "leaky": 80.0}  # 阈值 40°C
        }

    # 统一各启动子 mode/type 字段
    for name, p in promoter_params.items():
        # 优先已有规范 'type'
        mode = p.get('type') or p.get('_mode') or p.get('mode')
        if mode:
            mode_l = str(mode).lower()
            if mode_l.startswith('rep') or 'inh' in mode_l:
                p['type'] = 'rep'
            else:
                p['type'] = 'act'
        else:
            # 若未提供且名称含有常见抑制关键词，可在此扩展规则。目前默认激活。
            p.setdefault('type', 'act')

    # 分裂 T7 组装参数
    try:
        with open(base_path / splitT7_params_path, 'r', encoding='utf-8') as f:
            splitT7_params = json.load(f)
    except Exception:
        splitT7_params = {"alpha": 1.0, "Kd": 2.0e5, "leaky": 200.0}  # Kd 足够大避免早饱和
    return promoter_params, splitT7_params

# ------------------------------------------------------------------
# 简化 AND 门
# ------------------------------------------------------------------
class SimpleANDGate:
    def __init__(self, promoter_params=None, splitT7_params=None):
        if promoter_params is None or splitT7_params is None:
            promoter_params, splitT7_params = load_model_parameters()
        self.promoter_params = promoter_params
        self.splitT7_params = splitT7_params

    def _hill_function(self, x, params):
        # 支持标量/ndarray
        x = np.asarray(x, dtype=float)
        leaky = params['leaky']; beta = params['beta']; K = params['K']; n = params['n']
        ptype = params.get('type') or params.get('_mode') or params.get('mode') or 'act'
        ptype_l = str(ptype).lower()
        is_rep = (ptype_l.startswith('rep') or 'inh' in ptype_l)
        if not is_rep:  # 激活型
            return leaky + beta * (x**n) / (K**n + x**n)
        else:  # 抑制型 (x 增加输出降低)
            return leaky + beta * (K**n) / (K**n + x**n)

    def get_promoter_outputs(self, O2_percent, Temp_C):
        pPept_out = self._hill_function(O2_percent, self.promoter_params['pPept'])
        pLR_out   = self._hill_function(Temp_C,    self.promoter_params['pLR'])
        return pPept_out, pLR_out

    def get_t7_activity(self, O2_percent, Temp_C):
        A, B = self.get_promoter_outputs(O2_percent, Temp_C)
        # A,B 代表两片段表达量 (任意单位). 采用双底物样式饱和: 组装依赖乘积
        alpha = self.splitT7_params['alpha']
        Kd = self.splitT7_params['Kd']
        leaky = self.splitT7_params.get('leaky', 0.0)
        product = A * B
        return leaky + alpha * product / (Kd + product)

    # 调试辅助: 快速打印几个典型条件
    def quick_diagnose(self, O2_list=(1.0, 5.0, 21.0), Temp_list=(37.0, 42.0, 45.0)):
        print("\n[Quick Diagnose] Promoter parameters:")
        for k,v in self.promoter_params.items():
            print(f"  {k}: type={v.get('type')} K={v.get('K'):.3f} n={v.get('n'):.2f} beta={v.get('beta'):.1f} leaky={v.get('leaky'):.1f}")
        print("SplitT7:", self.splitT7_params)
        print("\nCondition Scan (T7 activity):")
        for o2 in O2_list:
            for T in Temp_list:
                t7 = self.get_t7_activity(o2, T)
                p1, p2 = self.get_promoter_outputs(o2, T)
                print(f"  O2={o2:5.1f}%  T={T:4.1f}°C  pPept={p1:8.1f}  pLR={p2:8.1f}  T7={t7:8.1f}")
        print("\n期望: 低 O2 + 高温 (例如 1% / 42°C 或以上) → T7 最高。若不满足，请调整 promoters.json 中 pLR 的 K (当前若高于目标温度会导致温度臂未开启)。")

# ------------------------------------------------------------------
# 详细模型 (占位, 保持接口)
# ------------------------------------------------------------------
class BaseDetailedANDGate:
    def __init__(self, **params):
        self.k_assembly = params.get('k_assembly', 1.0e-6)
        self.k_disassembly = params.get('k_disassembly', 1e-3)
        self.k_deg = params.get('k_deg', 0.05)
        self.simple = SimpleANDGate()

    def dydt(self, y, t, O2_percent, Temp_C):
        T7_active = y[0]
        A,B = self.simple.get_promoter_outputs(O2_percent, Temp_C)
        prod = A*B
        assemble = self.k_assembly * prod
        dis = self.k_disassembly * T7_active
        deg = self.k_deg * T7_active
        dT7 = assemble - dis - deg
        return [dT7]

    def simulate(self, O2_percent, Temp_C, t_end=24.0, dt=0.1):
        t = np.arange(0, t_end+dt, dt)
        y0 = [0.0]
        sol = odeint(self.dydt, y0, t, args=(O2_percent, Temp_C))
        return t, sol

class DetailedANDGate(BaseDetailedANDGate):
    def get_t7_activity(self, O2_percent, Temp_C, t_end=24.0):
        t, sol = self.simulate(O2_percent, Temp_C, t_end)
        return float(sol[-1,0])

# ------------------------------------------------------------------
# 自测
# ------------------------------------------------------------------
if __name__ == '__main__':
    gate = SimpleANDGate()
    gate.quick_diagnose()
    cases = [
        (1.0, 42.0, 'LOW O2 + HIGH T (ON)'),
        (21.0, 42.0, 'HIGH O2 + HIGH T (OFF via O2)'),
        (1.0, 37.0, 'LOW O2 + LOW T (OFF via Temp)'),
        (21.0, 37.0, 'HIGH O2 + LOW T (OFF both)')
    ]
    print('\n=== Simple AND Gate Test ===')
    for o2,temp,label in cases:
        act = gate.get_t7_activity(o2,temp)
        print(f'{label:30s} -> T7 = {act:.2f}')

    detailed = DetailedANDGate()
    print('\n=== Detailed (dynamic) AND Gate Test ===')
    for o2,temp,label in cases:
        act = detailed.get_t7_activity(o2,temp, t_end=12)
        print(f'{label:30s} -> T7_dynamic = {act:.2f}')
