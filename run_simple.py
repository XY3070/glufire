#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iGEM AND Gate Modeling Toolkit - Main Runner (Simple Version)
Simple runner script without complex encoding handling
"""

import sys
import subprocess
from pathlib import Path
import time

def run_script_simple(script_name, description):
    """Run a single script with simple error handling"""
    print(f"\n{'='*50}")
    print(f"Running: {description}")
    print(f"Script: {script_name}")
    print(f"{'='*50}")
    
    try:
        start_time = time.time()
        # Use system default encoding
        result = subprocess.run([sys.executable, script_name])
        end_time = time.time()
        
        if result.returncode == 0:
            print(f"✅ {description} - Success ({end_time-start_time:.1f}s)")
            return True
        else:
            print(f"❌ {description} - Failed (return code: {result.returncode})")
            return False
            
    except Exception as e:
        print(f"❌ {description} - Exception: {e}")
        return False

def main():
    print("""
    🧬 iGEM AND Gate Modeling Toolkit
    =================================
    
    Running analysis pipeline:
    1. Promoter transfer function fitting
    2. Split T7 AND gate modeling
    3. TX-TL to glutamate production simulation
    """)
    
    # Script execution order
    scripts = [
        ("01_promoter_fit.py", "Promoter Transfer Function Fitting"),
        ("02_splitT7_AND_model.py", "Split T7 AND Gate Modeling"),
        ("03_tx_tl_to_glu.py", "TX-TL Glutamate Production Simulation")
    ]
    
    success_count = 0
    
    for script, description in scripts:
        if Path(script).exists():
            if run_script_simple(script, description):
                success_count += 1
            else:
                print(f"Script {script} failed, continuing with next...")
        else:
            print(f"⚠️  Script {script} not found, skipping")
    
    print(f"\n{'='*50}")
    print(f"Summary: {success_count}/{len(scripts)} scripts completed successfully")
    
    if success_count == len(scripts):
        print("🎉 All analyses completed!")
        print("\nGenerated files:")
        print("📁 params/ - Model parameter files")
        print("🖼️  *.png - Analysis plots")
    else:
        print("⚠️  Some scripts failed, please check error messages")
    
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
