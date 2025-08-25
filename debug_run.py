#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Debug version of the main script with ASCII-only output
"""

import sys
import subprocess
from pathlib import Path
import time

def run_script_debug(script_name):
    """Run script with debug output"""
    print(f"\n{'='*50}")
    print(f"Running: {script_name}")
    print(f"{'='*50}")
    
    try:
        start_time = time.time()
        # Run without capturing output to see real-time results
        result = subprocess.run([sys.executable, script_name])
        end_time = time.time()
        
        print(f"\nScript finished with return code: {result.returncode}")
        print(f"Execution time: {end_time-start_time:.1f}s")
        return result.returncode == 0
        
    except Exception as e:
        print(f"Exception: {e}")
        return False

def main():
    print("iGEM AND Gate Modeling Toolkit - Debug Version")
    print("=" * 50)
    
    scripts = [
        "01_promoter_fit.py",
        "02_splitT7_AND_model.py", 
        "03_tx_tl_to_glu.py"
    ]
    
    for script in scripts:
        if Path(script).exists():
            print(f"\nFound script: {script}")
            success = run_script_debug(script)
            if success:
                print(f"SUCCESS: {script}")
            else:
                print(f"FAILED: {script}")
        else:
            print(f"NOT FOUND: {script}")

if __name__ == "__main__":
    main()
