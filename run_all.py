#!/usr/bin/env python3
"""
iGEM AND门谷氨酸生产系统建模工具包
主运行脚本

使用方法:
python run_all.py [--skip-plots] [--quick]

选项:
--skip-plots : 跳过图像生成
--quick      : 快速运行模式（减少数据点）
"""

import sys
import argparse
import subprocess
from pathlib import Path
import time

def run_script(script_name, description):
    """运行单个脚本并处理错误"""
    print(f"\n{'='*50}")
    print(f"运行: {description}")
    print(f"脚本: {script_name}")
    print(f"{'='*50}")
    
    try:
        start_time = time.time()
        # 使用utf-8编码并忽略编码错误
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=True, text=True, 
                              encoding='utf-8', errors='replace')
        end_time = time.time()
        
        if result.returncode == 0:
            print(f"✅ {description} - Success ({end_time-start_time:.1f}s)")
            if result.stdout:
                print("Output:")
                print(result.stdout)
        else:
            print(f"❌ {description} - Failed")
            print("Error:")
            if result.stderr:
                print(result.stderr)
            else:
                print("No error message available")
            return False
            
    except Exception as e:
        print(f"❌ {description} - Exception: {e}")
        return False
    
    return True

def check_dependencies():
    """Check required packages"""
    print("Checking dependencies...")
    required_packages = ['numpy', 'scipy', 'matplotlib', 'pandas', 'lmfit', 'sklearn']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"✅ {package}")
        except ImportError:
            print(f"❌ {package} - Not installed")
            missing_packages.append(package)
    
    if missing_packages:
        print(f"\nMissing packages: {', '.join(missing_packages)}")
        print("Please run: pip install -r requirements.txt")
        return False
    
    print("All dependencies installed ✅")
    return True

def main():
    parser = argparse.ArgumentParser(description='运行iGEM AND门建模工具包')
    parser.add_argument('--skip-plots', action='store_true', help='跳过图像生成')
    parser.add_argument('--quick', action='store_true', help='快速运行模式')
    parser.add_argument('--check-deps', action='store_true', help='只检查依赖包')
    
    args = parser.parse_args()
    
    # 检查依赖
    if not check_dependencies():
        sys.exit(1)
    
    if args.check_deps:
        print("依赖检查完成")
        return
    
    print("""
    🧬 iGEM AND门 → 谷氨酸 → 铁死亡建模工具包
    ==========================================
    
    这个工具包将依次运行以下分析:
    1. 启动子传递函数拟合 (pPept vs O2, pL/pR vs 温度)
    2. 分裂T7 AND门建模
    3. TX-TL到谷氨酸生产动力学模拟
    
    注意: 确保data/目录中有相应的CSV数据文件
    """)
    
    # 脚本运行顺序
    scripts = [
        ("01_promoter_fit.py", "启动子传递函数拟合"),
        ("02_splitT7_AND_model.py", "分裂T7 AND门建模"),
        ("03_tx_tl_to_glu.py", "TX-TL谷氨酸生产模拟")
    ]
    
    success_count = 0
    
    for script, description in scripts:
        if Path(script).exists():
            if run_script(script, description):
                success_count += 1
            else:
                print(f"脚本 {script} 运行失败，继续运行下一个...")
        else:
            print(f"⚠️  脚本 {script} 不存在，跳过")
    
    print(f"\n{'='*50}")
    print(f"运行总结: {success_count}/{len(scripts)} 个脚本成功完成")
    
    if success_count == len(scripts):
        print("🎉 所有分析完成！")
        print("\n生成的文件:")
        print("📁 params/ - 拟合参数文件")
        print("🖼️  *.png - 分析图像")
    else:
        print("⚠️  部分脚本运行失败，请检查错误信息")
    
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
