#!/usr/bin/env python3
"""
Indel Pipeline - 主运行脚本
执行 Indel 变异的注释与筛选流水线
"""

import os
import sys
import subprocess
from datetime import datetime

# 项目路径
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Indel 流水线步骤
STEPS = [
    ("step0_extract_indel.py", "提取 Indel 变异"),
    ("step1_gnomad_annotation.py", "gnomAD 注释与频率筛选"),
]


def run_step(step_num, script_name, description):
    """运行单个步骤"""
    script_path = os.path.join(SCRIPT_DIR, script_name)
    
    if not os.path.exists(script_path):
        print(f"❌ 错误: 脚本不存在: {script_path}")
        return False
    
    print(f"\n{'='*70}")
    print(f"Step {step_num}: {description}")
    print(f"脚本: {script_name}")
    print("=" * 70)
    
    start_time = datetime.now()
    
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            cwd=PROJECT_DIR,
            capture_output=False
        )
        
        elapsed = datetime.now() - start_time
        
        if result.returncode == 0:
            print(f"\n✅ Step {step_num} 完成 (耗时: {elapsed})")
            return True
        else:
            print(f"\n❌ Step {step_num} 失败")
            return False
            
    except Exception as e:
        print(f"\n❌ Step {step_num} 执行出错: {e}")
        return False


def count_variants(vcf_file):
    """统计 VCF 文件中的变异数量"""
    if not os.path.exists(vcf_file):
        return 0
    count = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count


def main():
    print("=" * 70)
    print("🧬 Indel 变异筛选流水线")
    print("=" * 70)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"项目目录: {PROJECT_DIR}")
    print("=" * 70)
    
    # 检查输入文件
    input_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001.vcf")
    if not os.path.exists(input_vcf):
        print(f"\n❌ 错误: 输入文件不存在: {input_vcf}")
        sys.exit(1)
    
    print(f"\n📄 输入文件: {input_vcf}")
    print(f"   总变异数: {count_variants(input_vcf):,}")
    
    # 创建结果目录
    results_dir = os.path.join(PROJECT_DIR, "results", "indel")
    os.makedirs(results_dir, exist_ok=True)
    
    pipeline_start = datetime.now()
    
    # 执行各步骤
    for i, (script_name, description) in enumerate(STEPS):
        success = run_step(i, script_name, description)
        if not success:
            print(f"\n❌ 流水线在 Step {i} 失败")
            sys.exit(1)
    
    total_elapsed = datetime.now() - pipeline_start
    
    # 打印摘要
    print("\n" + "=" * 70)
    print("流水线执行摘要")
    print("=" * 70)
    
    files_to_check = [
        ("input_data/HG001.vcf", "输入文件"),
        ("results/indel/indel_raw.vcf", "Step 0: Indel 提取"),
        ("results/indel/indel_gnomad_filtered.vcf", "Step 1: gnomAD 筛选后"),
    ]
    
    print(f"\n{'文件':<45} {'变异数':>12}")
    print("-" * 60)
    
    for filepath, description in files_to_check:
        full_path = os.path.join(PROJECT_DIR, filepath)
        if os.path.exists(full_path):
            count = count_variants(full_path)
            print(f"{description:<45} {count:>12,}")
        else:
            print(f"{description:<45} {'(不存在)':>12}")
    
    print("\n" + "=" * 70)
    print(f"总耗时: {total_elapsed}")
    print("=" * 70)
    print("\n✅ Indel 流水线执行完成!")


if __name__ == "__main__":
    main()
