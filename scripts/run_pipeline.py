#!/usr/bin/env python3
"""
HERITA变异筛选流水线 - 自动化执行脚本
=====================================
按顺序执行 step1 ~ step8 的所有步骤

使用方法:
    python run_pipeline.py [input.vcf] [--clean]

参数:
    input.vcf  - 输入的VCF文件路径（可选，默认使用 input_data/HG001.vcf）
    --clean    - 运行前清理之前的结果文件

输出:
    results/final_integrated.vcf - 最终整合的变异结果文件
"""

import subprocess
import os
import sys
import shutil
import argparse
from datetime import datetime

# 获取项目根目录
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# 步骤脚本列表
STEPS = [
    ("step1_gnomad_annotation_filter.py", "gnomAD频率注释与过滤"),
    ("step2_gnomad_flag_annotation.py", "gnomAD Flag注释"),
    ("step3_sequencing_quality_filter.py", "测序质量过滤"),
    ("step4_clinvar_annotation.py", "ClinVar注释"),
    ("step5_dbscSNV_annotation.py", "dbscSNV注释"),
    ("step6_dbNSFP_annotation_filter.py", "dbNSFP多评分注释与过滤"),
    ("step7_final_top50_selection.py", "最终Top50筛选"),
    ("step8_final_integration.py", "最终整合"),
    ("step9_intervar_acmg.py", "InterVar ACMG分类"),
    ("step10_acmg_detailed_report.py", "生成详细ACMG分析日志"),
]

# 结果文件列表（用于清理）
RESULT_FILES = [
    "vcf1.vcf", "vcf2.vcf", "vcf3.vcf", "vcf4.vcf", "vcf5.vcf",
    "vcf5_annotated.vcf", "vcf6.vcf",
    "final_clinvar.vcf", "final_dbscSNV.vcf", "final_top.vcf",
    "merged_final.vcf", "merged_annotated.vcf", "merged_cleaned.vcf",
    "final_integrated.vcf", "dbscsnv_extracted.vcf",
    "final_acmg_classified.vcf", "final_intervar_classified.vcf",
    "STEP10_ACMG_DETAILED_LOG.txt"
]


def print_header():
    """打印流水线头部信息"""
    print("=" * 70)
    print("HERITA 变异筛选流水线")
    print("=" * 70)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"项目目录: {PROJECT_ROOT}")
    print("=" * 70)


def clean_results():
    """清理之前的结果文件"""
    print("\n📁 清理之前的结果文件...")
    results_dir = os.path.join(PROJECT_ROOT, "results")
    temp_dir = os.path.join(PROJECT_ROOT, "temp")
    
    cleaned_count = 0
    for filename in RESULT_FILES:
        filepath = os.path.join(results_dir, filename)
        if os.path.exists(filepath):
            os.remove(filepath)
            print(f"   删除: {filename}")
            cleaned_count += 1
    
    # 清理temp目录
    if os.path.exists(temp_dir):
        for f in os.listdir(temp_dir):
            fpath = os.path.join(temp_dir, f)
            if os.path.isfile(fpath):
                os.remove(fpath)
                cleaned_count += 1
    
    print(f"   共清理 {cleaned_count} 个文件")


def copy_input_vcf(input_vcf):
    """如果指定了输入VCF，复制到input_data目录"""
    if input_vcf:
        input_vcf = os.path.abspath(input_vcf)
        if not os.path.exists(input_vcf):
            print(f"❌ 错误: 输入文件不存在: {input_vcf}")
            sys.exit(1)
        
        # 目标路径
        target_vcf = os.path.join(PROJECT_ROOT, "input_data", "HG001.vcf")
        
        # 备份原文件（如果存在且不同）
        if os.path.exists(target_vcf) and input_vcf != target_vcf:
            backup_vcf = target_vcf + ".backup"
            if not os.path.exists(backup_vcf):
                shutil.copy2(target_vcf, backup_vcf)
                print(f"📦 已备份原输入文件到: HG001.vcf.backup")
        
        # 复制新文件
        if input_vcf != target_vcf:
            shutil.copy2(input_vcf, target_vcf)
            print(f"📥 已复制输入文件: {os.path.basename(input_vcf)} -> HG001.vcf")


def run_step(step_num, script_name, description, ignore_error=False):
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
        # 运行脚本
        result = subprocess.run(
            [sys.executable, script_path],
            cwd=PROJECT_ROOT,
            capture_output=False,  # 直接输出到终端
            text=True
        )
        
        elapsed = datetime.now() - start_time
        
        if result.returncode == 0:
            print(f"\n✅ Step {step_num} 完成 (耗时: {elapsed})")
            return True
        else:
            if ignore_error:
                print(f"\n⚠️ Step {step_num} 返回非零状态码 ({result.returncode})，但继续执行 (耗时: {elapsed})")
                return True
            else:
                print(f"\n❌ Step {step_num} 失败 (返回码: {result.returncode})")
                return False
            
    except Exception as e:
        print(f"\n❌ Step {step_num} 执行出错: {e}")
        return False


def count_variants(vcf_file):
    """统计VCF文件中的变异数量"""
    if not os.path.exists(vcf_file):
        return 0
    count = 0
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count


def print_summary():
    """打印流水线执行摘要"""
    results_dir = os.path.join(PROJECT_ROOT, "results")
    
    print("\n" + "=" * 70)
    print("流水线执行摘要")
    print("=" * 70)
    
    # 统计各阶段变异数
    files_to_check = [
        ("input_data/HG001.vcf", "输入文件"),
        ("results/vcf1.vcf", "Step 1 (gnomAD过滤后)"),
        ("results/vcf2.vcf", "Step 2 (gnomAD Flag后)"),
        ("results/vcf3.vcf", "Step 3 (质量过滤后)"),
        ("results/vcf4.vcf", "Step 4 (ClinVar注释后)"),
        ("results/vcf5.vcf", "Step 5 (dbscSNV注释后)"),
        ("results/vcf6.vcf", "Step 6 (评分过滤后)"),
        ("results/final_clinvar.vcf", "ClinVar致病变异"),
        ("results/final_dbscSNV.vcf", "dbscSNV高分变异"),
        ("results/final_top.vcf", "Top50高分变异"),
        ("results/final_integrated.vcf", "最终整合结果"),
        ("results/final_intervar_classified.vcf", "ACMG分类结果"),
    ]
    
    print(f"\n{'文件':<40} {'变异数':>15}")
    print("-" * 60)
    
    for filepath, description in files_to_check:
        full_path = os.path.join(PROJECT_ROOT, filepath)
        if os.path.exists(full_path):
            count = count_variants(full_path)
            print(f"{description:<40} {count:>15,}")
        else:
            print(f"{description:<40} {'(不存在)':>15}")
    
    # 显示最终结果位置
    final_vcf = os.path.join(results_dir, "final_intervar_classified.vcf")
    if os.path.exists(final_vcf):
        print("\n" + "=" * 70)
        print(f"✅ 最终结果: {final_vcf}")
        print("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description='HERITA变异筛选流水线 - 自动执行 step1~step9',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
    # 使用默认输入文件 (input_data/HG001.vcf)
    python run_pipeline.py
    
    # 指定输入VCF文件
    python run_pipeline.py /path/to/your.vcf
    
    # 清理之前的结果后运行
    python run_pipeline.py --clean
    
    # 指定输入文件并清理
    python run_pipeline.py /path/to/your.vcf --clean
    
    # 只运行部分步骤
    python run_pipeline.py --from-step 6 --to-step 9
        """
    )
    
    parser.add_argument('input_vcf', nargs='?', default=None,
                       help='输入VCF文件路径（可选，默认使用 input_data/HG001.vcf）')
    parser.add_argument('--clean', action='store_true',
                       help='运行前清理之前的结果文件')
    parser.add_argument('--from-step', type=int, default=1,
                       help='从指定步骤开始运行（默认从step1开始）')
    parser.add_argument('--to-step', type=int, default=10,
                       help='运行到指定步骤停止（默认到step10）')
    
    parser.add_argument('--ignore-errors', action='store_true',
                       help='忽略步骤执行失败继续运行')
    
    args = parser.parse_args()
    
    # 打印头部
    print_header()
    
    # 清理结果（如果指定）
    if args.clean:
        clean_results()
    
    # 复制输入文件（如果指定）
    if args.input_vcf:
        copy_input_vcf(args.input_vcf)
    
    # 检查输入文件是否存在
    input_vcf = os.path.join(PROJECT_ROOT, "input_data", "HG001.vcf")
    if not os.path.exists(input_vcf):
        print(f"\n❌ 错误: 输入文件不存在: {input_vcf}")
        print("请指定输入VCF文件: python run_pipeline.py /path/to/your.vcf")
        sys.exit(1)
    
    print(f"\n📄 输入文件: {input_vcf}")
    print(f"   变异数量: {count_variants(input_vcf):,}")
    
    # 验证步骤范围
    from_step = max(1, min(10, args.from_step))
    to_step = max(from_step, min(10, args.to_step))
    
    print(f"\n🚀 将执行步骤: Step {from_step} ~ Step {to_step}")
    
    # 记录开始时间
    pipeline_start = datetime.now()
    
    # 执行各步骤
    failed_step = None
    for i, (script_name, description) in enumerate(STEPS, 1):
        if i < from_step:
            continue
        if i > to_step:
            break
            
        success = run_step(i, script_name, description, ignore_error=args.ignore_errors)
        if not success:
            failed_step = i
            print(f"\n❌ 流水线在 Step {i} 失败，停止执行")
            break
    
    # 计算总耗时
    total_elapsed = datetime.now() - pipeline_start
    
    # 打印摘要
    print_summary()
    
    # 打印最终状态
    print(f"\n总耗时: {total_elapsed}")
    
    if failed_step:
        print(f"\n❌ 流水线执行失败（在 Step {failed_step}）")
        sys.exit(1)
    else:
        print(f"\n✅ 流水线执行完成！")
        sys.exit(0)


if __name__ == '__main__':
    main()
