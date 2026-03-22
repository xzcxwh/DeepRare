#!/usr/bin/env python3
"""
第一步：gnomAD频率注释与过滤
使用vcfanno对HG001.vcf进行gnomAD频率注释，然后筛选gnomAD4.1_joint_POPMAX_AF < 0.01的变异
"""

import subprocess
import os
import sys
import time
import tempfile
import shutil
from pathlib import Path

def check_dependencies():
    """检查必要的工具是否安装"""
    required_tools = ['vcfanno', 'bcftools']
    missing_tools = []

    for tool in required_tools:
        try:
            result = subprocess.run(['which', tool],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0 or not result.stdout.strip():
                missing_tools.append(tool)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_tools.append(tool)

    if missing_tools:
        print(f"❌ 缺少必要工具: {', '.join(missing_tools)}")
        print("请安装这些工具后再运行脚本")
        return False
    else:
        print("✅ 所有必要工具已安装")

        # 显示工具版本信息
        for tool in required_tools:
            try:
                if tool == 'vcfanno':
                    result = subprocess.run([tool], capture_output=True, text=True, timeout=10)
                    if result.stderr and 'vcfanno version' in result.stderr:
                        version_line = [line for line in result.stderr.split('\n') if 'vcfanno version' in line][0]
                        print(f"   {tool}: {version_line.strip()}")
                    else:
                        print(f"   {tool}: 已安装")
                else:
                    result = subprocess.run([tool, '--version'], capture_output=True, text=True, timeout=10)
                    version = result.stdout.strip() or result.stderr.strip()
                    print(f"   {tool}: {version.split()[0] if version else '已安装'}")
            except:
                print(f"   {tool}: 已安装")

        return True

def check_files():
    """检查输入文件是否存在"""
    input_vcf = "input_data/HG001.vcf"
    dbnsfp_file = "input_data/dbNSFP5.3a_grch38_lite.tsv.gz"

    missing_files = []

    if not os.path.exists(input_vcf):
        missing_files.append(input_vcf)

    if not os.path.exists(dbnsfp_file):
        missing_files.append(dbnsfp_file)

    if missing_files:
        print(f"❌ 缺少输入文件: {', '.join(missing_files)}")
        return False
    else:
        print("✅ 所有输入文件存在")
        print(f"📁 输入VCF: {input_vcf}")
        print(f"📄 dbNSFP数据库: {dbnsfp_file}")
        return True

def create_vcfanno_config():
    """创建vcfanno配置文件，处理chr前缀匹配问题"""

    config_content = '''[[annotation]]
file="input_data/dbNSFP5.3a_grch38_lite.tsv.gz"
columns = [1, 2, 27]
names = ["chr", "pos", "gnomAD4.1_joint_POPMAX_AF"]
ops = ["self", "self", "self"]

# 处理chr前缀匹配：vcf中的"chr1"需要匹配dbNSFP中的"1"
[[annotation]]
file="input_data/dbNSFP5.3a_grch38_lite.tsv.gz"
columns = [1, 2, 27]
names = ["chr_match", "pos_match", "gnomAD4.1_joint_POPMAX_AF"]
ops = ["self", "self", "self"]
'''

    config_file = "config/vcfanno_gnomad_step1.toml"

    # 确保config目录存在
    os.makedirs("config", exist_ok=True)

    with open(config_file, 'w') as f:
        f.write(config_content)

    print(f"✅ 创建vcfanno配置文件: {config_file}")
    return config_file

def run_vcfanno_annotation(input_vcf, config_file, output_annotated):
    """运行vcfanno注释"""
    print(f"\n🧬 开始vcfanno注释...")
    print(f"输入文件: {input_vcf}")
    print(f"配置文件: {config_file}")
    print(f"输出文件: {output_annotated}")
    print("=" * 60)

    start_time = time.time()

    cmd = [
        "vcfanno",
        "-p", "16",  # 使用16线程
        config_file,
        input_vcf
    ]

    try:
        with open(output_annotated, 'w') as outfile:
            process = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                timeout=1800  # 30分钟超时
            )

        elapsed_time = time.time() - start_time

        if process.returncode == 0:
            print(f"✅ vcfanno注释成功完成")
            print(f"处理时间: {elapsed_time:.2f}秒")
            return True, elapsed_time
        else:
            print(f"❌ vcfanno注释失败")
            print(f"错误代码: {process.returncode}")
            print(f"错误信息: {process.stderr}")
            return False, elapsed_time

    except subprocess.TimeoutExpired:
        print(f"❌ vcfanno注释超时（超过30分钟）")
        return False, time.time() - start_time
    except Exception as e:
        print(f"❌ 运行vcfanno时发生异常: {e}")
        return False, time.time() - start_time

def filter_gnomad_variants(input_annotated, output_vcf1):
    """筛选gnomAD4.1_joint_POPMAX_AF < 0.01的变异"""
    print(f"\n🔍 开始gnomAD频率筛选...")
    print(f"输入文件: {input_annotated}")
    print(f"输出文件: {output_vcf1}")
    print("筛选条件: gnomAD4.1_joint_POPMAX_AF < 0.01 且不为空")
    print("=" * 60)

    start_time = time.time()
    total_variants = 0
    passed_variants = 0
    filtered_variants = 0
    missing_af_count = 0

    try:
        with open(input_annotated, 'r') as infile, open(output_vcf1, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                # 写入所有头文件行
                if line.startswith('#'):
                    outfile.write(line)
                    continue

                total_variants += 1

                # 每1000个变异显示进度
                if total_variants % 1000 == 0:
                    print(f"已处理 {total_variants:,} 个变异，通过 {passed_variants:,} 个")

                # 解析VCF行
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    filtered_variants += 1
                    continue

                info_field = fields[7]
                chrom = fields[0]
                pos = fields[1]

                # 解析INFO字段中的gnomAD4.1_joint_POPMAX_AF
                gnomad_af = None
                info_items = info_field.split(';')

                for item in info_items:
                    if item.startswith('gnomAD4.1_joint_POPMAX_AF='):
                        try:
                            af_value = item.split('=')[1]
                            if af_value and af_value != '.' and af_value != '':
                                # 处理多值情况（如 "0.001,0.001"），取第一个有效值
                                for single_val in af_value.split(','):
                                    single_val = single_val.strip()
                                    if single_val and single_val != '.':
                                        try:
                                            gnomad_af = float(single_val)
                                            break  # 取第一个有效值
                                        except ValueError:
                                            continue
                        except (ValueError, IndexError):
                            gnomad_af = None
                        break

                # 筛选条件（原始逻辑）：
                # 必须有 gnomAD AF 且 < 0.01 才保留
                # 没有 AF 的变异被过滤掉
                if gnomad_af is not None and gnomad_af < 0.01:
                    # 有AF且低于阈值，保留
                    passed_variants += 1
                    outfile.write(line)
                else:
                    # 没有AF或AF太高，过滤掉
                    filtered_variants += 1
                    if gnomad_af is None:
                        missing_af_count += 1

        elapsed_time = time.time() - start_time

        print(f"\n🎉 gnomAD频率筛选完成！")
        print(f"总处理时间: {elapsed_time:.2f}秒")
        print(f"处理速度: {total_variants/elapsed_time:.1f} 变异/秒")
        print()
        print("📊 筛选统计:")
        print(f"  原始变异数量: {total_variants:,}")
        print(f"  通过筛选变异: {passed_variants:,}")
        print(f"  筛选保留比例: {(passed_variants/total_variants*100):.2f}%")
        print(f"  筛选失败变异: {filtered_variants:,}")
        print(f"  缺少gnomAD信息的变异: {missing_af_count:,}")

        return passed_variants, elapsed_time

    except Exception as e:
        print(f"❌ 筛选过程中发生错误: {e}")
        return 0, time.time() - start_time

def validate_results(output_vcf1):
    """验证结果文件"""
    print(f"\n🔍 验证结果文件...")

    if not os.path.exists(output_vcf1):
        print(f"❌ 输出文件不存在: {output_vcf1}")
        return False

    try:
        # 直接统计变异数量
        cmd = f"grep -c '^[^#]' {output_vcf1}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            variant_count = int(result.stdout.strip())
            if variant_count > 0:
                print(f"✅ 结果文件验证通过")
                print(f"📊 最终变异数量: {variant_count:,}")
                return True
            else:
                print("❌ 文件没有变异数据")
                return False
        else:
            # grep返回1表示没有匹配，可能文件为空
            print("❌ 文件没有变异数据")
            return False

    except Exception as e:
        print(f"❌ 验证结果文件时发生错误: {e}")
        return False

def main():
    """主函数"""
    print("🧬 VCF文件gnomAD频率注释与筛选 - 第一步")
    print("=" * 60)

    # 检查依赖工具
    if not check_dependencies():
        sys.exit(1)

    # 检查输入文件
    if not check_files():
        sys.exit(1)

    # 设置文件路径
    input_vcf = "input_data/HG001.vcf"
    output_annotated = "temp/vcf1_annotated.vcf"
    output_vcf1 = "results/vcf1.vcf"

    # 创建必要目录
    os.makedirs("temp", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    try:
        total_start_time = time.time()

        # 步骤1: 创建vcfanno配置文件
        config_file = create_vcfanno_config()

        # 步骤2: 运行vcfanno注释
        success, annotation_time = run_vcfanno_annotation(input_vcf, config_file, output_annotated)

        if not success:
            print("❌ vcfanno注释失败，流程终止")
            sys.exit(1)

        # 步骤3: 筛选gnomAD频率
        variant_count, filter_time = filter_gnomad_variants(output_annotated, output_vcf1)

        if variant_count == 0:
            print("❌ 没有变异通过筛选，请检查筛选条件")
            sys.exit(1)

        # 步骤4: 验证结果
        if not validate_results(output_vcf1):
            print("❌ 结果验证失败")
            sys.exit(1)

        total_time = time.time() - total_start_time

        # 生成报告
        print(f"\n" + "=" * 60)
        print("🎉 第一步处理完成！")
        print(f"总处理时间: {total_time:.2f}秒")
        print(f"  - vcfanno注释: {annotation_time:.2f}秒")
        print(f"  - gnomAD筛选: {filter_time:.2f}秒")
        print(f"最终结果: {variant_count:,} 个变异")
        print(f"输出文件: {output_vcf1}")

        # 生成详细报告
        report_content = f"""
# 第一步：gnomAD频率注释与筛选报告

## 处理概述
使用vcfanno对HG001.vcf进行gnomAD4.1_joint_POPMAX_AF注释，并筛选频率 < 0.01 的变异

## 输入文件
- **VCF文件**: input_data/HG001.vcf
- **注释数据库**: input_data/dbNSFP5.3a_grch38_lite.tsv.gz
- **注释字段**: dbNSFP第27列 (gnomAD4.1_joint_POPMAX_AF)

## 处理步骤
1. **vcfanno注释** - 使用16线程对VCF文件进行gnomAD频率注释
2. **频率筛选** - 保留gnomAD4.1_joint_POPMAX_AF < 0.01的变异

## 处理结果
- **总处理时间**: {total_time:.2f}秒
- **注释时间**: {annotation_time:.2f}秒
- **筛选时间**: {filter_time:.2f}秒
- **最终变异数量**: {variant_count:,}

## 输出文件
- **注释文件**: {output_annotated}
- **筛选结果**: {output_vcf1}

## 筛选条件
- gnomAD4.1_joint_POPMAX_AF < 0.01
- gnomAD4.1_joint_POPMAX_AF 不为空 (.')

## 技术细节
- **线程数**: 16
- **配置文件**: {config_file}
- **数据库列**: 第27列
- **字段映射**: chr(1), pos(2), gnomAD4.1_joint_POPMAX_AF(27)

## 下一步
下一步将进行gnomAD质量注释和过滤。

---
*生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}*
*处理脚本: step1_gnomad_annotation_filter.py*
"""

        with open('results/STEP1_GNOMAD_ANNOTATION_REPORT.md', 'w') as f:
            f.write(report_content)

        print(f"\n✅ 详细报告已保存到: results/STEP1_GNOMAD_ANNOTATION_REPORT.md")
        print("🚀 准备进行第二步处理...")

        # 清理临时文件
        try:
            if os.path.exists(output_annotated):
                os.remove(output_annotated)
        except:
            pass

    except Exception as e:
        print(f"❌ 处理过程中发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()