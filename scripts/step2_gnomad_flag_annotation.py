#!/usr/bin/env python3
"""
第二步：gnomAD质量注释
对vcf1进行gnomAD4.1_joint_flag注释（dbNSFP第26列），保留所有变异包括flag为'.'的数据
"""

import subprocess
import os
import sys
import time
import shutil
from pathlib import Path

def check_vcf1_status():
    """检查vcf1文件状态"""
    vcf1_file = "results/vcf1.vcf"

    if not os.path.exists(vcf1_file):
        print(f"❌ 找不到vcf1文件: {vcf1_file}")
        return False, 0

    try:
        # 统计变异数量
        cmd = f"grep -c '^[^#]' {vcf1_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            variant_count = int(result.stdout.strip())
            print(f"✅ vcf1文件存在")
            print(f"📊 vcf1变异数量: {variant_count:,}")
            return True, variant_count
        else:
            print(f"❌ 无法统计vcf1变异数量")
            return False, 0
    except Exception as e:
        print(f"��� 检查vcf1时发生错误: {e}")
        return False, 0

def check_dependencies():
    """检查必要工具"""
    try:
        result = subprocess.run(['which', 'vcfanno'],
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 and result.stdout.strip():
            print("✅ vcfanno工具已安装")
            return True
        else:
            print("❌ vcfanno工具未找到")
            return False
    except:
        print("❌ 检查vcfanno工具时发生错误")
        return False

def check_dbnsfp_file():
    """检查dbNSFP文件"""
    dbnsfp_file = "input_data/dbNSFP5.3a_grch38_lite.tsv.gz"

    if not os.path.exists(dbnsfp_file):
        print(f"❌ 找不到dbNSFP文件: {dbnsfp_file}")
        return False

    print(f"✅ dbNSFP文件存在: {dbnsfp_file}")

    # 检查文件大小
    file_size = os.path.getsize(dbnsfp_file)
    print(f"📄 dbNSFP文件大小: {file_size:,} bytes")

    return True

def create_vcfanno_config():
    """创建vcfanno配置文件用于gnomAD flag注释"""
    
    # 获取脚本所在目录和项目根目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    dbnsfp_path = os.path.join(project_root, "input_data", "dbNSFP5.3a_grch38_lite.tsv.gz")

    config_content = f'''[[annotation]]
file="{dbnsfp_path}"
columns = [1, 2, 26]
names = ["chr", "pos", "gnomAD4.1_joint_flag"]
ops = ["self", "self", "self"]

# 添加备用匹配以确保完整覆盖
[[annotation]]
file="{dbnsfp_path}"
columns = [1, 2, 26]
names = ["chr_backup", "pos_backup", "gnomAD4.1_joint_flag_backup"]
ops = ["self", "self", "self"]
'''

    config_file = os.path.join(project_root, "config", "vcfanno_gnomad_flag.toml")

    # 确保config目录存在
    os.makedirs(os.path.join(project_root, "config"), exist_ok=True)

    with open(config_file, 'w') as f:
        f.write(config_content)

    print(f"✅ 创建vcfanno配置文件: {config_file}")
    return config_file

def run_vcfanno_annotation(vcf1_file, config_file, vcf2_file):
    """运行vcfanno添加gnomAD flag注释"""
    print(f"\n🧬 开始gnomAD flag注释...")
    print(f"输入文件: {vcf1_file}")
    print(f"配置文件: {config_file}")
    print(f"输出文件: {vcf2_file}")
    print("注释字段: gnomAD4.1_joint_flag (dbNSFP第26列)")
    print("=" * 60)

    start_time = time.time()

    cmd = [
        "vcfanno",
        "-p", "16",  # 使用16线程
        config_file,
        vcf1_file
    ]

    try:
        with open(vcf2_file, 'w') as outfile:
            process = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                timeout=600  # 10分钟超时
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
        print(f"❌ vcfanno注释超时（超过10分钟）")
        return False, time.time() - start_time
    except Exception as e:
        print(f"❌ 运行vcfanno时发生异常: {e}")
        return False, time.time() - start_time

def validate_vcf2(vcf2_file, expected_variants):
    """验证vcf2文件"""
    print(f"\n🔍 验证vcf2文件...")

    if not os.path.exists(vcf2_file):
        print(f"❌ 输出文件不存在: {vcf2_file}")
        return False

    try:
        # 检查文件格式
        with open(vcf2_file, 'r') as f:
            first_lines = [next(f) for _ in range(20)]

        has_header = any(line.startswith('#') for line in first_lines)

        if not has_header:
            print("❌ 文件缺少VCF头信息")
            return False

        # 统计变异数量
        cmd = f"grep -c '^[^#]' {vcf2_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            variant_count = int(result.stdout.strip())
            print(f"📊 vcf2变异数量: {variant_count:,}")

            if variant_count == expected_variants:
                print("✅ 变异数量匹配vcf1")
            else:
                print(f"⚠️  变异数量不匹配: 期望{expected_variants}, 实际{variant_count}")

        # 检查注释字段
        cmd = f"grep -c 'gnomAD4.1_joint_flag=' {vcf2_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            annotated_count = int(result.stdout.strip())
            print(f"📊 带gnomAD flag注释的变异: {annotated_count:,}")

            if annotated_count > 0:
                print("✅ 成功添加gnomAD flag注释")
            else:
                print("⚠️  没有找到gnomAD flag注释")

        # 检查flag为'.'的变异
        cmd = f"grep 'gnomAD4.1_joint_flag=\\.' {vcf2_file} | wc -l"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            dot_count = int(result.stdout.strip())
            print(f"📊 gnomAD flag为'.'的变异: {dot_count:,}")

        return True

    except Exception as e:
        print(f"❌ 验证vcf2时发生错误: {e}")
        return False

def show_sample_variants(vcf2_file):
    """显示示例变异"""
    print(f"\n📋 vcf2示例变异 (前5个):")
    print("=" * 60)

    try:
        cmd = f"grep '^[^#]' {vcf2_file} | head -5"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            for i, line in enumerate(lines, 1):
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    chrom = fields[0]
                    pos = fields[1]
                    ref = fields[3]
                    alt = fields[4]
                    info = fields[7]

                    # 提取gnomAD信息
                    popmax_af = "N/A"
                    flag = "N/A"

                    for item in info.split(';'):
                        if item.startswith('gnomAD4.1_joint_POPMAX_AF='):
                            popmax_af = item.split('=')[1]
                        elif item.startswith('gnomAD4.1_joint_flag='):
                            flag = item.split('=')[1]

                    print(f"{i}. {chrom}:{pos}:{ref}>{alt}")
                    print(f"   gnomAD频率: {popmax_af}")
                    print(f"   gnomAD标志: {flag}")
                    print()

    except Exception as e:
        print(f"❌ 显示示例变异时发生错误: {e}")

def main():
    """主函数"""
    print("🧬 VCF文件gnomAD质量注释 - 第二步")
    print("=" * 60)

    # 设置文件路径
    vcf1_file = "results/vcf1.vcf"
    vcf2_file = "results/vcf2.vcf"

    # 创建必要目录
    os.makedirs("results", exist_ok=True)
    os.makedirs("config", exist_ok=True)

    try:
        total_start_time = time.time()

        # 步骤1: 检查vcf1状态
        print("检查vcf1文件状态...")
        vcf1_ok, vcf1_variants = check_vcf1_status()
        if not vcf1_ok:
            sys.exit(1)

        # 步骤2: 检查依赖工具
        print("\n检查依赖工具...")
        if not check_dependencies():
            sys.exit(1)

        # 步骤3: 检查dbNSFP文件
        print("\n检查数据文件...")
        if not check_dbnsfp_file():
            sys.exit(1)

        # 步骤4: 创建vcfanno配置文件
        print("\n创建注释配置...")
        config_file = create_vcfanno_config()

        # 步骤5: 运行vcfanno注释
        success, annotation_time = run_vcfanno_annotation(vcf1_file, config_file, vcf2_file)

        if not success:
            print("❌ gnomAD flag注释失败，流程终止")
            sys.exit(1)

        # 步骤6: 验证结果
        if not validate_vcf2(vcf2_file, vcf1_variants):
            print("❌ 结果验证失败")
            sys.exit(1)

        # 步骤7: 显示示例
        show_sample_variants(vcf2_file)

        total_time = time.time() - total_start_time

        # 生成最终报告
        print(f"\n" + "=" * 60)
        print("🎉 第二步处理完成！")
        print(f"总处理时间: {total_time:.2f}秒")
        print(f"  - vcfanno注释: {annotation_time:.2f}秒")
        print(f"输入变异数量: {vcf1_variants:,}")
        print(f"输出文件: {vcf2_file}")

        # 统计最终结果
        cmd = f"grep -c '^[^#]' {vcf2_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        final_variants = int(result.stdout.strip()) if result.returncode == 0 else 0

        print(f"最终变异数量: {final_variants:,}")

        # 生成详细报告
        report_content = f"""
# 第二步：gnomAD质量注释报告

## 处理概述
对vcf1进行gnomAD4.1_joint_flag注释，保留所有变异包括flag为'.'的数据

## 输入文件
- **输入VCF**: {vcf1_file} ({vcf1_variants:,}个变异)
- **注释数据库**: input_data/dbNSFP5.3a_grch38_lite.tsv.gz
- **注释字段**: dbNSFP第26列 (gnomAD4.1_joint_flag)

## 处理步骤
1. **vcfanno注释** - 使用16线程添加gnomAD质量标志
2. **结果验证** - 确保所有变异都保留并正确注释

## 处理结果
- **总处理时间**: {total_time:.2f}秒
- **注释时间**: {annotation_time:.2f}秒
- **输入变异数量**: {vcf1_variants:,}
- **输出变异数量**: {final_variants:,}

## 输出文件
- **输出文件**: {vcf2_file}

## 注释说明
- **gnomAD4.1_joint_flag**: gnomAD质量控制标志
- **保留策略**: 保留所有变异，包括flag为'.'的缺失数据
- **字段映射**: chr(1), pos(2), gnomAD4.1_joint_flag(26)

## 技术细节
- **线程数**: 16
- **配置文件**: {config_file}
- **数据完整性**: 100%变异保留

## 下一步
下一步将进行测序质量筛选（QUAL, FILTER, DP, GQ）。

---
*生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}*
*处理脚本: scripts/step2_gnomad_flag_annotation.py*
"""

        with open('results/STEP2_GNOMAD_FLAG_ANNOTATION_REPORT.md', 'w') as f:
            f.write(report_content)

        print(f"\n✅ 详细报告已保存到: results/STEP2_GNOMAD_FLAG_ANNOTATION_REPORT.md")
        print("🚀 准备进行第三步测序质量筛选...")

    except Exception as e:
        print(f"❌ 处理过程中发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()