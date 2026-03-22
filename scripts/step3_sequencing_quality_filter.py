#!/usr/bin/env python3
"""
第三步：测序质量筛选
对vcf2进行测序质量筛选，保留满足以下条件的变异：
- QUAL >= 30
- FILTER = PASS
- DP >= 30
- GQ >= 60
"""

import subprocess
import os
import sys
import time
import re
from pathlib import Path

def check_vcf2_status():
    """检查vcf2文件状态"""
    vcf2_file = "results/vcf2.vcf"

    if not os.path.exists(vcf2_file):
        print(f"❌ 找不到vcf2文件: {vcf2_file}")
        return False, 0

    try:
        # 统计变异数量
        cmd = f"grep -c '^[^#]' {vcf2_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            variant_count = int(result.stdout.strip())
            print(f"✅ vcf2文件存在")
            print(f"📊 vcf2变异数量: {variant_count:,}")
            return True, variant_count
        else:
            print(f"❌ 无法统计vcf2变异数量")
            return False, 0
    except Exception as e:
        print(f"❌ 检查vcf2时发生错误: {e}")
        return False, 0

def analyze_quality_metrics(vcf2_file):
    """分析vcf2中的质量指标分布"""
    print(f"\n📊 分析vcf2质量指标分布...")
    print("=" * 50)

    try:
        # 使用awk分析质量指标
        cmd = f"""awk '
        BEGIN {{
            count = 0
            qual_count = 0
            filter_count = 0
            dp_count = 0
            gq_count = 0
            pass_count = 0
        }}
        $0 !~ /^#/ {{
            count++
            qual = $5
            filter_status = $6
            format_str = $8
            sample_str = $9

            # 统计QUAL分布
            if (qual != ".") {{
                qual_count++
            }}

            # 统计FILTER分布
            if (filter_status != ".") {{
                filter_count++
                if (filter_status == "PASS") {{
                    pass_count++
                }}
            }}

            # 解析FORMAT字段
            split(format_str, format_fields, ":")
            split(sample_str, sample_fields, ":")

            # 获取DP和GQ索引
            dp_idx = 0
            gq_idx = 0
            for (i = 1; i <= length(format_fields); i++) {{
                if (format_fields[i] == "DP") dp_idx = i
                if (format_fields[i] == "GQ") gq_idx = i
            }}

            # 统计DP
            if (length(sample_fields) > dp_idx && sample_fields[dp_idx+1] != ".") {{
                dp_count++
            }}

            # 统计GQ
            if (length(sample_fields) > gq_idx && sample_fields[gq_idx+1] != ".") {{
                gq_count++
            }}
        }}
        END {{
            printf "变异总数: %d\\n", count
            printf "有QUAL信息的变异: %d\\n", qual_count
            printf "有FILTER信息的变异: %d\\n", filter_count
            printf "FILTER=PASS的变异: %d\\n", pass_count
            printf "有DP信息的变异: %d\\n", dp_count
            printf "有GQ信息的变异: %d\\n", gq_count
        }}
        ' {vcf2_file}"""

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(result.stdout)
        else:
            print("⚠️  无法分析质量指标分布")

    except Exception as e:
        print(f"⚠️  分析质量指标时发生错误: {e}")

def quality_filter_variants(vcf2_file, vcf3_file):
    """执行测序质量筛选"""
    print(f"\n🧬 开始测序质量筛选...")
    print(f"输入文件: {vcf2_file}")
    print(f"输出文件: {vcf3_file}")
    print("筛选条件:")
    print("  - QUAL >= 30")
    print("  - FILTER = PASS")
    print("  - DP >= 30")
    print("  - GQ >= 60")
    print("=" * 50)

    start_time = time.time()

    # 筛选统计
    total_variants = 0
    passed_variants = 0
    filtered_variants = 0

    # 失败原因统计
    qual_fail = 0
    filter_fail = 0
    dp_fail = 0
    gq_fail = 0
    missing_info = 0

    try:
        with open(vcf2_file, 'r') as infile, open(vcf3_file, 'w') as outfile:
            for line in infile:
                # 写入所有头文件行
                if line.startswith('#'):
                    outfile.write(line)
                    continue

                total_variants += 1

                # 每50个变异显示进度
                if total_variants % 50 == 0:
                    print(f"已处理 {total_variants:,} 个变异，通过 {passed_variants:,} 个")

                # 分割VCF行
                fields = line.strip().split('\t')
                if len(fields) < 10:
                    missing_info += 1
                    filtered_variants += 1
                    continue

                # 提取质量信息
                chrom = fields[0]
                pos = fields[1]
                qual_str = fields[5]
                filter_status = fields[6]
                format_str = fields[8]
                sample_str = fields[9]

                # 筛选1: QUAL >= 30
                try:
                    qual = float(qual_str)
                except ValueError:
                    qual = 0.0

                if qual < 30.0:
                    qual_fail += 1
                    filtered_variants += 1
                    continue

                # 筛选2: FILTER = PASS
                if filter_status != "PASS":
                    filter_fail += 1
                    filtered_variants += 1
                    continue

                # 解析FORMAT字段
                format_fields = format_str.split(':')
                sample_fields = sample_str.split(':')

                # 获取DP和GQ的索引
                dp_index = None
                gq_index = None

                for i, format_field in enumerate(format_fields):
                    if format_field == "DP":
                        dp_index = i
                    elif format_field == "GQ":
                        gq_index = i

                # 筛选3: DP >= 30
                if dp_index is not None and len(sample_fields) > dp_index:
                    try:
                        dp = int(sample_fields[dp_index])
                        if dp < 30:
                            dp_fail += 1
                            filtered_variants += 1
                            continue
                    except (ValueError, IndexError):
                        dp_fail += 1
                        filtered_variants += 1
                        continue
                else:
                    dp_fail += 1
                    filtered_variants += 1
                    continue

                # 筛选4: GQ >= 60
                if gq_index is not None and len(sample_fields) > gq_index:
                    try:
                        gq = int(sample_fields[gq_index])
                        if gq < 60:
                            gq_fail += 1
                            filtered_variants += 1
                            continue
                    except (ValueError, IndexError):
                        gq_fail += 1
                        filtered_variants += 1
                        continue
                else:
                    gq_fail += 1
                    filtered_variants += 1
                    continue

                # 所有筛选条件都满足，写入输出文件
                outfile.write(line)
                passed_variants += 1

        elapsed_time = time.time() - start_time

        print(f"\n🎉 测序质量筛选完成！")
        print(f"总处理时间: {elapsed_time:.2f}秒")
        print(f"处理速度: {total_variants/elapsed_time:.1f} 变异/秒")
        print()
        print("📊 筛选统计:")
        print(f"  原始变异数量: {total_variants:,}")
        print(f"  通过筛选变异: {passed_variants:,}")
        print(f"  筛选保留比例: {(passed_variants/total_variants*100):.2f}%")
        print(f"  筛选失败变异: {filtered_variants:,}")
        print()
        print("📊 筛选失败原因:")
        print(f"  QUAL < 30: {qual_fail:,}")
        print(f"  FILTER ≠ PASS: {filter_fail:,}")
        print(f"  DP < 30: {dp_fail:,}")
        print(f"  GQ < 60: {gq_fail:,}")
        print(f"  信息缺失: {missing_info:,}")

        return passed_variants, elapsed_time, {
            'qual_fail': qual_fail,
            'filter_fail': filter_fail,
            'dp_fail': dp_fail,
            'gq_fail': gq_fail,
            'missing_info': missing_info
        }

    except Exception as e:
        print(f"❌ 筛选过程中发生错误: {e}")
        return 0, time.time() - start_time, {}

def validate_vcf3(vcf3_file, expected_min=0):
    """验证vcf3文件"""
    print(f"\n🔍 验证vcf3文件...")

    if not os.path.exists(vcf3_file):
        print(f"❌ 输出文件不存在: {vcf3_file}")
        return False, 0

    try:
        # 直接统计变异数量
        cmd = f"grep -c '^[^#]' {vcf3_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            variant_count = int(result.stdout.strip())
            if variant_count > 0:
                print(f"✅ 结果文件验证通过")
                print(f"📊 vcf3变异数量: {variant_count:,}")

                if variant_count >= expected_min:
                    print(f"✅ 变异数量符合预期 (≥ {expected_min})")
                else:
                    print(f"⚠️  变异数量少于预期")

                return True, variant_count
            else:
                print("❌ 文件没有变异数据")
                return False, 0
        else:
            print("❌ 文件没有变异数据")
            return False, 0

    except Exception as e:
        print(f"❌ 验证vcf3时发生错误: {e}")
        return False, 0

def show_sample_variants(vcf3_file):
    """显示通过筛选的示例变异"""
    print(f"\n📋 vcf3示例变异 (前5个通过筛选的变异):")
    print("=" * 60)

    try:
        cmd = f"grep '^[^#]' {vcf3_file} | head -5"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            lines = result.stdout.strip().split('\n')
            for i, line in enumerate(lines, 1):
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        chrom = fields[0]
                        pos = fields[1]
                        ref = fields[3]
                        alt = fields[4]
                        qual = fields[5]
                        filter_status = fields[6]
                        format_str = fields[8]
                        sample_str = fields[9]

                        # 解析DP和GQ
                        dp = "N/A"
                        gq = "N/A"

                        format_fields = format_str.split(':')
                        sample_fields = sample_str.split(':')

                        for j, format_field in enumerate(format_fields):
                            if format_field == "DP" and j < len(sample_fields):
                                dp = sample_fields[j]
                            elif format_field == "GQ" and j < len(sample_fields):
                                gq = sample_fields[j]

                        print(f"{i}. {chrom}:{pos}:{ref}>{alt}")
                        print(f"   QUAL: {qual}")
                        print(f"   FILTER: {filter_status}")
                        print(f"   DP: {dp}")
                        print(f"   GQ: {gq}")
                        print()

    except Exception as e:
        print(f"❌ 显示示例变异时发生错误: {e}")

def main():
    """主函数"""
    print("🧬 VCF文件测序质量筛选 - 第三步")
    print("=" * 60)

    # 设置文件路径
    vcf2_file = "results/vcf2.vcf"
    vcf3_file = "results/vcf3.vcf"

    # 创建必要目录
    os.makedirs("results", exist_ok=True)

    try:
        total_start_time = time.time()

        # 步骤1: 检查vcf2状态
        print("检查vcf2文件状态...")
        vcf2_ok, vcf2_variants = check_vcf2_status()
        if not vcf2_ok:
            sys.exit(1)

        # 步骤2: 分析质量指标
        analyze_quality_metrics(vcf2_file)

        # 步骤3: 执行质量筛选
        passed_variants, filter_time, fail_stats = quality_filter_variants(vcf2_file, vcf3_file)

        if passed_variants == 0:
            print("❌ 没有变异通过质量筛选")
            sys.exit(1)

        # 步骤4: 验证结果
        success, vcf3_variants = validate_vcf3(vcf3_file, 0)

        if not success:
            print("❌ 结果验证失败")
            sys.exit(1)

        # 步骤5: 显示示例
        show_sample_variants(vcf3_file)

        total_time = time.time() - total_start_time

        # 生成最终报告
        print(f"\n" + "=" * 60)
        print("🎉 第三步处理完成！")
        print(f"总处理时间: {total_time:.2f}秒")
        print(f"  - 质量指标分析: 完成")
        print(f"  - 质量筛选: {filter_time:.2f}秒")
        print(f"输入变异数量: {vcf2_variants:,}")
        print(f"输出变异数量: {vcf3_variants:,}")
        print(f"筛选保留比例: {(vcf3_variants/vcf2_variants*100):.2f}%")
        print(f"输出文件: {vcf3_file}")

        # 生成详细报告
        report_content = f"""
# 第三步：测序质量筛选报告

## 处理概述
对vcf2进行严格的测序质量筛选，保留满足以下条件的高质量变异：
- QUAL ≥ 30 (变异质量分数)
- FILTER = PASS (通过所有质量过滤)
- DP ≥ 30 (测序深度)
- GQ ≥ 60 (基因型质量分数)

## 输入文件
- **输入VCF**: {vcf2_file} ({vcf2_variants:,}个变异)
- **输出VCF**: {vcf3_file} ({vcf3_variants:,}个变异)

## 筛选标准
1. **变异质量分数 (QUAL)**: 衡量变异检测的可信度
2. **过滤状态 (FILTER)**: 变异是否通过质量控制
3. **测序深度 (DP)**: 支持变异的reads数量
4. **基因型质量 (GQ)**: 基因型判定的可信度

## 处理结果
- **总处理时间**: {total_time:.2f}秒
- **质量筛选时间**: {filter_time:.2f}秒
- **输入变异数量**: {vcf2_variants:,}
- **输出变异数量**: {vcf3_variants:,}
- **筛选保留比例**: {(vcf3_variants/vcf2_variants*100):.2f}%

## 筛选失败原因分析
| 筛选条件 | 失败数量 | 失败比例 |
|----------|----------|----------|
| QUAL < 30 | {fail_stats.get('qual_fail', 0):,} | {(fail_stats.get('qual_fail', 0)/vcf2_variants*100):.2f}% |
| FILTER ≠ PASS | {fail_stats.get('filter_fail', 0):,} | {(fail_stats.get('filter_fail', 0)/vcf2_variants*100):.2f}% |
| DP < 30 | {fail_stats.get('dp_fail', 0):,} | {(fail_stats.get('dp_fail', 0)/vcf2_variants*100):.2f}% |
| GQ < 60 | {fail_stats.get('gq_fail', 0):,} | {(fail_stats.get('gq_fail', 0)/vcf2_variants*100):.2f}% |
| 信息缺失 | {fail_stats.get('missing_info', 0):,} | {(fail_stats.get('missing_info', 0)/vcf2_variants*100):.2f}% |

## 筛选说明
### 质量标准定义
- **QUAL**: Phred质量分数，数值越高表示变异检测越可信
- **FILTER**: "."或"PASS"表示通过所有质控过滤
- **DP**: 支持该变异的测序reads总数量
- **GQ**: 基因型判定的Phred质量分数

### 筛选策略
采用四重质量标准确保：
1. 高可信度的变异检测
2. 完整的质量控制通过
3. 足够的测序覆盖度
4. 准确的基因型判定

## 技术细节
- **处理速度**: {vcf2_variants/filter_time:.1f} 变异/秒
- **筛选算法**: 逐条解析VCF记录，应用四重质量标准
- **保留策略**: 严格筛选，只保留高质量变异

## 下一步
下一步将进行ClinVar数据库注释和过滤。

---
*生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}*
*处理脚本: scripts/step3_sequencing_quality_filter.py*
"""

        with open('results/STEP3_SEQUENCING_QUALITY_FILTER_REPORT.md', 'w') as f:
            f.write(report_content)

        print(f"\n✅ 详细报告已保存到: results/STEP3_SEQUENCING_QUALITY_FILTER_REPORT.md")
        print("🚀 准备进行第四步ClinVar注释...")

    except Exception as e:
        print(f"❌ 处理过程中发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()