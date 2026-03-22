#!/usr/bin/env python3
"""
步骤6: 使用dbNSFP数据库对变异进行多个有害性评分注释，并根据条件筛选变异
"""

import subprocess
import os
import sys
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 文件路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
INPUT_VCF = os.path.join(PROJECT_ROOT, "results", "vcf5.vcf")
OUTPUT_VCF = os.path.join(PROJECT_ROOT, "results", "vcf6.vcf")
ANNOTATED_VCF = os.path.join(PROJECT_ROOT, "results", "vcf5_annotated.vcf")
VCFANNO_CONFIG = os.path.join(PROJECT_ROOT, "config", "vcfanno_multiple_scores_config.toml")

def check_files_exist():
    """检查必要的文件是否存在"""
    files_to_check = [INPUT_VCF, VCFANNO_CONFIG]
    for file_path in files_to_check:
        if not os.path.exists(file_path):
            logger.error(f"文件不存在: {file_path}")
            return False
    return True

def run_vcfanno_annotation():
    """使用vcfanno进行注释"""
    logger.info("开始使用vcfanno注释变异...")

    cmd = [
        "vcfanno",
        "-p", "8",  # 使用8个线程
        VCFANNO_CONFIG,
        INPUT_VCF
    ]

    try:
        with open(ANNOTATED_VCF, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)

        logger.info(f"vcfanno注释完成，结果保存到: {ANNOTATED_VCF}")
        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"vcfanno注释失败: {e}")
        logger.error(f"错误输出: {e.stderr}")
        return False

def filter_variants():
    """根据多个有害性评分筛选变异"""
    logger.info("开始根据有害性评分筛选变异...")

    try:
        # 首先提取所有需要的信息字段
        cmd = [
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/MetaRNN_score\t%INFO/REVEL_score\t%INFO/PrimateAI_score\t%INFO/ClinPred_score\t%INFO/ESM1b_score\t%INFO/AlphaMissense_score\t%INFO/CADD_phred\n",
            ANNOTATED_VCF
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')

        # 筛选条件：
        # MetaRNN_score ≥ 0.5
        # REVEL_score ≥ 0.5
        # PrimateAI_score ≥ 0.7
        # ClinPred_score ≥ 0.5
        # ESM1b_score < -3
        # AlphaMissense_score ≥ 0.564
        # CADD_phred ≥ 20
        # 至少满足上述2个及以上的条件

        passing_variants = []

        for line in lines:
            if not line:
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            chrom, pos, ref, alt = fields[:4]
            meta_rnn, revel, primate_ai, clinpred, esm1b, alphamissense, cadd = fields[4:11]

            # 计算满足的条件数量
            conditions_met = 0

            # 处理可能包含多个值的字段（以逗号分隔）
            # MetaRNN_score ≥ 0.5
            if meta_rnn != ".":
                meta_rnn_values = [float(x) for x in meta_rnn.split(',') if x != '.']
                if meta_rnn_values and max(meta_rnn_values) >= 0.5:
                    conditions_met += 1

            # REVEL_score ≥ 0.5
            if revel != ".":
                revel_values = [float(x) for x in revel.split(',') if x != '.']
                if revel_values and max(revel_values) >= 0.5:
                    conditions_met += 1

            # PrimateAI_score ≥ 0.7
            if primate_ai != ".":
                primate_ai_values = [float(x) for x in primate_ai.split(',') if x != '.']
                if primate_ai_values and max(primate_ai_values) >= 0.7:
                    conditions_met += 1

            # ClinPred_score ≥ 0.5
            if clinpred != ".":
                clinpred_values = [float(x) for x in clinpred.split(',') if x != '.']
                if clinpred_values and max(clinpred_values) >= 0.5:
                    conditions_met += 1

            # ESM1b_score < -3
            if esm1b != ".":
                esm1b_values = [float(x) for x in esm1b.split(',') if x != '.']
                if esm1b_values and min(esm1b_values) < -3:
                    conditions_met += 1

            # AlphaMissense_score ≥ 0.564
            if alphamissense != ".":
                alphamissense_values = [float(x) for x in alphamissense.split(',') if x != '.']
                if alphamissense_values and max(alphamissense_values) >= 0.564:
                    conditions_met += 1

            # CADD_phred ≥ 20
            if cadd != ".":
                cadd_values = [float(x) for x in cadd.split(',') if x != '.']
                if cadd_values and max(cadd_values) >= 20:
                    conditions_met += 1

            # 至少满足2个条件
            if conditions_met >= 2:
                passing_variants.append(f"{chrom}_{pos}_{ref}_{alt}")

        logger.info(f"找到 {len(passing_variants)} 个满足至少2个条件的变异")

        # 使用bcftools提取满足条件的变异
        if passing_variants:
            # 创建位置列表用于筛选
            positions_file = os.path.join(PROJECT_ROOT, "results", "passing_positions.txt")
            with open(positions_file, 'w') as f:
                for variant in passing_variants:
                    parts = variant.split('_')
                    if len(parts) >= 4:
                        chrom, pos, ref, alt = parts[:4]
                        f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

            # 使用bcftools view -T文件来提取变异
            cmd = [
                "bcftools", "view",
                "-T", positions_file,
                "-Ov",
                "-o", OUTPUT_VCF,
                ANNOTATED_VCF
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # 清理临时文件
            os.remove(positions_file)

            logger.info(f"变异筛选完成，结果保存到: {OUTPUT_VCF}")
        else:
            # 如果没有变异满足条件，创建一个空的VCF文件
            logger.warning("没有变异满足筛选条件，创建空的VCF文件")
            cmd = ["bcftools", "view", "-i", "0==1", "-Ov", "-o", OUTPUT_VCF, ANNOTATED_VCF]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"创建空VCF文件: {OUTPUT_VCF}")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"变异筛选失败: {e}")
        logger.error(f"错误输出: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"变异筛选过程中发生错误: {e}")
        return False

def count_variants():
    """统计变异数量"""
    try:
        # 统计原始vcf5中的变异数量
        cmd_original = ["bcftools", "view", "-H", INPUT_VCF]
        result_original = subprocess.run(cmd_original, capture_output=True, text=True, check=True)
        original_count = len(result_original.stdout.strip().split('\n')) if result_original.stdout.strip() else 0

        # 统计注释后vcf5_annotated中的变异数量
        cmd_annotated = ["bcftools", "view", "-H", ANNOTATED_VCF]
        result_annotated = subprocess.run(cmd_annotated, capture_output=True, text=True, check=True)
        annotated_count = len(result_annotated.stdout.strip().split('\n')) if result_annotated.stdout.strip() else 0

        # 统计筛选后vcf6中的变异数量
        cmd_filtered = ["bcftools", "view", "-H", OUTPUT_VCF]
        result_filtered = subprocess.run(cmd_filtered, capture_output=True, text=True, check=True)
        filtered_count = len(result_filtered.stdout.strip().split('\n')) if result_filtered.stdout.strip() else 0

        logger.info(f"原始变异数量 (vcf5): {original_count}")
        logger.info(f"注释后变异数量 (vcf5_annotated): {annotated_count}")
        logger.info(f"筛选后变异数量 (vcf6): {filtered_count}")
        logger.info(f"筛选比例: {filtered_count/original_count*100:.2f}%" if original_count > 0 else "N/A")

        return original_count, annotated_count, filtered_count

    except subprocess.CalledProcessError as e:
        logger.error(f"统计变异数量失败: {e}")
        return None, None, None

def main():
    """主函数"""
    logger.info("开始执行步骤6: dbNSFP注释和变异筛选")

    # 检查文件是否存在
    if not check_files_exist():
        sys.exit(1)

    # 创建输出目录
    os.makedirs(os.path.dirname(OUTPUT_VCF), exist_ok=True)

    # 步骤1: 使用vcfanno进行注释
    if not run_vcfanno_annotation():
        sys.exit(1)

    # 步骤2: 根据条件筛选变异
    if not filter_variants():
        sys.exit(1)

    # 步骤3: 统计变异数量
    original_count, annotated_count, filtered_count = count_variants()

    if filtered_count is not None:
        logger.info("步骤6执行完成!")
        logger.info(f"成功从 {original_count} 个变异中筛选出 {filtered_count} 个满足条件的变异")
    else:
        logger.error("步骤6执行失败")
        sys.exit(1)

if __name__ == "__main__":
    main()