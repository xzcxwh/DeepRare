#!/usr/bin/env python3
"""
步骤7: 对vcf6中的变异进行打分，选出得分最高的50个变异进入决赛圈
"""

import subprocess
import os
import sys
import logging
import csv
from collections import defaultdict

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 文件路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
INPUT_VCF = os.path.join(PROJECT_ROOT, "results", "vcf6.vcf")
OUTPUT_VCF = os.path.join(PROJECT_ROOT, "results", "final_top.vcf")
REPORT_FILE = os.path.join(PROJECT_ROOT, "docs", "final_top_report.md")

def calculate_variant_score(meta_rnn, revel, primate_ai, clinpred, esm1b, alphamissense, cadd):
    """
    计算变异得分

    评分规则：
    - AlphaMissense_score ≥ 0.564: +2分
    - ClinPred_score ≥ 0.5: +2分
    - MetaRNN_score ≥ 0.5: +2分
    - REVEL_score ≥ 0.5: +2分
    - CADD_phred ≥ 20: +1分
    - PrimateAI_score ≥ 0.7: +1分
    - ESM1b_score < -3: +1分
    """
    score = 0
    score_details = []

    # AlphaMissense_score ≥ 0.564: +2分
    if alphamissense != ".":
        alphamissense_values = [float(x) for x in alphamissense.split(',') if x != '.']
        if alphamissense_values and max(alphamissense_values) >= 0.564:
            score += 2
            score_details.append(f"AlphaMissense≥0.564(+2分, max={max(alphamissense_values):.3f})")

    # ClinPred_score ≥ 0.5: +2分
    if clinpred != ".":
        clinpred_values = [float(x) for x in clinpred.split(',') if x != '.']
        if clinpred_values and max(clinpred_values) >= 0.5:
            score += 2
            score_details.append(f"ClinPred≥0.5(+2分, max={max(clinpred_values):.3f})")

    # MetaRNN_score ≥ 0.5: +2分
    if meta_rnn != ".":
        meta_rnn_values = [float(x) for x in meta_rnn.split(',') if x != '.']
        if meta_rnn_values and max(meta_rnn_values) >= 0.5:
            score += 2
            score_details.append(f"MetaRNN≥0.5(+2分, max={max(meta_rnn_values):.3f})")

    # REVEL_score ≥ 0.5: +2分
    if revel != ".":
        revel_values = [float(x) for x in revel.split(',') if x != '.']
        if revel_values and max(revel_values) >= 0.5:
            score += 2
            score_details.append(f"REVEL≥0.5(+2分, max={max(revel_values):.3f})")

    # CADD_phred ≥ 20: +1分
    if cadd != ".":
        cadd_values = [float(x) for x in cadd.split(',') if x != '.']
        if cadd_values and max(cadd_values) >= 20:
            score += 1
            score_details.append(f"CADD≥20(+1分, max={max(cadd_values):.1f})")

    # PrimateAI_score ≥ 0.7: +1分
    if primate_ai != ".":
        primate_ai_values = [float(x) for x in primate_ai.split(',') if x != '.']
        if primate_ai_values and max(primate_ai_values) >= 0.7:
            score += 1
            score_details.append(f"PrimateAI≥0.7(+1分, max={max(primate_ai_values):.3f})")

    # ESM1b_score < -3: +1分
    if esm1b != ".":
        esm1b_values = [float(x) for x in esm1b.split(',') if x != '.']
        if esm1b_values and min(esm1b_values) < -3:
            score += 1
            score_details.append(f"ESM1b<-3(+1分, min={min(esm1b_values):.3f})")

    return score, score_details

def score_variants():
    """对所有变异进行打分"""
    logger.info("开始对变异进行打分...")

    cmd = [
        "bcftools", "query",
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/MetaRNN_score\t%INFO/REVEL_score\t%INFO/PrimateAI_score\t%INFO/ClinPred_score\t%INFO/ESM1b_score\t%INFO/AlphaMissense_score\t%INFO/CADD_phred\n",
        INPUT_VCF
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')

        scored_variants = []

        for line in lines:
            if not line:
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue

            chrom, pos, ref, alt = fields[:4]
            meta_rnn, revel, primate_ai, clinpred, esm1b, alphamissense, cadd = fields[4:11]

            score, score_details = calculate_variant_score(meta_rnn, revel, primate_ai, clinpred, esm1b, alphamissense, cadd)

            variant_data = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'score': score,
                'score_details': score_details,
                'meta_rnn': meta_rnn,
                'revel': revel,
                'primate_ai': primate_ai,
                'clindpred': clinpred,
                'esm1b': esm1b,
                'alphamissense': alphamissense,
                'cadd': cadd
            }

            scored_variants.append(variant_data)

        logger.info(f"成功为 {len(scored_variants)} 个变异打分")
        return scored_variants

    except subprocess.CalledProcessError as e:
        logger.error(f"获取变异信息失败: {e}")
        return None

def select_top_variants(scored_variants, min_n=50):
    """选择得分最高的变异，包含与第N名得分相同的所有变异"""
    logger.info(f"选择得分最高的变异，确保至少包含 {min_n} 个变异...")

    # 按得分降序排列，得分相同则按位置排序
    sorted_variants = sorted(scored_variants, key=lambda x: (-x['score'], int(x['pos'])))

    if len(sorted_variants) < min_n:
        logger.warning(f"警告: 只有 {len(sorted_variants)} 个变异，少于要求的 {min_n} 个")
        return sorted_variants

    # 获取第N名的得分
    nth_score = sorted_variants[min_n - 1]['score']

    # 选择得分 >= 第N名得分的所有变异
    top_variants = [v for v in sorted_variants if v['score'] >= nth_score]

    logger.info(f"成功选出 {len(top_variants)} 个变异")
    logger.info(f"第{min_n}名得分为: {nth_score}分")
    logger.info(f"得分范围: {top_variants[0]['score']} - {top_variants[-1]['score']}")

    return top_variants

def create_final_vcf(top_variants):
    """创建final_top.vcf文件"""
    logger.info("创建final_top.vcf文件...")

    try:
        # 创建位置列表用于筛选
        positions_file = os.path.join(PROJECT_ROOT, "results", "top_positions.txt")
        with open(positions_file, 'w') as f:
            for variant in top_variants:
                f.write(f"{variant['chrom']}\t{variant['pos']}\t{variant['ref']}\t{variant['alt']}\n")

        # 使用bcftools view -T文件来提取变异
        cmd = [
            "bcftools", "view",
            "-T", positions_file,
            "-Ov",
            "-o", OUTPUT_VCF,
            INPUT_VCF
        ]

        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        # 清理临时文件
        os.remove(positions_file)

        logger.info(f"成功创建 {OUTPUT_VCF}")

        # 验证文件
        cmd_count = ["bcftools", "view", "-H", OUTPUT_VCF]
        result_count = subprocess.run(cmd_count, capture_output=True, text=True, check=True)
        variant_count = len(result_count.stdout.strip().split('\n')) if result_count.stdout.strip() else 0

        logger.info(f"final_top.vcf包含 {variant_count} 个变异")

        if variant_count != len(top_variants):
            logger.warning(f"警告: 期望 {len(top_variants)} 个变异，实际得到 {variant_count} 个")

        return True

    except subprocess.CalledProcessError as e:
        logger.error(f"创建final_top50.vcf失败: {e}")
        logger.error(f"错误输出: {e.stderr}")
        return False

def generate_report(top_variants):
    """生成详细的报告"""
    logger.info("生成详细报告...")

    # 确保docs目录存在
    os.makedirs(os.path.dirname(REPORT_FILE), exist_ok=True)

    try:
        with open(REPORT_FILE, 'w', encoding='utf-8') as f:
            f.write("# Final Top 变异报告\n\n")
            f.write(f"本报告详细列出了从149个候选变异中筛选出的得分最高的变异。\n")
            f.write(f"筛选规则：选择得分最高的前50个变异，同时包含与第50名得分相同的所有变异。\n\n")

            # 评分标准说明
            f.write("## 评分标准\n\n")
            f.write("每个变异根据以下标准进行打分：\n\n")
            f.write("| 评分条件 | 分值 | 说明 |\n")
            f.write("|---------|------|------|\n")
            f.write("| AlphaMissense_score ≥ 0.564 | +2分 | AlphaMissense深度学习预测分数 |\n")
            f.write("| ClinPred_score ≥ 0.5 | +2分 | ClinPred临床预测分数 |\n")
            f.write("| MetaRNN_score ≥ 0.5 | +2分 | MetaRNN集成学习分数 |\n")
            f.write("| REVEL_score ≥ 0.5 | +2分 | REVEL集成预测分数 |\n")
            f.write("| CADD_phred ≥ 20 | +1分 | CADD有害性预测分数 |\n")
            f.write("| PrimateAI_score ≥ 0.7 | +1分 | PrimateAI灵长类保守性分数 |\n")
            f.write("| ESM1b_score < -3 | +1分 | ESM1b蛋白质结构影响分数 |\n")
            f.write("\n")
            f.write("最高可能得分：11分\n\n")

            # 统计信息
            f.write("## 统计信息\n\n")
            f.write(f"- 候选变异数量：149个\n")
            f.write(f"- 选中变异数量：{len(top_variants)}个\n")
            f.write(f"- 得分范围：{top_variants[0]['score']} - {top_variants[-1]['score']}\n")

            # 计算第50名的得分（如果总数>=50）
            if len(top_variants) >= 50:
                f.write(f"- 第50名得分：{top_variants[49]['score']}分\n")
            f.write("\n")

            # 得分分布
            score_distribution = defaultdict(int)
            for variant in top_variants:
                score_distribution[variant['score']] += 1

            f.write("### 得分分布\n\n")
            f.write("| 得分 | 变异数量 |\n")
            f.write("|------|----------|\n")
            for score in sorted(score_distribution.keys(), reverse=True):
                f.write(f"| {score} | {score_distribution[score]} |\n")
            f.write("\n")

            # 详细变异列表
            f.write(f"## 前{len(top_variants)}个变异详细信息\n\n")
            f.write("按得分降序排列：\n\n")

            for i, variant in enumerate(top_variants, 1):
                f.write(f"### {i}. {variant['chrom']}:{variant['pos']} {variant['ref']}>{variant['alt']} (得分: {variant['score']})\n\n")

                f.write("**满足的评分条件：**\n")
                if variant['score_details']:
                    for detail in variant['score_details']:
                        f.write(f"- {detail}\n")
                else:
                    f.write("无\n")
                f.write("\n")

                f.write("**详细分数：**\n")
                f.write(f"- MetaRNN_score: {variant['meta_rnn']}\n")
                f.write(f"- REVEL_score: {variant['revel']}\n")
                f.write(f"- PrimateAI_score: {variant['primate_ai']}\n")
                f.write(f"- ClinPred_score: {variant['clindpred']}\n")
                f.write(f"- ESM1b_score: {variant['esm1b']}\n")
                f.write(f"- AlphaMissense_score: {variant['alphamissense']}\n")
                f.write(f"- CADD_phred: {variant['cadd']}\n\n")

                f.write("---\n\n")

            # 总结
            f.write("## 总结\n\n")
            f.write(f"本分析从149个通过初步筛选的变异中，基于多个有害性预测算法的综合评分，\n")
            f.write(f"选出了得分最高的{len(top_variants)}个变异进入决赛圈。筛选规则确保了至少包含50个变异，\n")
            f.write(f"同时包含与第50名得分相同的所有变异。这些变异在多个预测工具中都表现出\n")
            f.write(f"较高的有害性预测，值得进一步的实验验证和临床研究。\n\n")

            f.write(f"分析完成时间：{subprocess.run(['date'], capture_output=True, text=True).stdout.strip()}\n")

        logger.info(f"成功生成报告: {REPORT_FILE}")
        return True

    except Exception as e:
        logger.error(f"生成报告失败: {e}")
        return False

def main():
    """主函数"""
    logger.info("开始执行步骤7: 最终Top变异选择")

    # 检查输入文件
    if not os.path.exists(INPUT_VCF):
        logger.error(f"输入文件不存在: {INPUT_VCF}")
        sys.exit(1)

    # 步骤1: 对所有变异进行打分
    scored_variants = score_variants()
    if scored_variants is None:
        sys.exit(1)

    # 步骤2: 选择得分最高的变异（至少50个）
    top_variants = select_top_variants(scored_variants, 50)

    # 步骤3: 创建final_top.vcf文件
    if not create_final_vcf(top_variants):
        sys.exit(1)

    # 步骤4: 生成详细报告
    if not generate_report(top_variants):
        sys.exit(1)

    logger.info("步骤7执行完成!")
    logger.info(f"成功选出{len(top_variants)}个最高得分变异保存到: {OUTPUT_VCF}")
    logger.info(f"详细报告保存到: {REPORT_FILE}")

if __name__ == "__main__":
    main()