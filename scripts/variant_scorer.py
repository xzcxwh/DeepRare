#!/usr/bin/env python3
"""
变异打分系统 v2.0 - 根据多种注释对变异进行综合打分
使用简化的打分逻辑，选取分数最高的50个变异进入决赛圈

打分规则：
- AlphaMissense_score ≥ 0.564：+2
- REVEL_score ≥ 0.5：+2
- CADD_phred ≥ 20：+2
- MetaRNN_score ≥ 0.5：+1
- PrimateAI_score ≥ 0.7：+1
- ClinPred_score ≥ 0.5：+1
- ESM1b_score < -3：+1
- PhyloP ≥ 2：+1

注意：AlphaMissense、CADD、REVEL 可能同时被 VEP 和 dbNSFP 注释，
      取两者中的最大值进行打分。
"""

import re
import sys
import argparse
from collections import defaultdict

def parse_info(info_str):
    """解析 VCF INFO 字段"""
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info

def get_max_score(value_str):
    """从逗号分隔的多值字段中获取最大有效数值"""
    if not value_str or value_str == '.':
        return None
    max_val = None
    for v in value_str.split(','):
        v = v.strip()
        if v and v != '.':
            try:
                val = float(v)
                if max_val is None or val > max_val:
                    max_val = val
            except ValueError:
                continue
    return max_val

def get_min_score(value_str):
    """从逗号分隔的多值字段中获取最小有效数值（用于 ESM1b）"""
    if not value_str or value_str == '.':
        return None
    min_val = None
    for v in value_str.split(','):
        v = v.strip()
        if v and v != '.':
            try:
                val = float(v)
                if min_val is None or val < min_val:
                    min_val = val
            except ValueError:
                continue
    return min_val

def get_merged_score(info, *keys):
    """
    从多个可能的字段中获取最大值
    用于处理 VEP 和 dbNSFP 同时注释的情况
    """
    max_val = None
    for key in keys:
        val = get_max_score(info.get(key, ''))
        if val is not None:
            if max_val is None or val > max_val:
                max_val = val
    return max_val

def score_variant(info):
    """
    综合打分逻辑 v2.0
    
    打分规则（最高 11 分）：
    - AlphaMissense_score ≥ 0.564：+2
    - REVEL_score ≥ 0.5：+2
    - CADD_phred ≥ 20：+2
    - MetaRNN_score ≥ 0.5：+1
    - PrimateAI_score ≥ 0.7：+1
    - ClinPred_score ≥ 0.5：+1
    - ESM1b_score < -3：+1
    - PhyloP ≥ 2：+1
    """
    score = 0
    score_breakdown = {}
    
    # AlphaMissense (VEP 和 dbNSFP 可能都有，取最大值)
    am = get_merged_score(info, 'dbNSFP_AlphaMissense', 'VEP_AlphaMissense', 'AlphaMissense_score')
    if am is not None and am >= 0.564:
        score += 2
        score_breakdown['AlphaMissense'] = f"{am:.3f} (≥0.564) +2"
    else:
        score_breakdown['AlphaMissense'] = f"{am:.3f} (<0.564) +0" if am else "N/A"
    
    # REVEL (VEP 和 dbNSFP 可能都有，取最大值)
    revel = get_merged_score(info, 'dbNSFP_REVEL', 'VEP_REVEL', 'REVEL_score')
    if revel is not None and revel >= 0.5:
        score += 2
        score_breakdown['REVEL'] = f"{revel:.3f} (≥0.5) +2"
    else:
        score_breakdown['REVEL'] = f"{revel:.3f} (<0.5) +0" if revel else "N/A"
    
    # CADD (VEP 和 dbNSFP 可能都有，取最大值)
    cadd = get_merged_score(info, 'dbNSFP_CADD', 'VEP_CADD', 'CADD_phred')
    if cadd is not None and cadd >= 20:
        score += 2
        score_breakdown['CADD'] = f"{cadd:.1f} (≥20) +2"
    else:
        score_breakdown['CADD'] = f"{cadd:.1f} (<20) +0" if cadd else "N/A"
    
    # MetaRNN
    metarnn = get_max_score(info.get('dbNSFP_MetaRNN', '') or info.get('MetaRNN_score', ''))
    if metarnn is not None and metarnn >= 0.5:
        score += 1
        score_breakdown['MetaRNN'] = f"{metarnn:.3f} (≥0.5) +1"
    else:
        score_breakdown['MetaRNN'] = f"{metarnn:.3f} (<0.5) +0" if metarnn else "N/A"
    
    # PrimateAI
    primateai = get_max_score(info.get('dbNSFP_PrimateAI', '') or info.get('PrimateAI_score', ''))
    if primateai is not None and primateai >= 0.7:
        score += 1
        score_breakdown['PrimateAI'] = f"{primateai:.3f} (≥0.7) +1"
    else:
        score_breakdown['PrimateAI'] = f"{primateai:.3f} (<0.7) +0" if primateai else "N/A"
    
    # ClinPred
    clinpred = get_max_score(info.get('dbNSFP_ClinPred', '') or info.get('ClinPred_score', ''))
    if clinpred is not None and clinpred >= 0.5:
        score += 1
        score_breakdown['ClinPred'] = f"{clinpred:.3f} (≥0.5) +1"
    else:
        score_breakdown['ClinPred'] = f"{clinpred:.3f} (<0.5) +0" if clinpred else "N/A"
    
    # ESM1b (更负的值更有害)
    esm1b = get_min_score(info.get('dbNSFP_ESM1b', '') or info.get('ESM1b_score', ''))
    if esm1b is not None and esm1b < -3:
        score += 1
        score_breakdown['ESM1b'] = f"{esm1b:.3f} (<-3) +1"
    else:
        score_breakdown['ESM1b'] = f"{esm1b:.3f} (≥-3) +0" if esm1b else "N/A"
    
    # PhyloP
    phylop = get_max_score(info.get('VEP_PhyloP', '') or info.get('phyloP', ''))
    if phylop is not None and phylop >= 2:
        score += 1
        score_breakdown['PhyloP'] = f"{phylop:.2f} (≥2) +1"
    else:
        score_breakdown['PhyloP'] = f"{phylop:.2f} (<2) +0" if phylop else "N/A"
    
    return score, score_breakdown

def main():
    parser = argparse.ArgumentParser(description='变异打分系统 v2.0')
    parser.add_argument('-i', '--input', required=True, help='输入 VCF 文件')
    parser.add_argument('-o', '--output', help='输出报告文件 (默认: stdout)')
    parser.add_argument('-n', '--top', type=int, default=50, help='输出 TOP N 变异 (默认: 50)')
    parser.add_argument('--vcf-output', help='输出 TOP N 变异的 VCF 文件')
    parser.add_argument('--verbose', action='store_true', help='显示详细分数分解')
    args = parser.parse_args()
    
    variants = []
    header_lines = []
    
    with open(args.input, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.rstrip())
                continue
            
            fields = line.rstrip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]
            info = parse_info(info_str)
            
            score, breakdown = score_variant(info)
            
            # 提取基因名 - 尝试多个可能的字段
            gene = info.get('VEP_Gene', '')
            if not gene:
                # 尝试从 ANN 字段提取
                ann = info.get('ANN', '')
                if ann:
                    ann_parts = ann.split('|')
                    if len(ann_parts) > 3:
                        gene = ann_parts[3]
            if not gene:
                gene = 'Unknown'
            
            # 提取功能影响
            consequence = info.get('VEP_Consequence', '')
            if not consequence:
                ann = info.get('ANN', '')
                if ann:
                    ann_parts = ann.split('|')
                    if len(ann_parts) > 1:
                        consequence = ann_parts[1]
            if not consequence:
                consequence = 'Unknown'
            
            variants.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'gene': gene,
                'consequence': consequence,
                'score': score,
                'breakdown': breakdown,
                'line': line.rstrip()
            })
    
    # 按分数排序（降序）
    variants.sort(key=lambda x: (-x['score'], x['chrom'], int(x['pos'])))
    
    # 确定 TOP N 边界，处理同分情况
    if len(variants) > 0:
        top_n = args.top
        if len(variants) > top_n:
            # 找到第 N 个变异的分数
            threshold_score = variants[top_n - 1]['score']
            # 包含所有同分变异
            actual_top_n = top_n
            while actual_top_n < len(variants) and variants[actual_top_n]['score'] == threshold_score:
                actual_top_n += 1
        else:
            actual_top_n = len(variants)
    else:
        actual_top_n = 0
    
    # 输出结果
    out = open(args.output, 'w') if args.output else sys.stdout
    
    print("=" * 100, file=out)
    print(f"变异打分结果 v2.0 - TOP {args.top}", file=out)
    print("=" * 100, file=out)
    print(f"\n总变异数: {len(variants)}", file=out)
    if len(variants) > 0:
        print(f"分数范围: {variants[-1]['score']} - {variants[0]['score']} (最高 11 分)\n", file=out)
    
    print("\n打分规则:", file=out)
    print("  - AlphaMissense ≥ 0.564: +2", file=out)
    print("  - REVEL ≥ 0.5: +2", file=out)
    print("  - CADD ≥ 20: +2", file=out)
    print("  - MetaRNN ≥ 0.5: +1", file=out)
    print("  - PrimateAI ≥ 0.7: +1", file=out)
    print("  - ClinPred ≥ 0.5: +1", file=out)
    print("  - ESM1b < -3: +1", file=out)
    print("  - PhyloP ≥ 2: +1", file=out)
    print(f"\n选取 TOP {args.top} 变异 (含同分): 实际选取 {actual_top_n} 个\n", file=out)
    
    print("-" * 100, file=out)
    print(f"{'排名':<4} {'位置':<20} {'基因':<15} {'变异':<12} {'功能影响':<30} {'分数':<6}", file=out)
    print("-" * 100, file=out)
    
    for i, v in enumerate(variants[:actual_top_n], 1):
        var_str = f"{v['ref']}>{v['alt']}"
        if len(var_str) > 10:
            var_str = var_str[:10] + "..."
        
        conseq = v['consequence']
        if len(conseq) > 28:
            conseq = conseq[:28] + ".."
        
        gene = v['gene']
        if len(gene) > 13:
            gene = gene[:13] + ".."
        
        print(f"{i:<4} {v['chrom']}:{v['pos']:<12} {gene:<15} {var_str:<12} {conseq:<30} {v['score']:<6}", file=out)
        
        if args.verbose:
            b = v['breakdown']
            print(f"     AM:{b['AlphaMissense']}  REVEL:{b['REVEL']}  CADD:{b['CADD']}", file=out)
            print(f"     MetaRNN:{b['MetaRNN']}  PrimateAI:{b['PrimateAI']}  ClinPred:{b['ClinPred']}", file=out)
            print(f"     ESM1b:{b['ESM1b']}  PhyloP:{b['PhyloP']}", file=out)
            print("", file=out)
    
    print("-" * 100, file=out)
    
    # 分数分布统计
    print("\n分数分布统计:", file=out)
    score_counts = defaultdict(int)
    for v in variants:
        score_counts[v['score']] += 1
    
    for score in sorted(score_counts.keys(), reverse=True):
        count = score_counts[score]
        bar = "█" * min(count, 50)
        print(f"  分数 {score:>2}: {count:>4} 个 {bar}", file=out)
    
    if args.output:
        out.close()
        print(f"结果已保存到: {args.output}")
    
    # 输出 TOP N 的 VCF 文件（包含同分变异）
    if args.vcf_output:
        with open(args.vcf_output, 'w') as vcf_out:
            # 写入原始头部
            for h in header_lines:
                vcf_out.write(h + '\n')
            
            # 写入 TOP N 变异（包含同分）
            for v in variants[:actual_top_n]:
                # 在 INFO 字段中添加 HERITA_SCORE
                fields = v['line'].split('\t')
                fields[7] = f"HERITA_SCORE={v['score']};" + fields[7]
                vcf_out.write('\t'.join(fields) + '\n')
        
        print(f"TOP {actual_top_n} VCF 已保存到: {args.vcf_output}")

if __name__ == '__main__':
    main()
