#!/usr/bin/env python3
"""
HERITA ACMG Classifier
根据ACMG/AMP指南对变异进行致病性分类
参考InterVar的评级逻辑和知识库
"""

import sys
import os
import re
import argparse
from collections import defaultdict

# 知识库路径
INTERVAR_DB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                           "InterVar-master", "intervardb")

# 全局知识库字典
lof_genes = set()           # PVS1: Loss-of-function genes
pp2_genes = set()           # PP2: Missense pathogenic genes
bp1_genes = set()           # BP1: Truncating genes
ps1_aa_changes = {}         # PS1: Same AA change as pathogenic
pm1_domains = {}            # PM1: Functional domains without benign
ps4_variants = set()        # PS4: GWAS significant variants
mim_recessive = set()       # Recessive disease genes
mim_dominant = set()        # Dominant disease genes
bs2_variants = {}           # BS2: Observed in healthy individuals

def load_knowledge_bases():
    """加载InterVar知识库"""
    global lof_genes, pp2_genes, bp1_genes, ps1_aa_changes, pm1_domains
    global ps4_variants, mim_recessive, mim_dominant, bs2_variants
    
    # PVS1: LOF genes
    lof_file = os.path.join(INTERVAR_DB, "PVS1.LOF.genes.hg38")
    if os.path.exists(lof_file):
        with open(lof_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    lof_genes.add(gene)
    
    # PP2: Missense pathogenic genes
    pp2_file = os.path.join(INTERVAR_DB, "PP2.genes.hg38")
    if os.path.exists(pp2_file):
        with open(pp2_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    pp2_genes.add(gene)
    
    # BP1: Truncating genes
    bp1_file = os.path.join(INTERVAR_DB, "BP1.genes.hg38")
    if os.path.exists(bp1_file):
        with open(bp1_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    bp1_genes.add(gene)
    
    # PS1: Same AA change as pathogenic
    ps1_file = os.path.join(INTERVAR_DB, "PS1.AA.change.patho.hg38")
    if os.path.exists(ps1_file):
        with open(ps1_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    # chr_pos_ref_aa
                    key = f"{parts[0]}_{parts[1]}_{parts[3]}_{parts[5]}"
                    ps1_aa_changes[key] = parts[5]  # pathogenic AA
    
    # PM1: Functional domains
    pm1_file = os.path.join(INTERVAR_DB, "PM1_domains_with_benigns.hg38")
    if os.path.exists(pm1_file):
        with open(pm1_file, 'r') as f:
            next(f, None)  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    key = f"{parts[0]}_{parts[1]}: {parts[2]}"
                    pm1_domains[key] = True
    
    # PS4: GWAS significant variants
    ps4_file = os.path.join(INTERVAR_DB, "PS4.variants.hg38")
    if os.path.exists(ps4_file):
        with open(ps4_file, 'r') as f:
            next(f, None)  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    key = f"{parts[0]}_{parts[1]}_{parts[3]}_{parts[4]}"
                    ps4_variants.add(key)
    
    # MIM recessive
    mim_rec_file = os.path.join(INTERVAR_DB, "mim_recessive.txt")
    if os.path.exists(mim_rec_file):
        with open(mim_rec_file, 'r') as f:
            for line in f:
                mim = line.strip()
                if mim:
                    mim_recessive.add(mim)
    
    # MIM dominant
    mim_dom_file = os.path.join(INTERVAR_DB, "mim_domin.txt")
    if os.path.exists(mim_dom_file):
        with open(mim_dom_file, 'r') as f:
            for line in f:
                mim = line.strip()
                if mim:
                    mim_dominant.add(mim)
    
    # BS2: Observed in healthy individuals
    bs2_file = os.path.join(INTERVAR_DB, "BS2_hom_het.hg38")
    if os.path.exists(bs2_file):
        with open(bs2_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 6:
                    key = f"{parts[0]}_{parts[1]}_{parts[2]}_{parts[3]}"
                    # parts[4] = hom count, parts[5] = het count
                    bs2_variants[key] = {'hom': int(parts[4]), 'het': int(parts[5])}
    
    print(f"知识库加载完成:", file=sys.stderr)
    print(f"  LOF genes: {len(lof_genes)}", file=sys.stderr)
    print(f"  PP2 genes: {len(pp2_genes)}", file=sys.stderr)
    print(f"  BP1 genes: {len(bp1_genes)}", file=sys.stderr)
    print(f"  PS1 AA changes: {len(ps1_aa_changes)}", file=sys.stderr)
    print(f"  PM1 domains: {len(pm1_domains)}", file=sys.stderr)
    print(f"  PS4 variants: {len(ps4_variants)}", file=sys.stderr)
    print(f"  BS2 variants: {len(bs2_variants)}", file=sys.stderr)


class Variant:
    """变异信息类"""
    def __init__(self, chrom, pos, ref, alt, info, vcf_line):
        self.chrom = chrom.replace('chr', '')
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info
        self.vcf_line = vcf_line
        self.parsed_info = self._parse_info(info)
        
    def _parse_info(self, info_str):
        """解析INFO字段"""
        result = {}
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                result[key] = value
            else:
                result[item] = True
        return result
    
    def get_info(self, key, default=None):
        """获取INFO字段值"""
        return self.parsed_info.get(key, default)
    
    def get_gene(self):
        """从ANN字段提取基因名"""
        ann = self.get_info('ANN', '')
        if ann:
            parts = ann.split('|')
            if len(parts) > 3:
                return parts[3]
        return None
    
    def get_effect(self):
        """从ANN字段提取功能影响"""
        ann = self.get_info('ANN', '')
        if ann:
            parts = ann.split('|')
            if len(parts) > 1:
                return parts[1]
        return None
    
    def get_impact(self):
        """从ANN字段提取影响等级"""
        ann = self.get_info('ANN', '')
        if ann:
            parts = ann.split('|')
            if len(parts) > 2:
                return parts[2]
        return None
    
    def get_aa_change(self):
        """从ANN字段提取氨基酸变化"""
        ann = self.get_info('ANN', '')
        if ann:
            parts = ann.split('|')
            if len(parts) > 10:
                return parts[10]  # p.XXX
        return None
    
    def get_gnomad_af(self):
        """获取gnomAD频率"""
        af = self.get_info('gnomAD_AF_grpmax')
        if af and af != '.':
            try:
                return float(af)
            except:
                pass
        return None
    
    def get_clinvar_sig(self):
        """获取ClinVar意义"""
        return self.get_info('ClinVar_CLNSIG')
    
    def get_spliceai_max(self):
        """获取SpliceAI最大分数"""
        max_score = 0
        for key in ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL']:
            val = self.get_info(key)
            if val and val != '.':
                try:
                    max_score = max(max_score, float(val))
                except:
                    pass
        return max_score
    
    def get_score(self, field):
        """获取数值型分数，同时检查VEP和dbNSFP注释"""
        # 映射字段名
        field_map = {
            'AlphaMissense': ['AlphaMissense', 'AlphaMissense_score', 'VEP_AlphaMissense'],
            'REVEL': ['REVEL', 'REVEL_score', 'VEP_REVEL'],
            'CADD': ['CADD', 'CADD_phred', 'VEP_CADD'],
            'MetaRNN': ['MetaRNN_score'],
            'PrimateAI': ['PrimateAI_score'],
            'ClinPred': ['ClinPred_score'],
            'ESM1b': ['ESM1b_score'],
            'phyloP': ['phyloP', 'VEP_phyloP'],
        }
        
        fields_to_check = field_map.get(field, [field])
        values = []
        
        for f in fields_to_check:
            val = self.get_info(f)
            if val and val != '.' and val != 'NA':
                try:
                    values.append(float(val))
                except:
                    pass
        
        if values:
            return max(values)  # 取最大值
        return None


class ACMGClassifier:
    """ACMG分类器"""
    
    def __init__(self, verbose=False):
        self.verbose = verbose
    
    def check_PVS1(self, variant):
        """
        PVS1: 在LOF是已知致病机制的基因中，出现以下变异类型：
        - 无义突变 (nonsense/stop_gained)
        - 移码突变 (frameshift)
        - 经典剪接位点变异 (splice_donor/acceptor)
        - 起始密码子突变 (start_lost)
        """
        effect = variant.get_effect() or ''
        gene = variant.get_gene()
        
        lof_effects = ['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 
                       'start_lost', 'nonsense']
        
        is_lof_effect = any(eff in effect.lower() for eff in lof_effects)
        is_lof_gene = gene in lof_genes if gene else False
        
        if is_lof_effect and is_lof_gene:
            return 1
        return 0
    
    def check_PS1(self, variant):
        """
        PS1: 与已知致病变异产生相同的氨基酸变化
        """
        effect = variant.get_effect() or ''
        if 'missense' not in effect.lower():
            return 0
        
        # 检查是否存在相同位置的已知致病AA变化
        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        if key in ps1_aa_changes:
            return 1
        return 0
    
    def check_PS3(self, variant):
        """
        PS3: 功能研究支持致病效应
        (通过计算预测器分数来近似)
        """
        score = 0
        alpha = variant.get_score('AlphaMissense')
        revel = variant.get_score('REVEL')
        
        # 如果多个高分预测器一致预测致病
        if alpha and alpha >= 0.9:
            score += 1
        if revel and revel >= 0.9:
            score += 1
        
        return 1 if score >= 2 else 0
    
    def check_PS4(self, variant):
        """
        PS4: 在患者中显著富集 (GWAS)
        """
        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        if key in ps4_variants:
            return 1
        return 0
    
    def check_PM1(self, variant):
        """
        PM1: 位于突变热点或功能域，且无良性变异
        """
        effect = variant.get_effect() or ''
        if 'missense' not in effect.lower():
            return 0
        
        gene = variant.get_gene()
        # 简化：检查基因是否在功能域列表中
        for key in pm1_domains:
            if gene and gene in key:
                return 1
        return 0
    
    def check_PM2(self, variant):
        """
        PM2: 在对照人群中缺失或极低频率
        """
        af = variant.get_gnomad_af()
        if af is None or af == 0:
            return 1
        if af < 0.0001:  # 极低频率
            return 1
        return 0
    
    def check_PM4(self, variant):
        """
        PM4: 非重复区的框内缺失/插入导致蛋白长度变化，或终止密码子丢失
        """
        effect = variant.get_effect() or ''
        if 'inframe_deletion' in effect.lower() or 'inframe_insertion' in effect.lower():
            return 1
        if 'stop_lost' in effect.lower():
            return 1
        return 0
    
    def check_PM5(self, variant):
        """
        PM5: 在已知致病变异位点发生不同的错义变化
        """
        effect = variant.get_effect() or ''
        if 'missense' not in effect.lower():
            return 0
        
        # 检查相同位置是否有其他致病变异
        for key in ps1_aa_changes:
            parts = key.split('_')
            if len(parts) >= 2:
                if parts[0] == variant.chrom and parts[1] == str(variant.pos):
                    # 同位置不同变异
                    if parts[3] != variant.alt:
                        return 1
        return 0
    
    def check_PP2(self, variant):
        """
        PP2: 错义变异位于良性错义变异率低的基因
        """
        effect = variant.get_effect() or ''
        gene = variant.get_gene()
        
        if 'missense' in effect.lower() and gene in pp2_genes:
            return 1
        return 0
    
    def check_PP3(self, variant):
        """
        PP3: 多个计算证据支持致病效应
        """
        evidence_count = 0
        
        # AlphaMissense >= 0.564
        alpha = variant.get_score('AlphaMissense')
        if alpha and alpha >= 0.564:
            evidence_count += 1
        
        # REVEL >= 0.5
        revel = variant.get_score('REVEL')
        if revel and revel >= 0.5:
            evidence_count += 1
        
        # CADD >= 20
        cadd = variant.get_score('CADD')
        if cadd and cadd >= 20:
            evidence_count += 1
        
        # MetaRNN >= 0.5
        metarnn = variant.get_score('MetaRNN')
        if metarnn and metarnn >= 0.5:
            evidence_count += 1
        
        # phyloP >= 2 (保守性)
        phylop = variant.get_score('phyloP')
        if phylop and phylop >= 2:
            evidence_count += 1
        
        # SpliceAI > 0.5
        spliceai = variant.get_spliceai_max()
        if spliceai > 0.5:
            evidence_count += 1
        
        return 1 if evidence_count >= 2 else 0
    
    def check_PP5(self, variant):
        """
        PP5: ClinVar报告为致病
        """
        clinsig = variant.get_clinvar_sig()
        if clinsig:
            clinsig_lower = clinsig.lower()
            if 'pathogenic' in clinsig_lower and 'conflicting' not in clinsig_lower:
                return 1
        return 0
    
    def check_BA1(self, variant):
        """
        BA1: 等位基因频率 > 5%
        """
        af = variant.get_gnomad_af()
        if af and af > 0.05:
            return 1
        return 0
    
    def check_BS1(self, variant):
        """
        BS1: 等位基因频率高于预期 (> 1%)
        """
        af = variant.get_gnomad_af()
        if af and af > 0.01:
            return 1
        return 0
    
    def check_BS2(self, variant):
        """
        BS2: 在健康成人中观察到
        """
        key = f"{variant.chrom}_{variant.pos}_{variant.ref}_{variant.alt}"
        if key in bs2_variants:
            return 1
        return 0
    
    def check_BP1(self, variant):
        """
        BP1: 在主要由截断变异致病的基因中的错义变异
        """
        effect = variant.get_effect() or ''
        gene = variant.get_gene()
        
        if 'missense' in effect.lower() and gene in bp1_genes:
            return 1
        return 0
    
    def check_BP4(self, variant):
        """
        BP4: 多个计算证据支持无影响
        """
        evidence_count = 0
        
        # AlphaMissense < 0.34 (benign)
        alpha = variant.get_score('AlphaMissense')
        if alpha is not None and alpha < 0.34:
            evidence_count += 1
        
        # REVEL < 0.25
        revel = variant.get_score('REVEL')
        if revel is not None and revel < 0.25:
            evidence_count += 1
        
        # CADD < 15
        cadd = variant.get_score('CADD')
        if cadd is not None and cadd < 15:
            evidence_count += 1
        
        # phyloP < 1 (不保守)
        phylop = variant.get_score('phyloP')
        if phylop is not None and phylop < 1:
            evidence_count += 1
        
        return 1 if evidence_count >= 2 else 0
    
    def check_BP6(self, variant):
        """
        BP6: ClinVar报告为良性
        """
        clinsig = variant.get_clinvar_sig()
        if clinsig:
            clinsig_lower = clinsig.lower()
            if 'benign' in clinsig_lower:
                return 1
        return 0
    
    def check_BP7(self, variant):
        """
        BP7: 同义变异且不影响剪接
        """
        effect = variant.get_effect() or ''
        if 'synonymous' in effect.lower():
            spliceai = variant.get_spliceai_max()
            if spliceai < 0.2:
                return 1
        return 0
    
    def classify(self, variant):
        """
        根据ACMG标准对变异进行分类
        返回: (分类结果, 证据字典)
        """
        # 收集所有证据
        evidence = {
            'PVS1': self.check_PVS1(variant),
            'PS1': self.check_PS1(variant),
            'PS3': self.check_PS3(variant),
            'PS4': self.check_PS4(variant),
            'PM1': self.check_PM1(variant),
            'PM2': self.check_PM2(variant),
            'PM4': self.check_PM4(variant),
            'PM5': self.check_PM5(variant),
            'PP2': self.check_PP2(variant),
            'PP3': self.check_PP3(variant),
            'PP5': self.check_PP5(variant),
            'BA1': self.check_BA1(variant),
            'BS1': self.check_BS1(variant),
            'BS2': self.check_BS2(variant),
            'BP1': self.check_BP1(variant),
            'BP4': self.check_BP4(variant),
            'BP6': self.check_BP6(variant),
            'BP7': self.check_BP7(variant),
        }
        
        # 计算证据数量
        PVS1 = evidence['PVS1']
        PS_sum = evidence['PS1'] + evidence['PS3'] + evidence['PS4']
        PM_sum = evidence['PM1'] + evidence['PM2'] + evidence['PM4'] + evidence['PM5']
        PP_sum = evidence['PP2'] + evidence['PP3'] + evidence['PP5']
        
        BA1 = evidence['BA1']
        BS_sum = evidence['BS1'] + evidence['BS2']
        BP_sum = evidence['BP1'] + evidence['BP4'] + evidence['BP6'] + evidence['BP7']
        
        # ACMG分类规则
        classification = "VUS"  # 默认为意义不明
        
        # 致病性判定
        pathogenic = False
        likely_pathogenic = False
        
        # Pathogenic
        if PVS1 == 1:
            if PS_sum >= 1:
                pathogenic = True
            elif PM_sum >= 2:
                pathogenic = True
            elif PM_sum == 1 and PP_sum >= 1:
                pathogenic = True
            elif PP_sum >= 2:
                pathogenic = True
        
        if PS_sum >= 2:
            pathogenic = True
        
        if PS_sum == 1:
            if PM_sum >= 3:
                pathogenic = True
            elif PM_sum == 2 and PP_sum >= 2:
                pathogenic = True
            elif PM_sum == 1 and PP_sum >= 4:
                pathogenic = True
        
        # Likely Pathogenic
        if not pathogenic:
            if PVS1 == 1 and PM_sum == 1:
                likely_pathogenic = True
            elif PS_sum == 1 and (PM_sum == 1 or PM_sum == 2):
                likely_pathogenic = True
            elif PS_sum == 1 and PP_sum >= 2:
                likely_pathogenic = True
            elif PM_sum >= 3:
                likely_pathogenic = True
            elif PM_sum == 2 and PP_sum >= 2:
                likely_pathogenic = True
            elif PM_sum == 1 and PP_sum >= 4:
                likely_pathogenic = True
        
        # 良性判定
        benign = False
        likely_benign = False
        
        if BA1 == 1 or BS_sum >= 2:
            benign = True
        
        if not benign:
            if BS_sum == 1 and BP_sum >= 1:
                likely_benign = True
            elif BP_sum >= 2:
                likely_benign = True
        
        # 最终分类
        if pathogenic and not benign:
            classification = "Pathogenic"
        elif likely_pathogenic and not benign and not likely_benign:
            classification = "Likely_pathogenic"
        elif benign and not pathogenic:
            classification = "Benign"
        elif likely_benign and not pathogenic and not likely_pathogenic:
            classification = "Likely_benign"
        elif (pathogenic or likely_pathogenic) and (benign or likely_benign):
            classification = "Conflicting"
        else:
            classification = "VUS"
        
        return classification, evidence


def format_evidence(evidence):
    """格式化证据为字符串"""
    active = []
    for key, value in evidence.items():
        if value == 1:
            active.append(key)
    return ','.join(active) if active else 'None'


def process_vcf(input_vcf, output_vcf, verbose=False):
    """处理VCF文件，添加ACMG分类"""
    
    classifier = ACMGClassifier(verbose=verbose)
    
    stats = defaultdict(int)
    
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('##'):
                fout.write(line)
            elif line.startswith('#CHROM'):
                # 添加新的INFO字段定义
                fout.write('##INFO=<ID=ACMG_Classification,Number=1,Type=String,Description="ACMG classification: Pathogenic, Likely_pathogenic, VUS, Likely_benign, Benign">\n')
                fout.write('##INFO=<ID=ACMG_Evidence,Number=.,Type=String,Description="ACMG evidence codes">\n')
                fout.write(line)
            else:
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    fout.write(line)
                    continue
                
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                info = fields[7]
                
                # 创建变异对象
                variant = Variant(chrom, pos, ref, alt, info, line.strip())
                
                # 进行ACMG分类
                classification, evidence = classifier.classify(variant)
                
                # 统计
                stats[classification] += 1
                
                # 添加ACMG信息到INFO字段
                evidence_str = format_evidence(evidence)
                new_info = f"{info};ACMG_Classification={classification};ACMG_Evidence={evidence_str}"
                fields[7] = new_info
                
                fout.write('\t'.join(fields) + '\n')
    
    return stats


def main():
    parser = argparse.ArgumentParser(description='HERITA ACMG Classifier - 根据ACMG/AMP指南对变异进行分类')
    parser.add_argument('-i', '--input', required=True, help='输入VCF文件')
    parser.add_argument('-o', '--output', required=True, help='输出VCF文件')
    parser.add_argument('-r', '--report', help='输出报告文件')
    parser.add_argument('--intervar-db', help='InterVar知识库路径')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细输出')
    
    args = parser.parse_args()
    
    # 设置InterVar知识库路径
    global INTERVAR_DB
    if args.intervar_db:
        INTERVAR_DB = args.intervar_db
    
    print("=" * 60, file=sys.stderr)
    print("HERITA ACMG Classifier", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"输入文件: {args.input}", file=sys.stderr)
    print(f"输出文件: {args.output}", file=sys.stderr)
    print(f"知识库路径: {INTERVAR_DB}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    # 加载知识库
    print("\n加载InterVar知识库...", file=sys.stderr)
    load_knowledge_bases()
    
    # 处理VCF
    print("\n处理VCF文件...", file=sys.stderr)
    stats = process_vcf(args.input, args.output, args.verbose)
    
    # 输出统计
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("HERITA ACMG Classification Report")
    report_lines.append("=" * 60)
    report_lines.append(f"输入文件: {args.input}")
    report_lines.append(f"输出文件: {args.output}")
    report_lines.append("")
    report_lines.append("ACMG分类统计")
    report_lines.append("-" * 40)
    total = sum(stats.values())
    for cls in ['Pathogenic', 'Likely_pathogenic', 'VUS', 'Likely_benign', 'Benign', 'Conflicting']:
        count = stats.get(cls, 0)
        pct = count / total * 100 if total > 0 else 0
        report_lines.append(f"  {cls}: {count} ({pct:.1f}%)")
    report_lines.append(f"  总计: {total}")
    report_lines.append("=" * 60)
    
    # 打印到stderr
    for line in report_lines:
        print(line, file=sys.stderr)
    
    # 写入报告文件
    if args.report:
        with open(args.report, 'w') as f:
            f.write('\n'.join(report_lines))
        print(f"\n报告已写入: {args.report}", file=sys.stderr)
    
    print(f"\n完成! 输出已写入: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
