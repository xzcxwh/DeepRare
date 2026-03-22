#!/usr/bin/env python3
"""
HERITA ACMG Classifier v2
根据ACMG/AMP 2015指南对变异进行致病性分类
参考InterVar的评级逻辑和知识库

输出4个字段:
- ACMG_InterVar: 分类结果 (Pathogenic/Likely_pathogenic/Uncertain_significance/Likely_benign/Benign)
- ACMG_PathEvidence: 致病证据代码
- ACMG_BenEvidence: 良性证据代码
- ACMG_Reason: 分类原因
"""

import sys
import os
import re
import argparse
from datetime import datetime
from collections import defaultdict

# 默认知识库路径
DEFAULT_INTERVAR_DB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                                    "InterVar-master", "intervardb")

# 全局知识库
knowledge_base = {}

# ============================================================================
# 加载 InterVar 知识库
# ============================================================================

def load_knowledge_bases(intervar_db):
    """加载所有InterVar知识库"""
    global knowledge_base
    
    print("加载InterVar知识库...", file=sys.stderr)
    
    # PVS1: LOF不耐受基因
    lof_genes = set()
    lof_file = os.path.join(intervar_db, "PVS1.LOF.genes.hg38")
    if os.path.exists(lof_file):
        with open(lof_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    lof_genes.add(gene.upper())
    print(f"  LOF genes: {len(lof_genes)}", file=sys.stderr)
    
    # PS1: 已知致病氨基酸改变
    ps1_changes = {}
    ps1_file = os.path.join(intervar_db, "PS1.AA.change.patho.hg38")
    if os.path.exists(ps1_file):
        with open(ps1_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    chrom, start, end, ref, alt = parts[:5]
                    aa_change = parts[7] if len(parts) > 7 else ''
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"
                    key = f"{chrom}_{start}_{ref}_{alt}"
                    ps1_changes[key] = aa_change
    print(f"  PS1 AA changes: {len(ps1_changes)}", file=sys.stderr)
    
    # PM1: 功能域
    pm1_domains = {}
    pm1_file = os.path.join(intervar_db, "PM1_domains_with_benigns.hg38")
    if os.path.exists(pm1_file):
        with open(pm1_file, 'r') as f:
            next(f, None)  # 跳过header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom, gene, domain = parts[:3]
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"
                    key = f"{chrom}_{gene.upper()}"
                    pm1_domains[key] = domain
    print(f"  PM1 domains: {len(pm1_domains)}", file=sys.stderr)
    
    # PP2: 错义致病基因
    pp2_genes = set()
    pp2_file = os.path.join(intervar_db, "PP2.genes.hg38")
    if os.path.exists(pp2_file):
        with open(pp2_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    pp2_genes.add(gene.upper())
    print(f"  PP2 genes: {len(pp2_genes)}", file=sys.stderr)
    
    # BP1: 截断良性基因
    bp1_genes = set()
    bp1_file = os.path.join(intervar_db, "BP1.genes.hg38")
    if os.path.exists(bp1_file):
        with open(bp1_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    bp1_genes.add(gene.upper())
    print(f"  BP1 genes: {len(bp1_genes)}", file=sys.stderr)
    
    # PS4: 病例对照显著变异
    ps4_variants = set()
    ps4_file = os.path.join(intervar_db, "PS4.variants.hg38")
    if os.path.exists(ps4_file):
        with open(ps4_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom, pos, ref, alt = parts[:4]
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"
                    key = f"{chrom}_{pos}_{ref}_{alt}"
                    ps4_variants.add(key)
    print(f"  PS4 variants: {len(ps4_variants)}", file=sys.stderr)
    
    # BS2: 健康人群变异
    bs2_variants = set()
    bs2_file = os.path.join(intervar_db, "BS2_hom_het.hg38")
    if os.path.exists(bs2_file):
        try:
            with open(bs2_file, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        chrom, pos, ref, alt = parts[:4]
                        if not chrom.startswith('chr'):
                            chrom = f"chr{chrom}"
                        key = f"{chrom}_{pos}_{ref}_{alt}"
                        bs2_variants.add(key)
        except:
            pass
    print(f"  BS2 variants: {len(bs2_variants)}", file=sys.stderr)
    
    knowledge_base = {
        'lof_genes': lof_genes,
        'ps1_changes': ps1_changes,
        'pm1_domains': pm1_domains,
        'pp2_genes': pp2_genes,
        'bp1_genes': bp1_genes,
        'ps4_variants': ps4_variants,
        'bs2_variants': bs2_variants,
    }
    
    return knowledge_base

# ============================================================================
# INFO字段解析
# ============================================================================

def parse_info(info_str):
    """解析VCF INFO字段"""
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info

def get_float(info, keys, default=None):
    """安全获取浮点数，支持多个候选key，处理多值字段"""
    if isinstance(keys, str):
        keys = [keys]
    
    for key in keys:
        val = info.get(key, '.')
        if val not in ['.', '', 'NA', 'None', '-']:
            try:
                val_str = str(val)
                # 处理多值情况 - 找到第一个有效的数值
                if '\n' in val_str or '\t' in val_str:
                    # 多行值，取第一个有效数值
                    for v in val_str.replace('\n', ',').replace('\t', ',').split(','):
                        v = v.strip()
                        if v and v not in ['.', '', 'NA', 'None', '-']:
                            try:
                                return float(v)
                            except:
                                continue
                elif ',' in val_str:
                    for v in val_str.split(','):
                        v = v.strip()
                        if v and v not in ['.', '', 'NA', 'None', '-']:
                            try:
                                return float(v)
                            except:
                                continue
                else:
                    return float(val_str)
            except:
                continue
    return default

def get_str(info, keys, default='.'):
    """安全获取字符串，支持多个候选key"""
    if isinstance(keys, str):
        keys = [keys]
    
    for key in keys:
        val = info.get(key, '')
        if val and val not in ['.', '', 'NA', 'None']:
            return val
    return default

def get_gene_from_info(info):
    """从INFO中提取基因名"""
    # 尝试多种字段 - 优先使用明确的字段
    gene = get_str(info, ['genename', 'VEP_Gene', 'Gene', 'SYMBOL'], '.')
    
    # 如果基因名包含多个值（换行或逗号），取第一个
    if gene != '.':
        gene = gene.split('\n')[0].split(',')[0].split('\t')[0].strip()
    
    # 从ANN字段提取
    if gene == '.' and 'ANN' in info:
        ann = info['ANN']
        # ANN格式: ALT|Effect|Impact|Gene|...
        parts = ann.split('|')
        if len(parts) >= 4:
            gene = parts[3]
    
    return gene.upper() if gene != '.' else '.'

def get_effect_from_info(info):
    """从INFO中提取变异效果"""
    # 从ANN字段提取
    if 'ANN' in info:
        ann = info['ANN']
        parts = ann.split('|')
        if len(parts) >= 2:
            return parts[1].lower()
    
    # 从VEP_Consequence
    consequence = get_str(info, ['VEP_Consequence', 'Consequence'], '.')
    return consequence.lower() if consequence != '.' else ''

def get_aa_change(info):
    """获取氨基酸改变"""
    aaref = get_str(info, ['aaref', 'VEP_AAref'], '.')
    aaalt = get_str(info, ['aaalt', 'VEP_AAalt'], '.')
    
    # 处理多值字段
    if aaref != '.':
        aaref = aaref.split('\n')[0].split(',')[0].split('\t')[0].strip()
    if aaalt != '.':
        aaalt = aaalt.split('\n')[0].split(',')[0].split('\t')[0].strip()
    
    return aaref, aaalt

# ============================================================================
# ACMG 证据评估
# ============================================================================

def evaluate_evidence(chrom, pos, ref, alt, info):
    """评估单个变异的ACMG证据"""
    
    var_key = f"{chrom}_{pos}_{ref}_{alt}"
    gene = get_gene_from_info(info)
    effect = get_effect_from_info(info)
    aaref, aaalt = get_aa_change(info)
    
    # 获取各种注释分数
    af = get_float(info, ['gnomAD4.1_joint_POPMAX_AF', 'gnomAD_AF', 'AF', 'gnomAD_AF_grpmax'])
    clinvar_sig = get_str(info, ['clinvar_clnsig', 'ClinVar_CLNSIG', 'CLNSIG'], '.').lower()
    
    # 获取 ClinVar review status - 检查多种字段名
    clinvar_review = get_str(info, ['clinvar_review', 'ClinVar_CLNREVSTAT', 'CLNREVSTAT'], '.').lower()
    
    revel = get_float(info, ['REVEL_score', 'REVEL', 'VEP_REVEL'])
    cadd = get_float(info, ['CADD_phred', 'CADD', 'VEP_CADD'])
    metarnn = get_float(info, ['MetaRNN_score', 'MetaRNN'])
    alphamissense = get_float(info, ['AlphaMissense_score', 'AlphaMissense', 'VEP_AlphaMissense'])
    clinpred = get_float(info, ['ClinPred_score', 'ClinPred'])
    phylop = get_float(info, ['phyloP100way_vertebrate', 'phyloP', 'VEP_PhyloP'])
    primateai = get_float(info, ['PrimateAI_score', 'PrimateAI'])
    esm1b = get_float(info, ['ESM1b_score', 'ESM1b'])
    
    # SpliceAI
    spliceai_max = 0
    for key in ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL']:
        val = get_float(info, [key])
        if val is not None and val > spliceai_max:
            spliceai_max = val
    
    # 从VEP_SpliceAI_Details解析
    spliceai_details = get_str(info, ['VEP_SpliceAI_Details'], '.')
    if spliceai_details != '.':
        for part in spliceai_details.split(','):
            if '=' in part:
                try:
                    val = float(part.split('=')[1])
                    if val > spliceai_max:
                        spliceai_max = val
                except:
                    pass
    
    pathogenic = []
    benign = []
    details = []
    
    # ========================================================================
    # 致病证据
    # ========================================================================
    
    # PVS1: LOF变异在LOF不耐受基因中
    # 检测方式1: 通过effect字段
    lof_effects = ['frameshift', 'stop_gained', 'splice_donor', 'splice_acceptor', 
                   'start_lost', 'nonsense', 'stop_lost']
    is_lof_by_effect = any(eff in effect for eff in lof_effects)
    
    # 检测方式2: 通过氨基酸改变 (aaalt=X 表示终止密码子/无义突变)
    is_lof_by_aa = (aaalt.upper() == 'X' or aaalt == '*')
    
    is_lof = is_lof_by_effect or is_lof_by_aa
    
    if is_lof and gene in knowledge_base['lof_genes']:
        pathogenic.append('PVS1')
        if is_lof_by_aa:
            details.append(f"PVS1: 无义突变(aaalt={aaalt})在LOF不耐受基因{gene}")
        else:
            details.append(f"PVS1: LOF变异({effect})在LOF不耐受基因{gene}")
    
    # PS1: 与已知致病变异相同的氨基酸改变
    if var_key in knowledge_base['ps1_changes']:
        pathogenic.append('PS1')
        details.append("PS1: 与已知致病变异相同的氨基酸改变")
    
    # PS3: ClinVar专家评审 (从ClinVar注释推断)
    # 检测 "reviewed_by_expert_panel" 等关键词
    has_expert_review = ('expert' in clinvar_review or 
                         'reviewed_by_expert' in clinvar_review or
                         'practice_guideline' in clinvar_review)
    
    if 'pathogenic' in clinvar_sig and 'conflicting' not in clinvar_sig:
        if has_expert_review:
            pathogenic.append('PS3')
            details.append("PS3: ClinVar专家评审的致病证据")
    
    # PS4: 病例对照研究显著
    if var_key in knowledge_base['ps4_variants']:
        pathogenic.append('PS4')
        details.append("PS4: 病例对照研究显著关联")
    
    # PM1: 位于突变热点或功能域
    # 如果基因在pm1_domains知识库中，且是错义变异，则加PM1
    is_missense = ('missense' in effect or 
                   (aaref != '.' and aaalt != '.' and aaalt.upper() != 'X' and aaref != aaalt))
    
    domain_key = f"{chrom}_{gene}"
    # 检查该基因是否在功能域列表中
    gene_in_domains = any(gene in key for key in knowledge_base['pm1_domains'].keys())
    
    if is_missense and gene_in_domains:
        pathogenic.append('PM1')
        details.append(f"PM1: 位于{gene}功能域")
    # 如果基因在LOF列表但不在良性域列表中，也可能是PM1
    elif is_missense and gene in knowledge_base['lof_genes']:
        pathogenic.append('PM1')
        details.append(f"PM1: 位于{gene}功能重要区域")
    
    # PM2: 人群频率极低
    if af is None:
        pathogenic.append('PM2')
        details.append("PM2: 在gnomAD中未发现")
    elif af < 0.0001:
        pathogenic.append('PM2')
        details.append(f"PM2: 人群频率极低(AF={af:.2e})")
    
    # PM4: 框内缺失/插入或终止密码子丢失
    if 'inframe_deletion' in effect or 'inframe_insertion' in effect:
        pathogenic.append('PM4')
        details.append("PM4: 框内缺失/插入")
    if 'stop_lost' in effect:
        if 'PM4' not in pathogenic:
            pathogenic.append('PM4')
        details.append("PM4: 终止密码子丢失")
    
    # PM5: 与已知致病变异相同位置的新错义突变
    for known_key in knowledge_base['ps1_changes']:
        if known_key.startswith(f"{chrom}_{pos}_") and known_key != var_key:
            pathogenic.append('PM5')
            details.append("PM5: 与已知致病变异相同位置的新错义突变")
            break
    
    # PP2: 错义变异在错义致病基因中
    if is_missense and gene in knowledge_base['pp2_genes']:
        pathogenic.append('PP2')
        details.append(f"PP2: {gene}基因中错义突变常致病")
    
    # PP3: 多个计算预测有害
    comp_damaging = 0
    comp_details = []
    
    if revel is not None and revel > 0.5:
        comp_damaging += 1
        comp_details.append(f"REVEL={revel:.3f}")
    if cadd is not None and cadd > 20:
        comp_damaging += 1
        comp_details.append(f"CADD={cadd:.1f}")
    if metarnn is not None and metarnn > 0.5:
        comp_damaging += 1
        comp_details.append(f"MetaRNN={metarnn:.3f}")
    if alphamissense is not None and alphamissense > 0.564:
        comp_damaging += 1
        comp_details.append(f"AlphaMissense={alphamissense:.3f}")
    if clinpred is not None and clinpred > 0.5:
        comp_damaging += 1
        comp_details.append(f"ClinPred={clinpred:.3f}")
    if phylop is not None and phylop > 2.5:
        comp_damaging += 1
        comp_details.append(f"phyloP={phylop:.2f}")
    if primateai is not None and primateai > 0.7:
        comp_damaging += 1
        comp_details.append(f"PrimateAI={primateai:.3f}")
    if esm1b is not None and esm1b < -7:
        comp_damaging += 1
        comp_details.append(f"ESM1b={esm1b:.2f}")
    if spliceai_max > 0.5:
        comp_damaging += 1
        comp_details.append(f"SpliceAI={spliceai_max:.2f}")
    
    if comp_damaging >= 3:
        pathogenic.append('PP3')
        details.append(f"PP3: 多个计算预测有害({', '.join(comp_details)})")
    
    # PP5: ClinVar报告致病
    if 'pathogenic' in clinvar_sig or 'likely_pathogenic' in clinvar_sig:
        if 'conflicting' not in clinvar_sig:
            pathogenic.append('PP5')
            details.append(f"PP5: ClinVar报告为致病")
    
    # ========================================================================
    # 良性证据
    # ========================================================================
    
    # BA1: 人群频率 > 5%
    if af is not None and af > 0.05:
        benign.append('BA1')
        details.append(f"BA1: 人群频率>5%(AF={af:.4f})")
    
    # BS1: 人群频率 > 1%
    if af is not None and af > 0.01:
        benign.append('BS1')
        details.append(f"BS1: 人群频率较高(AF={af:.4f})")
    
    # BS2: 在健康人群中观察到
    if var_key in knowledge_base['bs2_variants']:
        benign.append('BS2')
        details.append("BS2: 在健康人群中观察到")
    
    # BP1: 截断良性基因中的错义变异
    # 使用 is_missense 判断 (已经在上面计算过)
    if is_missense and gene in knowledge_base['bp1_genes']:
        benign.append('BP1')
        details.append(f"BP1: {gene}基因中错义突变通常良性")
    
    # BP4: 多个计算预测良性
    comp_benign = 0
    ben_details = []
    
    if revel is not None and revel < 0.3:
        comp_benign += 1
        ben_details.append(f"REVEL={revel:.3f}")
    if cadd is not None and cadd < 15:
        comp_benign += 1
        ben_details.append(f"CADD={cadd:.1f}")
    if metarnn is not None and metarnn < 0.3:
        comp_benign += 1
        ben_details.append(f"MetaRNN={metarnn:.3f}")
    if alphamissense is not None and alphamissense < 0.34:
        comp_benign += 1
        ben_details.append(f"AlphaMissense={alphamissense:.3f}")
    
    if comp_benign >= 3:
        benign.append('BP4')
        details.append(f"BP4: 多个计算预测良性({', '.join(ben_details)})")
    
    # BP6: ClinVar报告良性
    if 'benign' in clinvar_sig or 'likely_benign' in clinvar_sig:
        if 'conflicting' not in clinvar_sig:
            benign.append('BP6')
            details.append(f"BP6: ClinVar报告为良性")
    
    # BP7: 同义变异且不影响剪接
    if 'synonymous' in effect and spliceai_max < 0.2:
        benign.append('BP7')
        details.append("BP7: 同义变异且不影响剪接")
    
    return pathogenic, benign, details

# ============================================================================
# ACMG 分类
# ============================================================================

def classify_acmg(pathogenic_evidence, benign_evidence):
    """根据ACMG规则进行分类"""
    
    # 计数各级别证据
    pvs = sum(1 for e in pathogenic_evidence if e.startswith('PVS'))
    ps = sum(1 for e in pathogenic_evidence if e.startswith('PS'))
    pm = sum(1 for e in pathogenic_evidence if e.startswith('PM'))
    pp = sum(1 for e in pathogenic_evidence if e.startswith('PP'))
    
    ba = sum(1 for e in benign_evidence if e.startswith('BA'))
    bs = sum(1 for e in benign_evidence if e.startswith('BS'))
    bp = sum(1 for e in benign_evidence if e.startswith('BP'))
    
    # Stand-alone Benign
    if ba >= 1:
        return "Benign", "BA1独立证据"
    
    # Benign
    if bs >= 2:
        return "Benign", f"{bs}个强良性证据"
    
    # Likely Benign
    if bs >= 1 and bp >= 1:
        return "Likely_benign", f"{bs}强+{bp}支持良性证据"
    if bp >= 2:
        return "Likely_benign", f"{bp}个支持良性证据"
    
    # Pathogenic
    if pvs >= 1:
        if ps >= 1:
            return "Pathogenic", f"PVS1+{ps}个强致病证据"
        if pm >= 2:
            return "Pathogenic", f"PVS1+{pm}个中等致病证据"
        if pm >= 1 and pp >= 1:
            return "Pathogenic", f"PVS1+{pm}中等+{pp}支持证据"
        if pp >= 2:
            return "Pathogenic", f"PVS1+{pp}个支持证据"
    
    if ps >= 2:
        return "Pathogenic", f"{ps}个强致病证据"
    if ps >= 1 and pm >= 3:
        return "Pathogenic", f"{ps}强+{pm}中等证据"
    if ps >= 1 and pm >= 2 and pp >= 2:
        return "Pathogenic", f"{ps}强+{pm}中等+{pp}支持证据"
    if ps >= 1 and pm >= 1 and pp >= 4:
        return "Pathogenic", f"{ps}强+{pm}中等+{pp}支持证据"
    
    # Likely Pathogenic
    if pvs >= 1 and pm >= 1:
        return "Likely_pathogenic", f"PVS1+{pm}中等证据"
    if ps >= 1 and pm >= 1:
        return "Likely_pathogenic", f"{ps}强+{pm}中等证据"
    if ps >= 1 and pp >= 2:
        return "Likely_pathogenic", f"{ps}强+{pp}支持证据"
    if pm >= 3:
        return "Likely_pathogenic", f"{pm}个中等证据"
    if pm >= 2 and pp >= 2:
        return "Likely_pathogenic", f"{pm}中等+{pp}支持证据"
    if pm >= 1 and pp >= 4:
        return "Likely_pathogenic", f"{pm}中等+{pp}支持证据"
    
    # VUS - 有一些证据但不足以分类
    if pvs + ps + pm + pp > 0 or bs + bp > 0:
        p_score = pvs * 4 + ps * 3 + pm * 2 + pp * 1
        b_score = bs * 3 + bp * 1
        if p_score > b_score:
            return "Uncertain_significance", f"致病倾向(P-{p_score}_vs_B-{b_score})"
        elif b_score > p_score:
            return "Uncertain_significance", f"良性倾向(P-{p_score}_vs_B-{b_score})"
        else:
            return "Uncertain_significance", f"证据平衡(P-{p_score}_vs_B-{b_score})"
    
    return "Uncertain_significance", "证据不足"

# ============================================================================
# 主处理函数
# ============================================================================

def process_vcf(input_vcf, output_vcf, report_file=None, detailed_log_file=None, verbose=False):
    """处理VCF文件并进行ACMG分类"""
    
    stats = defaultdict(int)
    results = []
    
    # 读取并处理VCF
    header_lines = []
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    chrom = fields[0]
                    pos = fields[1]
                    ref = fields[3]
                    alt = fields[4]
                    info_str = fields[7]
                    
                    info = parse_info(info_str)
                    
                    # 评估证据
                    path_ev, ben_ev, details = evaluate_evidence(chrom, pos, ref, alt, info)
                    
                    # 分类
                    classification, reason = classify_acmg(path_ev, ben_ev)
                    
                    stats[classification] += 1
                    
                    results.append({
                        'line_fields': fields,
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'gene': get_gene_from_info(info),
                        'path_ev': path_ev,
                        'ben_ev': ben_ev,
                        'details': details,
                        'classification': classification,
                        'reason': reason,
                        'info': info,
                    })
    
    print(f"处理了 {len(results)} 个变异", file=sys.stderr)
    
    # 写入输出VCF
    new_info_headers = [
        '##INFO=<ID=ACMG_InterVar,Number=1,Type=String,Description="ACMG classification by InterVar-style analysis">\n',
        '##INFO=<ID=ACMG_PathEvidence,Number=.,Type=String,Description="Pathogenic evidence codes">\n',
        '##INFO=<ID=ACMG_BenEvidence,Number=.,Type=String,Description="Benign evidence codes">\n',
        '##INFO=<ID=ACMG_Reason,Number=1,Type=String,Description="Classification reason">\n',
    ]
    
    with open(output_vcf, 'w') as f:
        for line in header_lines:
            if line.startswith('#CHROM'):
                # 在#CHROM行之前插入新的INFO定义
                for new_header in new_info_headers:
                    f.write(new_header)
            f.write(line)
        
        for result in results:
            fields = result['line_fields']
            path_ev_str = '|'.join(result['path_ev']) if result['path_ev'] else '.'
            ben_ev_str = '|'.join(result['ben_ev']) if result['ben_ev'] else '.'
            reason_str = result['reason'].replace(' ', '_').replace(':', '-')
            
            # 添加ACMG字段到INFO
            new_info = fields[7]
            new_info += f";ACMG_InterVar={result['classification']}"
            new_info += f";ACMG_PathEvidence={path_ev_str}"
            new_info += f";ACMG_BenEvidence={ben_ev_str}"
            new_info += f";ACMG_Reason={reason_str}"
            
            fields[7] = new_info
            f.write('\t'.join(fields) + '\n')
    
    # 写入报告
    if report_file:
        with open(report_file, 'w') as f:
            f.write("# ACMG Classification Report (InterVar-style)\n\n")
            f.write(f"**输入文件**: {input_vcf}\n")
            f.write(f"**输出文件**: {output_vcf}\n\n")
            
            f.write("## 分类统计\n\n")
            f.write("| 分类 | 数量 | 比例 |\n")
            f.write("|------|------|------|\n")
            
            order = ['Pathogenic', 'Likely_pathogenic', 'Uncertain_significance', 'Likely_benign', 'Benign']
            total = len(results)
            for cls in order:
                count = stats.get(cls, 0)
                pct = count / total * 100 if total > 0 else 0
                f.write(f"| {cls} | {count} | {pct:.1f}% |\n")
            
            f.write(f"\n**总计**: {total} 个变异\n\n")
            
            f.write("## 详细结果\n\n")
            
            for cls in order:
                cls_results = [r for r in results if r['classification'] == cls]
                if not cls_results:
                    continue
                
                f.write(f"### {cls} ({len(cls_results)}个)\n\n")
                
                for r in cls_results:
                    f.write(f"**{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}** ({r['gene']})\n")
                    f.write(f"- 致病证据: {', '.join(r['path_ev']) if r['path_ev'] else '无'}\n")
                    f.write(f"- 良性证据: {', '.join(r['ben_ev']) if r['ben_ev'] else '无'}\n")
                    f.write(f"- 分类原因: {r['reason']}\n\n")

    if detailed_log_file:
        write_detailed_log(results, input_vcf, detailed_log_file)
    
    # 打印统计
    print("\n" + "=" * 60, file=sys.stderr)
    print("ACMG分类统计", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    order = ['Pathogenic', 'Likely_pathogenic', 'Uncertain_significance', 'Likely_benign', 'Benign']
    for cls in order:
        count = stats.get(cls, 0)
        pct = count / len(results) * 100 if results else 0
        print(f"  {cls}: {count} ({pct:.1f}%)", file=sys.stderr)
    print(f"  总计: {len(results)}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    return stats

def write_detailed_log(results, input_vcf, log_file):
    """生成详细的ACMG分析日志"""
    log_dir = os.path.dirname(log_file)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)
    with open(log_file, 'w') as f:
        f.write("HERITA Pipeline - Step 10: Detailed ACMG Analysis Log\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input: {input_vcf}\n")
        f.write("=" * 80 + "\n\n")

        for idx, r in enumerate(results, 1):
            info = r.get('info', {})
            gene = r.get('gene', '.')
            f.write(f"Variant #{idx}: {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}\n")
            f.write(f"Gene: {gene}\n")
            f.write("-" * 80 + "\n")
            f.write(f"ACMG Classification: {r['classification']}\n")
            f.write("-" * 80 + "\n")

            f.write("1. Basic Information:\n")
            f.write(f"   - Position: {r['chrom']}:{r['pos']}\n")
            f.write(f"   - Ref/Alt: {r['ref']}/{r['alt']}\n")
            f.write(f"   - Gene: {gene}\n")

            aa_ref = get_str(info, ['aaref'], '')
            aa_pos = get_str(info, ['aapos'], '')
            aa_alt = get_str(info, ['aaalt'], '')
            if not aa_ref and not aa_alt:
                aa_ref, aa_alt = get_aa_change(info)
            aa_change = ""
            if aa_ref or aa_alt:
                aa_change = f"{aa_ref}{aa_pos}{aa_alt}".strip()
            if aa_change:
                f.write(f"   - AA Change: {aa_change}\n")

            f.write("\n2. Pathogenic Evidence:\n")
            if r['path_ev']:
                for ev in sorted(set(r['path_ev'])):
                    ev_details = [d for d in r['details'] if d.startswith(ev)]
                    if ev_details:
                        for detail in ev_details:
                            f.write(f"   [+] {detail}\n")
                    else:
                        f.write(f"   [+] {ev}\n")
            else:
                f.write("   (None)\n")

            f.write("\n3. Benign Evidence:\n")
            if r['ben_ev']:
                for ev in sorted(set(r['ben_ev'])):
                    ev_details = [d for d in r['details'] if d.startswith(ev)]
                    if ev_details:
                        for detail in ev_details:
                            f.write(f"   [-] {detail}\n")
                    else:
                        f.write(f"   [-] {ev}\n")
            else:
                f.write("   (None)\n")

            f.write("\n4. Classification Reason:\n")
            f.write(f"   {r['reason']}\n")

            f.write("\n5. Key Annotations:\n")
            key_fields = [
                ('gnomAD AF', ['gnomAD4.1_joint_POPMAX_AF', 'gnomAD_AF', 'AF', 'gnomAD_AF_grpmax']),
                ('ClinVar', ['clinvar_clnsig', 'ClinVar_CLNSIG', 'CLNSIG']),
                ('ClinVar Review', ['clinvar_review', 'ClinVar_CLNREVSTAT', 'CLNREVSTAT']),
                ('REVEL', ['REVEL_score', 'REVEL', 'VEP_REVEL']),
                ('CADD', ['CADD_phred', 'CADD', 'VEP_CADD']),
                ('MetaRNN', ['MetaRNN_score', 'MetaRNN']),
                ('AlphaMissense', ['AlphaMissense_score', 'AlphaMissense', 'VEP_AlphaMissense']),
                ('ClinPred', ['ClinPred_score', 'ClinPred']),
                ('phyloP', ['phyloP100way_vertebrate', 'phyloP', 'VEP_PhyloP']),
                ('PrimateAI', ['PrimateAI_score', 'PrimateAI']),
                ('ESM1b', ['ESM1b_score', 'ESM1b']),
                ('SpliceAI', ['SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL']),
                ('dbscSNV ADA', ['dbscSNV_ADA_SCORE', 'DBSCSNV_ADA']),
                ('dbscSNV RF', ['dbscSNV_RF_SCORE', 'DBSCSNV_RF']),
            ]
            has_annotation = False
            for label, keys in key_fields:
                val = get_str(info, keys, '.')
                if val != '.':
                    has_annotation = True
                    f.write(f"   - {label}: {val}\n")
            if not has_annotation:
                f.write("   (None)\n")

            f.write("=" * 80 + "\n\n")

def main():
    parser = argparse.ArgumentParser(description='HERITA ACMG Classifier v2 - InterVar风格ACMG分类')
    parser.add_argument('-i', '--input', required=True, help='输入VCF文件')
    parser.add_argument('-o', '--output', required=True, help='输出VCF文件')
    parser.add_argument('-r', '--report', help='输出报告文件')
    parser.add_argument('--intervar-db', default=DEFAULT_INTERVAR_DB, help='InterVar知识库路径')
    parser.add_argument('--detailed-log', help='输出详细ACMG分析日志')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细输出')
    
    args = parser.parse_args()
    
    print("=" * 60, file=sys.stderr)
    print("HERITA ACMG Classifier v2", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"输入文件: {args.input}", file=sys.stderr)
    print(f"输出文件: {args.output}", file=sys.stderr)
    print(f"知识库路径: {args.intervar_db}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.intervar_db):
        print(f"错误: 知识库路径不存在: {args.intervar_db}", file=sys.stderr)
        sys.exit(1)
    
    # 加载知识库
    load_knowledge_bases(args.intervar_db)
    
    # 处理VCF
    print("\n处理VCF文件...", file=sys.stderr)
    process_vcf(
        args.input,
        args.output,
        report_file=args.report,
        detailed_log_file=args.detailed_log,
        verbose=args.verbose,
    )
    
    print(f"\n完成! 输出已写入: {args.output}", file=sys.stderr)
    if args.report:
        print(f"报告已写入: {args.report}", file=sys.stderr)
    if args.detailed_log:
        print(f"详细日志已写入: {args.detailed_log}", file=sys.stderr)

if __name__ == '__main__':
    main()
