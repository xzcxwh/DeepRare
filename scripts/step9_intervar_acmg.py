#!/usr/bin/env python3
"""
Step 9: InterVar-style ACMG Classification
使用 InterVar 的 intervardb 知识库 + 现有 VCF 注释进行 ACMG-AMP 2015 分类

ACMG 证据等级:
- PVS1: 非常强致病证据 (null variant in LOF-intolerant gene)
- PS1-4: 强致病证据
- PM1-6: 中等致病证据  
- PP1-5: 支持致病证据
- BA1: 独立良性证据 (AF > 5%)
- BS1-4: 强良性证据
- BP1-7: 支持良性证据
"""

import os
import sys
import re
import gzip
from collections import defaultdict

# 项目路径
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INTERVARDB_DIR = os.path.join(PROJECT_DIR, "InterVar-master", "intervardb")

# ============================================================================
# 加载 InterVar 知识库
# ============================================================================

def load_lof_genes(build="hg38"):
    """加载 PVS1 LOF不耐受基因列表"""
    filepath = os.path.join(INTERVARDB_DIR, f"PVS1.LOF.genes.{build}")
    genes = set()
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene.upper())
    print(f"  [intervardb] 加载 {len(genes)} 个 LOF 不耐受基因 (PVS1)")
    return genes

def load_ps1_aa_changes(build="hg38"):
    """加载 PS1 已知致病氨基酸改变位点"""
    filepath = os.path.join(INTERVARDB_DIR, f"PS1.AA.change.patho.{build}")
    changes = {}  # key: chr_pos_ref_alt, value: aa_change
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    chrom, start, end, ref, alt, aa_ref, aa_alt, aa_change = parts[:8]
                    # 标准化染色体名
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"
                    key = f"{chrom}_{start}_{ref}_{alt}"
                    changes[key] = aa_change
    print(f"  [intervardb] 加载 {len(changes)} 个已知致病氨基酸改变 (PS1)")
    return changes

def load_pm1_domains(build="hg38"):
    """加载 PM1 功能域（含良性变异的域需排除）"""
    filepath = os.path.join(INTERVARDB_DIR, f"PM1_domains_with_benigns.{build}")
    domains = {}  # key: chr_gene, value: domain
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            next(f, None)  # 跳过header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom, gene, domain = parts[:3]
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"
                    key = f"{chrom}_{gene.upper()}"
                    domains[key] = domain
    print(f"  [intervardb] 加载 {len(domains)} 个功能域记录 (PM1)")
    return domains

def load_pp2_genes(build="hg38"):
    """加载 PP2 错义变异致病基因"""
    filepath = os.path.join(INTERVARDB_DIR, f"PP2.genes.{build}")
    genes = set()
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene.upper())
    print(f"  [intervardb] 加载 {len(genes)} 个错义致病基因 (PP2)")
    return genes

def load_bp1_genes(build="hg38"):
    """加载 BP1 截断变异良性基因"""
    filepath = os.path.join(INTERVARDB_DIR, f"BP1.genes.{build}")
    genes = set()
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene.upper())
    print(f"  [intervardb] 加载 {len(genes)} 个截断良性基因 (BP1)")
    return genes

def load_bs2_variants(build="hg38"):
    """加载 BS2 健康人群纯合/杂合位点"""
    filepath = os.path.join(INTERVARDB_DIR, f"BS2_hom_het.{build}")
    variants = set()
    
    # 尝试多种文件路径
    possible_paths = [filepath, filepath + ".txt"]
    
    for fpath in possible_paths:
        if os.path.exists(fpath):
            try:
                with open(fpath, 'r', encoding='utf-8', errors='ignore') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 4:
                            chrom, pos, ref, alt = parts[:4]
                            if not chrom.startswith('chr'):
                                chrom = f"chr{chrom}"
                            key = f"{chrom}_{pos}_{ref}_{alt}"
                            variants.add(key)
                break
            except Exception as e:
                continue
    
    print(f"  [intervardb] 加载 {len(variants)} 个健康人群变异 (BS2)")
    return variants

def load_ps4_variants(build="hg38"):
    """加载 PS4 病例对照研究显著变异"""
    filepath = os.path.join(INTERVARDB_DIR, f"PS4.variants.{build}")
    variants = set()
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom, pos, ref, alt = parts[:4]
                    if not chrom.startswith('chr'):
                        chrom = f"chr{chrom}"
                    key = f"{chrom}_{pos}_{ref}_{alt}"
                    variants.add(key)
    print(f"  [intervardb] 加载 {len(variants)} 个病例对照显著变异 (PS4)")
    return variants

def load_mim_inheritance():
    """加载 OMIM 遗传模式"""
    dominant = set()
    recessive = set()
    
    dom_file = os.path.join(INTERVARDB_DIR, "mim_domin.txt")
    rec_file = os.path.join(INTERVARDB_DIR, "mim_recessive.txt")
    
    if os.path.exists(dom_file):
        with open(dom_file, 'r') as f:
            for line in f:
                mim = line.strip()
                if mim:
                    dominant.add(mim)
    
    if os.path.exists(rec_file):
        with open(rec_file, 'r') as f:
            for line in f:
                mim = line.strip()
                if mim:
                    recessive.add(mim)
    
    print(f"  [intervardb] 加载 {len(dominant)} 显性 / {len(recessive)} 隐性遗传模式")
    return dominant, recessive

# ============================================================================
# ACMG 证据评估函数
# ============================================================================

def parse_info_field(info_str):
    """解析 VCF INFO 字段"""
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info[key] = value
        else:
            info[item] = True
    return info

def get_float(info, key, default=None):
    """安全获取浮点数"""
    val = info.get(key, '.')
    if val in ['.', '', 'NA', 'None']:
        return default
    try:
        return float(val)
    except:
        return default

def get_str(info, key, default='.'):
    """安全获取字符串"""
    val = info.get(key, default)
    if val in ['.', '', 'NA', 'None']:
        return default
    return val

def evaluate_acmg_evidence(variant, knowledge_base):
    """
    评估单个变异的 ACMG 证据
    返回: (pathogenic_evidence, benign_evidence, details)
    """
    chrom = variant['chrom']
    pos = variant['pos']
    ref = variant['ref']
    alt = variant['alt']
    info = variant['info']
    
    var_key = f"{chrom}_{pos}_{ref}_{alt}"
    gene = get_str(info, 'genename', '.').upper()
    
    # 获取各种注释
    af = get_float(info, 'gnomAD4.1_joint_POPMAX_AF')
    clinvar_sig = get_str(info, 'clinvar_clnsig', '.')
    
    # 功能预测分数
    revel = get_float(info, 'REVEL_score')
    cadd = get_float(info, 'CADD_phred')
    metarnn = get_float(info, 'MetaRNN_score')
    alphamissense = get_float(info, 'AlphaMissense_score')
    clinpred = get_float(info, 'ClinPred_score')
    
    # dbscSNV 剪切预测
    dbscsnv_ada = get_float(info, 'DBSCSNV_ADA')
    dbscsnv_rf = get_float(info, 'DBSCSNV_RF')
    
    # 氨基酸改变
    aaref = get_str(info, 'aaref', '.')
    aaalt = get_str(info, 'aaalt', '.')
    
    # 知识库
    lof_genes = knowledge_base['lof_genes']
    ps1_changes = knowledge_base['ps1_changes']
    pp2_genes = knowledge_base['pp2_genes']
    bp1_genes = knowledge_base['bp1_genes']
    bs2_variants = knowledge_base['bs2_variants']
    ps4_variants = knowledge_base['ps4_variants']
    pm1_domains = knowledge_base['pm1_domains']
    
    pathogenic = []
    benign = []
    details = []
    
    # ========================================================================
    # 致病证据评估
    # ========================================================================
    
    # PVS1: Null variant in LOF-intolerant gene
    # 需要判断是否是无义突变/移码突变 (这里用氨基酸改变近似判断)
    if aaalt == 'X' and gene in lof_genes:  # X 表示终止密码子
        pathogenic.append('PVS1')
        details.append(f"PVS1: 无义突变在LOF不耐受基因 {gene}")
    
    # PS1: Same amino acid change as known pathogenic
    if var_key in ps1_changes:
        pathogenic.append('PS1')
        details.append(f"PS1: 与已知致病变异相同的氨基酸改变")
    
    # PS3: Well-established functional studies (从 ClinVar 推断)
    if 'pathogenic' in clinvar_sig.lower() and 'conflicting' not in clinvar_sig.lower():
        if 'reviewed_by_expert_panel' in get_str(info, 'clinvar_review', ''):
            pathogenic.append('PS3')
            details.append("PS3: ClinVar专家评审的致病证据")
    
    # PS4: Case-control study significant
    if var_key in ps4_variants:
        pathogenic.append('PS4')
        details.append("PS4: 病例对照研究显著关联")
    
    # PM1: Located in mutational hot spot / functional domain
    domain_key = f"{chrom}_{gene}"
    if domain_key not in pm1_domains and gene in lof_genes:
        # 在LOF敏感基因中，且不在已知良性域
        if aaref != '.' and aaalt != '.' and aaref != aaalt:
            pathogenic.append('PM1')
            details.append(f"PM1: 位于 {gene} 功能重要区域")
    
    # PM2: Absent or extremely low frequency in population
    if af is not None and af < 0.0001:  # < 0.01%
        pathogenic.append('PM2')
        details.append(f"PM2: 人群频率极低 (AF={af:.2e})")
    elif af is None:
        pathogenic.append('PM2')
        details.append("PM2: 在gnomAD中未发现")
    
    # PM5: Novel missense at position of known pathogenic missense
    # 简化版：如果有已知致病的同位置不同氨基酸改变
    for known_key in ps1_changes:
        if known_key.startswith(f"{chrom}_{pos}_") and known_key != var_key:
            pathogenic.append('PM5')
            details.append("PM5: 与已知致病变异相同位置的新错义突变")
            break
    
    # PP2: Missense in gene with low benign missense rate
    if gene in pp2_genes and aaref != '.' and aaalt != '.' and aaalt != 'X':
        pathogenic.append('PP2')
        details.append(f"PP2: {gene} 基因中错义突变常致病")
    
    # PP3: Multiple computational evidence supports deleterious
    computational_damaging = 0
    computational_details = []
    
    if revel is not None and revel > 0.5:
        computational_damaging += 1
        computational_details.append(f"REVEL={revel:.3f}")
    if cadd is not None and cadd > 20:
        computational_damaging += 1
        computational_details.append(f"CADD={cadd:.1f}")
    if metarnn is not None and metarnn > 0.5:
        computational_damaging += 1
        computational_details.append(f"MetaRNN={metarnn:.3f}")
    if alphamissense is not None and alphamissense > 0.5:
        computational_damaging += 1
        computational_details.append(f"AlphaMissense={alphamissense:.3f}")
    if clinpred is not None and clinpred > 0.5:
        computational_damaging += 1
        computational_details.append(f"ClinPred={clinpred:.3f}")
    if dbscsnv_ada is not None and dbscsnv_ada > 0.6:
        computational_damaging += 1
        computational_details.append(f"dbscSNV_ADA={dbscsnv_ada:.3f}")
    if dbscsnv_rf is not None and dbscsnv_rf > 0.6:
        computational_damaging += 1
        computational_details.append(f"dbscSNV_RF={dbscsnv_rf:.3f}")
    
    if computational_damaging >= 3:
        pathogenic.append('PP3')
        details.append(f"PP3: 多个计算预测有害 ({', '.join(computational_details)})")
    
    # PP5: Reputable source reports pathogenic (ClinVar)
    if 'pathogenic' in clinvar_sig.lower() or 'likely_pathogenic' in clinvar_sig.lower():
        if 'conflicting' not in clinvar_sig.lower():
            pathogenic.append('PP5')
            details.append(f"PP5: ClinVar报告为 {clinvar_sig}")
    
    # ========================================================================
    # 良性证据评估
    # ========================================================================
    
    # BA1: Allele frequency > 5% in any population
    if af is not None and af > 0.05:
        benign.append('BA1')
        details.append(f"BA1: 人群频率>5% (AF={af:.4f})")
    
    # BS1: Allele frequency greater than expected
    if af is not None and af > 0.01:  # > 1%
        benign.append('BS1')
        details.append(f"BS1: 人群频率较高 (AF={af:.4f})")
    
    # BS2: Observed in healthy individual (dominant) or homozygous (recessive)
    if var_key in bs2_variants:
        benign.append('BS2')
        details.append("BS2: 在健康人群中观察到")
    
    # BP1: Missense in gene where truncating cause disease
    if gene in bp1_genes and aaref != '.' and aaalt != '.' and aaalt != 'X':
        benign.append('BP1')
        details.append(f"BP1: {gene} 基因中错义突变通常良性")
    
    # BP4: Multiple computational evidence supports benign
    computational_benign = 0
    benign_details = []
    
    if revel is not None and revel < 0.3:
        computational_benign += 1
        benign_details.append(f"REVEL={revel:.3f}")
    if cadd is not None and cadd < 15:
        computational_benign += 1
        benign_details.append(f"CADD={cadd:.1f}")
    if metarnn is not None and metarnn < 0.3:
        computational_benign += 1
        benign_details.append(f"MetaRNN={metarnn:.3f}")
    if alphamissense is not None and alphamissense < 0.3:
        computational_benign += 1
        benign_details.append(f"AlphaMissense={alphamissense:.3f}")
    
    if computational_benign >= 3:
        benign.append('BP4')
        details.append(f"BP4: 多个计算预测良性 ({', '.join(benign_details)})")
    
    # BP6: Reputable source reports benign
    if 'benign' in clinvar_sig.lower() or 'likely_benign' in clinvar_sig.lower():
        if 'conflicting' not in clinvar_sig.lower():
            benign.append('BP6')
            details.append(f"BP6: ClinVar报告为 {clinvar_sig}")
    
    return pathogenic, benign, details

def classify_acmg(pathogenic_evidence, benign_evidence):
    """
    根据 ACMG-AMP 2015 规则进行最终分类
    
    Pathogenic:
    - 1 Very Strong (PVS1) + ≥1 Strong (PS1-4) OR
    - 1 Very Strong + ≥2 Moderate (PM1-6) OR
    - 1 Very Strong + 1 Moderate + 1 Supporting (PP1-5) OR
    - 1 Very Strong + ≥2 Supporting OR
    - ≥2 Strong OR
    - 1 Strong + ≥3 Moderate OR
    - 1 Strong + 2 Moderate + ≥2 Supporting OR
    - 1 Strong + 1 Moderate + ≥4 Supporting
    
    Likely Pathogenic:
    - 1 Very Strong + 1 Moderate OR
    - 1 Strong + 1-2 Moderate OR
    - 1 Strong + ≥2 Supporting OR
    - ≥3 Moderate OR
    - 2 Moderate + ≥2 Supporting OR
    - 1 Moderate + ≥4 Supporting
    
    Benign:
    - 1 Stand-alone (BA1) OR
    - ≥2 Strong (BS1-4)
    
    Likely Benign:
    - 1 Strong + 1 Supporting (BP1-7) OR
    - ≥2 Supporting
    """
    
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
            return "Uncertain_significance", f"致病倾向(P:{p_score} vs B:{b_score})"
        elif b_score > p_score:
            return "Uncertain_significance", f"良性倾向(P:{p_score} vs B:{b_score})"
        else:
            return "Uncertain_significance", f"证据平衡(P:{p_score} vs B:{b_score})"
    
    return "Uncertain_significance", "证据不足"

# ============================================================================
# 主处理函数
# ============================================================================

def process_vcf(input_vcf, output_vcf, output_report):
    """处理 VCF 文件并进行 ACMG 分类"""
    
    print("\n" + "="*70)
    print("Step 9: InterVar-style ACMG Classification")
    print("="*70)
    
    # 加载知识库
    print("\n[1/3] 加载 InterVar 知识库...")
    knowledge_base = {
        'lof_genes': load_lof_genes(),
        'ps1_changes': load_ps1_aa_changes(),
        'pm1_domains': load_pm1_domains(),
        'pp2_genes': load_pp2_genes(),
        'bp1_genes': load_bp1_genes(),
        'bs2_variants': load_bs2_variants(),
        'ps4_variants': load_ps4_variants(),
    }
    
    # 读取并处理 VCF
    print("\n[2/3] 评估 ACMG 证据...")
    
    variants = []
    header_lines = []
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # 只保留必要的header行（过滤掉contig等无用信息）
                if line.startswith('##fileformat') or line.startswith('##INFO') or line.startswith('#CHROM'):
                    header_lines.append(line)
            else:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    variant = {
                        'chrom': parts[0],
                        'pos': parts[1],
                        'id': parts[2],
                        'ref': parts[3],
                        'alt': parts[4],
                        'qual': parts[5],
                        'filter': parts[6],
                        'info': parse_info_field(parts[7]),
                        'rest': parts[8:] if len(parts) > 8 else [],
                        'raw_info': parts[7]
                    }
                    variants.append(variant)
    
    print(f"  读取 {len(variants)} 个变异")
    
    # 评估每个变异
    results = []
    classification_counts = defaultdict(int)
    
    for var in variants:
        path_ev, ben_ev, details = evaluate_acmg_evidence(var, knowledge_base)
        classification, reason = classify_acmg(path_ev, ben_ev)
        
        results.append({
            'variant': var,
            'pathogenic_evidence': path_ev,
            'benign_evidence': ben_ev,
            'details': details,
            'classification': classification,
            'reason': reason
        })
        classification_counts[classification] += 1
    
    # 写入输出 VCF
    print("\n[3/3] 写入结果...")
    
    # 添加新的 INFO 字段定义
    new_info_headers = [
        '##INFO=<ID=ACMG_InterVar,Number=1,Type=String,Description="ACMG classification by InterVar-style analysis">\n',
        '##INFO=<ID=ACMG_PathEvidence,Number=.,Type=String,Description="Pathogenic evidence codes">\n',
        '##INFO=<ID=ACMG_BenEvidence,Number=.,Type=String,Description="Benign evidence codes">\n',
        '##INFO=<ID=ACMG_Reason,Number=1,Type=String,Description="Classification reason">\n',
    ]
    
    with open(output_vcf, 'w') as f:
        # 写入 header
        for line in header_lines:
            if line.startswith('#CHROM'):
                # 在 #CHROM 行之前插入新的 INFO 定义
                for new_header in new_info_headers:
                    f.write(new_header)
            f.write(line)
        
        # 写入变异
        for result in results:
            var = result['variant']
            path_ev = '|'.join(result['pathogenic_evidence']) if result['pathogenic_evidence'] else '.'
            ben_ev = '|'.join(result['benign_evidence']) if result['benign_evidence'] else '.'
            classification = result['classification']
            reason = result['reason'].replace(' ', '_').replace(':', '-')
            
            # 构建新的 INFO 字段
            new_info = var['raw_info']
            new_info += f";ACMG_InterVar={classification}"
            new_info += f";ACMG_PathEvidence={path_ev}"
            new_info += f";ACMG_BenEvidence={ben_ev}"
            new_info += f";ACMG_Reason={reason}"
            
            # 写入行
            out_parts = [
                var['chrom'], var['pos'], var['id'], var['ref'], var['alt'],
                var['qual'], var['filter'], new_info
            ] + var['rest']
            f.write('\t'.join(out_parts) + '\n')
    
    # 生成报告
    with open(output_report, 'w') as f:
        f.write("# ACMG Classification Report (InterVar-style)\n\n")
        f.write(f"**分析日期**: 2025年12月3日\n")
        f.write(f"**输入文件**: {input_vcf}\n")
        f.write(f"**输出文件**: {output_vcf}\n\n")
        
        f.write("## 分类统计\n\n")
        f.write("| 分类 | 数量 | 比例 |\n")
        f.write("|------|------|------|\n")
        
        order = ['Pathogenic', 'Likely_pathogenic', 'Uncertain_significance', 'Likely_benign', 'Benign']
        for cls in order:
            count = classification_counts.get(cls, 0)
            pct = count / len(results) * 100 if results else 0
            f.write(f"| {cls} | {count} | {pct:.1f}% |\n")
        
        f.write(f"\n**总计**: {len(results)} 个变异\n\n")
        
        f.write("## 详细结果\n\n")
        
        # 按分类分组输出
        for cls in order:
            cls_results = [r for r in results if r['classification'] == cls]
            if not cls_results:
                continue
            
            f.write(f"### {cls} ({len(cls_results)}个)\n\n")
            
            for r in cls_results:
                var = r['variant']
                gene = get_str(var['info'], 'genename', '.')
                clinvar = get_str(var['info'], 'clinvar_clnsig', '.')
                
                f.write(f"**{var['chrom']}:{var['pos']} {var['ref']}>{var['alt']}** ({gene})\n")
                f.write(f"- ClinVar: {clinvar}\n")
                f.write(f"- 致病证据: {', '.join(r['pathogenic_evidence']) if r['pathogenic_evidence'] else '无'}\n")
                f.write(f"- 良性证据: {', '.join(r['benign_evidence']) if r['benign_evidence'] else '无'}\n")
                f.write(f"- 分类原因: {r['reason']}\n")
                if r['details']:
                    f.write(f"- 证据详情:\n")
                    for detail in r['details']:
                        f.write(f"  - {detail}\n")
                f.write("\n")
    
    # 打印统计
    print("\n" + "-"*50)
    print("ACMG 分类结果统计:")
    print("-"*50)
    for cls in order:
        count = classification_counts.get(cls, 0)
        pct = count / len(results) * 100 if results else 0
        bar = '█' * int(pct / 2)
        print(f"  {cls:30s}: {count:3d} ({pct:5.1f}%) {bar}")
    print("-"*50)
    print(f"  总计: {len(results)} 个变异")
    
    print(f"\n输出文件:")
    print(f"  - VCF: {output_vcf}")
    print(f"  - 报告: {output_report}")
    
    return classification_counts

def main():
    # 默认路径
    input_vcf = os.path.join(PROJECT_DIR, "results", "final_integrated.vcf")
    output_vcf = os.path.join(PROJECT_DIR, "results", "final_intervar_classified.vcf")
    output_report = os.path.join(PROJECT_DIR, "results", "STEP9_INTERVAR_ACMG_REPORT.md")
    
    # 命令行参数
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_vcf = sys.argv[2]
    if len(sys.argv) > 3:
        output_report = sys.argv[3]
    
    if not os.path.exists(input_vcf):
        print(f"错误: 输入文件不存在: {input_vcf}")
        sys.exit(1)
    
    process_vcf(input_vcf, output_vcf, output_report)
    print("\n✓ Step 9 完成!")

if __name__ == "__main__":
    main()
