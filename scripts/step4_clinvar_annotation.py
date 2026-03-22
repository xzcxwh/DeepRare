#!/usr/bin/env python3
"""
VCF ClinVar注释处理脚本 - 优化版
使用tabix按需查询，而不是加载整个dbNSFP数据库

优化说明：
- 原版本加载整个dbNSFP(1.1GB压缩, 14M+记录)到内存，耗时约25秒
- 优化版先收集vcf3中的变异位点（约276个），再用tabix按需查询
- 预计耗时：2-5秒
"""

import gzip
import os
import subprocess

# 文件路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

DBNSFP_FILE = os.path.join(PROJECT_ROOT, "input_data", "dbNSFP5.3a_grch38_lite.tsv.gz")
VCF3_PATH = os.path.join(PROJECT_ROOT, "results", "vcf3.vcf")
FINAL_CLINVAR_PATH = os.path.join(PROJECT_ROOT, "results", "final_clinvar.vcf")
VCF4_PATH = os.path.join(PROJECT_ROOT, "results", "vcf4.vcf")

# 缓存dbNSFP header
_dbnsfp_header = None
_dbnsfp_indices = None


def get_dbnsfp_header():
    """获取并缓存dbNSFP header和列索引"""
    global _dbnsfp_header, _dbnsfp_indices
    if _dbnsfp_header is None:
        with gzip.open(DBNSFP_FILE, 'rt') as f:
            _dbnsfp_header = f.readline().strip().split('\t')
        _dbnsfp_indices = {
            'chr': _dbnsfp_header.index('#chr'),
            'pos': _dbnsfp_header.index('pos(1-based)'),
            'ref': _dbnsfp_header.index('ref'),
            'alt': _dbnsfp_header.index('alt'),
            'clinvar_clnsig': _dbnsfp_header.index('clinvar_clnsig'),
            'clinvar_review': _dbnsfp_header.index('clinvar_review')
        }
    return _dbnsfp_header, _dbnsfp_indices


def query_clinvar_tabix(chrom, pos, ref, alt):
    """
    使用tabix快速查询dbNSFP获取ClinVar注释
    返回: (clinvar_clnsig, clinvar_review) 或 ('.', '.')
    """
    # dbNSFP使用不带chr前缀的染色体
    chrom_db = chrom[3:] if chrom.startswith('chr') else chrom
    region = f"{chrom_db}:{pos}-{pos}"
    
    try:
        result = subprocess.run(
            ['tabix', DBNSFP_FILE, region],
            capture_output=True, text=True, timeout=5
        )
        
        if result.returncode != 0 or not result.stdout.strip():
            return '.', '.'
        
        header, indices = get_dbnsfp_header()
        
        # 解析结果，查找匹配的ref和alt
        for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')
            if len(fields) <= max(indices['clinvar_clnsig'], indices['clinvar_review']):
                continue
            
            # 检查ref和alt是否匹配
            if fields[indices['ref']] == ref and fields[indices['alt']] == alt:
                clnsig = fields[indices['clinvar_clnsig']]
                review = fields[indices['clinvar_review']]
                return clnsig if clnsig else '.', review if review else '.'
        
        return '.', '.'
        
    except Exception as e:
        return '.', '.'


def should_be_in_final_clinvar(clnsig, review):
    """
    判断是否应该进入final_clinvar.vcf
    条件（满足任一即可）：
    1. pathogenic变异且review为practice_guideline或reviewed_by_expert_panel
    2. pathogenic变异且review包含multiple_submitters且no_conflicts
    """
    if clnsig == '.' or review == '.':
        return False

    pathogenic_terms = [
        'Pathogenic',
        'Pathogenic/Likely_pathogenic',
        'Pathogenic|drug_response'
    ]

    # 条件1：专家审核
    expert_reviews = [
        'practice_guideline',
        'reviewed_by_expert_panel'
    ]

    # 条件2：多个提交者且无冲突
    multiple_no_conflict = (
        'multiple_submitters' in review and 
        'no_conflicts' in review
    )

    is_pathogenic = clnsig in pathogenic_terms
    has_expert_review = any(r in review for r in expert_reviews)

    return is_pathogenic and (has_expert_review or multiple_no_conflict)


def should_be_deleted_from_vcf4(clnsig, review):
    """
    判断是否应该从vcf4中删除
    benign变异且review为practice_guideline或reviewed_by_expert_panel
    """
    if clnsig == '.' or review == '.':
        return False

    benign_terms = [
        'Likely_benign',
        'Benign',
        'Benign/Likely_benign'
    ]

    expert_reviews = [
        'practice_guideline',
        'reviewed_by_expert_panel'
    ]

    return (clnsig in benign_terms and
            any(r in review for r in expert_reviews))


def process_vcf3():
    """处理vcf3文件（优化版 - 使用tabix按需查询）"""
    import time
    start_time = time.time()
    
    print("🧬 ClinVar注释 - 优化版 (使用tabix按需查询)")
    print("=" * 60)
    
    # 检查tabix是否可用
    try:
        subprocess.run(['tabix', '--version'], capture_output=True, check=True)
    except Exception:
        print("❌ tabix未安装或不可用，请安装htslib: brew install htslib")
        return
    
    # 检查输入文件
    if not os.path.exists(VCF3_PATH):
        print(f"❌ 输入文件不存在: {VCF3_PATH}")
        return
    
    # 统计vcf3变异数量
    vcf3_count = 0
    with open(VCF3_PATH, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                vcf3_count += 1
    print(f"📊 vcf3变异数量: {vcf3_count}")
    print()
    
    final_clinvar_count = 0
    vcf4_count = 0
    deleted_count = 0
    total_processed = 0

    with open(VCF3_PATH, 'r') as infile, \
         open(FINAL_CLINVAR_PATH, 'w') as final_out, \
         open(VCF4_PATH, 'w') as vcf4_out:

        # 写入VCF头部
        header_lines = []
        for line in infile:
            if line.startswith('#'):
                header_lines.append(line)
                final_out.write(line)
                vcf4_out.write(line)
            else:
                # 处理完头部，回到第一个变异行
                infile.seek(0)
                # 跳过头部
                for _ in header_lines:
                    next(infile)
                break

        # 处理变异行
        for line in infile:
            if line.startswith('#') or line.strip() == '':
                continue

            total_processed += 1
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            # 处理多个ALT等位基因
            alts = alt.split(',')
            for single_alt in alts:
                # 使用tabix查询ClinVar注释
                clnsig, review = query_clinvar_tabix(chrom, pos, ref, single_alt)

                # 添加ClinVar注释到INFO字段
                updated_info = info
                if updated_info == '.':
                    updated_info = f'CLINVAR_CLNSIG={clnsig};CLINVAR_REVIEW={review}'
                else:
                    updated_info += f';CLINVAR_CLNSIG={clnsig};CLINVAR_REVIEW={review}'

                # 更新字段
                updated_fields = fields.copy()
                updated_fields[4] = single_alt
                updated_fields[7] = updated_info
                updated_line = '\t'.join(updated_fields) + '\n'

                # 判断分类
                if should_be_in_final_clinvar(clnsig, review):
                    final_out.write(updated_line)
                    final_clinvar_count += 1
                    print(f"变异 {chrom}:{pos} {ref}>{single_alt} 进入决赛圈: {clnsig}, {review}")

                elif should_be_deleted_from_vcf4(clnsig, review):
                    deleted_count += 1
                    print(f"变异 {chrom}:{pos} {ref}>{single_alt} 被删除: {clnsig}, {review}")
                else:
                    vcf4_out.write(updated_line)
                    vcf4_count += 1

    elapsed = time.time() - start_time
    
    print(f"\n处理完成！")
    print(f"总处理变异行数: {total_processed}")
    print(f"进入决赛圈 (final_clinvar.vcf): {final_clinvar_count}")
    print(f"保留在vcf4中: {vcf4_count}")
    print(f"从vcf4中删除: {deleted_count}")
    total_after = final_clinvar_count + vcf4_count + deleted_count
    print(f"验证: {final_clinvar_count} + {vcf4_count} + {deleted_count} = {total_after}")
    print(f"\n⏱️ 总耗时: {elapsed:.2f} 秒")


if __name__ == "__main__":
    process_vcf3()
