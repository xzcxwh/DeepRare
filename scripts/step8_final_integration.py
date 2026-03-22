#!/usr/bin/env python3
"""
步骤8: 整合三个决赛圈文件并补全缺失的注释 (优化版)
- 合并三个决赛圈VCF文件（去重）
- 使用tabix按需查询数据库（而不是加载整个数据库）
- 清理重复的注释值
- 对评分字段取最优值（max或min）

优化说明：
- 原版本加载整个dbNSFP(1.1GB)和dbscSNV(469MB)到内存，耗时约3分钟
- 优化版使用tabix索引按需查询，只查询实际需要的变异位点（约72个）
- 预计耗时：10-20秒
"""

import subprocess
import os
import sys
import logging
import gzip

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 文件路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# 输入文件
CLINVAR_VCF = os.path.join(PROJECT_ROOT, "results", "final_clinvar.vcf")
DBSCSNV_VCF = os.path.join(PROJECT_ROOT, "results", "final_dbscSNV.vcf")
TOP_VCF = os.path.join(PROJECT_ROOT, "results", "final_top.vcf")

# 中间文件和输出文件
MERGED_VCF = os.path.join(PROJECT_ROOT, "results", "merged_final.vcf")
ANNOTATED_VCF = os.path.join(PROJECT_ROOT, "results", "merged_annotated.vcf")
FINAL_VCF = os.path.join(PROJECT_ROOT, "results", "final_integrated.vcf")

# 数据库文件
DBNSFP_FILE = os.path.join(PROJECT_ROOT, "input_data", "dbNSFP5.3a_grch38_lite.tsv.gz")
DBSCSNV_FILE = os.path.join(PROJECT_ROOT, "dbscSNV1", "dbscSNV1.1_hg38_sorted.tsv")

# 定义需要取最大值的字段
MAX_SCORE_FIELDS = [
    'AlphaMissense_score', 'REVEL_score', 'MetaRNN_score', 
    'PrimateAI_score', 'CADD_phred', 'ClinPred_score',
    'DBSCSNV_ADA', 'DBSCSNV_RF'
]

# 定义需要取最小值的字段（ESM1b_score越小表示影响越大）
MIN_SCORE_FIELDS = ['ESM1b_score']

# 定义需要取唯一值（去重）的字段
UNIQUE_FIELDS = [
    'aaref', 'aaalt', 'rs_dbSNP', 'genename',
    'clinvar_id', 'clinvar_clnsig', 'clinvar_trait', 'clinvar_review',
    'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id', 
    'clinvar_OMIM_id', 'clinvar_Orphanet_id',
    'gnomAD4.1_joint', 'gnomAD4.1_joint_POPMAX_AF'
]

# 定义需要保留的INFO字段
FIELDS_TO_KEEP = [
    'gnomAD4.1_joint_POPMAX_AF', 'gnomAD4.1_joint',
    'DBSCSNV_ADA', 'DBSCSNV_RF',
    'MetaRNN_score', 'REVEL_score', 'PrimateAI_score', 'ClinPred_score',
    'ESM1b_score', 'AlphaMissense_score', 'CADD_phred',
    'aaref', 'aaalt', 'rs_dbSNP', 'genename',
    'clinvar_id', 'clinvar_clnsig', 'clinvar_trait', 'clinvar_review',
    'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id',
    'clinvar_OMIM_id', 'clinvar_Orphanet_id'
]

# 缓存dbNSFP header
_dbnsfp_header = None


def check_input_files():
    """检查输入文件是否存在"""
    files_to_check = [CLINVAR_VCF, DBSCSNV_VCF, TOP_VCF, DBNSFP_FILE, DBSCSNV_FILE]
    for file_path in files_to_check:
        if not os.path.exists(file_path):
            logger.error(f"输入文件不存在: {file_path}")
            return False
    return True


def clean_info_value(field_name, value):
    """
    清理INFO字段值，根据字段类型进行不同处理
    支持逗号、分号、竖线等多种分隔符
    """
    if value is None or value == '':
        return '.'
    
    value = str(value).strip()
    
    # 检查是否全是点号
    cleaned_for_check = value.replace(',', '').replace(';', '').replace('|', '').replace('.', '').strip()
    if cleaned_for_check == '':
        return '.'
    
    if value == '.':
        return '.'
    
    # 分割值（支持逗号、分号、竖线）
    if ';' in value:
        values = value.split(';')
    elif ',' in value:
        values = value.split(',')
    elif '|' in value:
        values = value.split('|')
    else:
        if value != '.':
            return value
        return '.'
    
    # 清理每个值
    cleaned_values = []
    for v in values:
        v = v.strip()
        if v and v != '.':
            cleaned_values.append(v)
    
    if not cleaned_values:
        return '.'
    
    # 处理需要取最大值的评分字段
    if field_name in MAX_SCORE_FIELDS:
        numeric_values = []
        for v in cleaned_values:
            try:
                numeric_values.append(float(v))
            except ValueError:
                pass
        if numeric_values:
            return str(max(numeric_values))
        return '.'
    
    # 处理需要取最小值的评分字段
    if field_name in MIN_SCORE_FIELDS:
        numeric_values = []
        for v in cleaned_values:
            try:
                numeric_values.append(float(v))
            except ValueError:
                pass
        if numeric_values:
            return str(min(numeric_values))
        return '.'
    
    # 处理需要取唯一值的字段
    if field_name in UNIQUE_FIELDS or field_name == 'genename':
        unique_values = []
        seen = set()
        for v in cleaned_values:
            if v not in seen:
                unique_values.append(v)
                seen.add(v)
        if unique_values:
            if len(unique_values) == 1:
                return unique_values[0]
            return '|'.join(unique_values)
        return '.'
    
    # 默认：取第一个非空值
    if cleaned_values:
        return cleaned_values[0]
    return '.'


def get_dbnsfp_header():
    """获取并缓存dbNSFP header"""
    global _dbnsfp_header
    if _dbnsfp_header is None:
        with gzip.open(DBNSFP_FILE, 'rt') as f:
            _dbnsfp_header = f.readline().strip().split('\t')
    return _dbnsfp_header


def query_dbnsfp_tabix(chrom, pos, ref, alt):
    """
    使用tabix快速查询dbNSFP数据库
    返回: {field: value, ...} 或 None
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
            return None
        
        header = get_dbnsfp_header()
        
        # 解析结果
        for line in result.stdout.strip().split('\n'):
            fields = line.split('\t')
            if len(fields) < 5:
                continue
            
            # 检查ref和alt是否匹配
            ref_idx = header.index('ref') if 'ref' in header else 2
            alt_idx = header.index('alt') if 'alt' in header else 3
            
            if fields[ref_idx] == ref and fields[alt_idx] == alt:
                # 构建注释字典
                annotations = {}
                field_names = [
                    'aaref', 'aaalt', 'rs_dbSNP', 'genename',
                    'clinvar_id', 'clinvar_clnsig', 'clinvar_trait', 'clinvar_review',
                    'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id', 
                    'clinvar_OMIM_id', 'clinvar_Orphanet_id',
                    'MetaRNN_score', 'REVEL_score', 'PrimateAI_score', 'ClinPred_score',
                    'ESM1b_score', 'AlphaMissense_score', 'CADD_phred',
                    'gnomAD4.1_joint', 'gnomAD4.1_joint_POPMAX_AF'
                ]
                
                for field_name in field_names:
                    if field_name in header:
                        idx = header.index(field_name)
                        if idx < len(fields) and fields[idx] != '.' and fields[idx] != '':
                            annotations[field_name] = fields[idx]
                
                return annotations
        
        return None
        
    except Exception as e:
        logger.debug(f"tabix查询失败 {region}: {e}")
        return None


def load_dbscsnv_for_positions(positions):
    """
    只加载需要的dbscSNV位点
    positions: set of (chrom, pos, ref, alt)
    返回: {(chrom, pos, ref, alt): {'ada_score': value, 'rf_score': value}}
    """
    logger.info(f"从dbscSNV查询 {len(positions)} 个位点...")
    
    # 转换为查找集合（不带chr前缀）
    positions_set = set()
    for chrom, pos, ref, alt in positions:
        chrom_db = chrom[3:] if chrom.startswith('chr') else chrom
        positions_set.add((chrom_db, pos, ref, alt))
    
    lookup = {}
    
    with open(DBSCSNV_FILE, 'r') as f:
        header = f.readline().strip().split('\t')
        
        chr_idx = header.index('chr')
        pos_idx = header.index('pos')
        ref_idx = header.index('ref')
        alt_idx = header.index('alt')
        ada_idx = header.index('ada_score')
        rf_idx = header.index('rf_score')
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) <= max(ada_idx, rf_idx):
                continue
            
            chrom = fields[chr_idx]
            pos = fields[pos_idx]
            ref = fields[ref_idx]
            alt = fields[alt_idx]
            
            key = (chrom, pos, ref, alt)
            if key in positions_set:
                ada_score = fields[ada_idx] if fields[ada_idx] != '.' else None
                rf_score = fields[rf_idx] if fields[rf_idx] != '.' else None
                
                if ada_score or rf_score:
                    lookup[key] = {
                        'ada_score': ada_score,
                        'rf_score': rf_score
                    }
                
                # 移除已找到的位点，减少后续比较
                positions_set.discard(key)
                
                # 如果所有位点都找到了，提前退出
                if not positions_set:
                    break
    
    logger.info(f"dbscSNV匹配到 {len(lookup)} 个位点")
    return lookup


def generate_clean_header():
    """生成干净的VCF header，只包含必要的信息"""
    header_lines = []
    
    # VCF版本
    header_lines.append("##fileformat=VCFv4.2\n")
    
    # 只保留必要的INFO字段定义
    info_definitions = [
        '##INFO=<ID=gnomAD4.1_joint_POPMAX_AF,Number=1,Type=Float,Description="gnomAD 4.1 joint POPMAX allele frequency">',
        '##INFO=<ID=gnomAD4.1_joint,Number=1,Type=Float,Description="gnomAD 4.1 joint allele frequency">',
        '##INFO=<ID=DBSCSNV_ADA,Number=1,Type=Float,Description="dbscSNV ada_score for splice site prediction">',
        '##INFO=<ID=DBSCSNV_RF,Number=1,Type=Float,Description="dbscSNV rf_score for splice site prediction">',
        '##INFO=<ID=MetaRNN_score,Number=1,Type=Float,Description="MetaRNN pathogenicity score">',
        '##INFO=<ID=REVEL_score,Number=1,Type=Float,Description="REVEL pathogenicity score">',
        '##INFO=<ID=PrimateAI_score,Number=1,Type=Float,Description="PrimateAI pathogenicity score">',
        '##INFO=<ID=ClinPred_score,Number=1,Type=Float,Description="ClinPred pathogenicity score">',
        '##INFO=<ID=ESM1b_score,Number=1,Type=Float,Description="ESM1b pathogenicity score (lower = more pathogenic)">',
        '##INFO=<ID=AlphaMissense_score,Number=1,Type=Float,Description="AlphaMissense pathogenicity score">',
        '##INFO=<ID=CADD_phred,Number=1,Type=Float,Description="CADD phred-scaled score">',
        '##INFO=<ID=aaref,Number=1,Type=String,Description="Reference amino acid">',
        '##INFO=<ID=aaalt,Number=1,Type=String,Description="Alternate amino acid">',
        '##INFO=<ID=rs_dbSNP,Number=1,Type=String,Description="dbSNP rsID">',
        '##INFO=<ID=genename,Number=1,Type=String,Description="Gene name">',
        '##INFO=<ID=clinvar_id,Number=1,Type=String,Description="ClinVar variation ID">',
        '##INFO=<ID=clinvar_clnsig,Number=1,Type=String,Description="ClinVar clinical significance">',
        '##INFO=<ID=clinvar_trait,Number=1,Type=String,Description="ClinVar associated trait/disease">',
        '##INFO=<ID=clinvar_review,Number=1,Type=String,Description="ClinVar review status">',
        '##INFO=<ID=clinvar_hgvs,Number=1,Type=String,Description="ClinVar HGVS notation">',
        '##INFO=<ID=clinvar_var_source,Number=1,Type=String,Description="ClinVar variant source">',
        '##INFO=<ID=clinvar_MedGen_id,Number=1,Type=String,Description="ClinVar MedGen ID">',
        '##INFO=<ID=clinvar_OMIM_id,Number=1,Type=String,Description="ClinVar OMIM ID">',
        '##INFO=<ID=clinvar_Orphanet_id,Number=1,Type=String,Description="ClinVar Orphanet ID">',
    ]
    
    for info_def in info_definitions:
        header_lines.append(info_def + "\n")
    
    return header_lines


def merge_vcf_files():
    """合并三个VCF文件（去重），收集需要查询的位点"""
    logger.info("开始合并三个决赛圈VCF文件...")

    all_variants = {}
    positions_to_query = set()
    
    # 读取第一个文件的样本列信息（#CHROM行）
    chrom_header_line = None
    with open(CLINVAR_VCF, 'r') as f:
        for line in f:
            if line.startswith('#CHROM'):
                chrom_header_line = line
                break
    
    # 生成干净的header
    header_lines = generate_clean_header()
    # 添加#CHROM行
    if chrom_header_line:
        header_lines.append(chrom_header_line)

    # 从所有文件读取变异
    for vcf_file in [CLINVAR_VCF, DBSCSNV_VCF, TOP_VCF]:
        logger.info(f"读取文件: {vcf_file}")
        with open(vcf_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if len(fields) >= 8:
                        chrom = fields[0]
                        pos = fields[1]
                        var_id = fields[2]
                        ref = fields[3]
                        alt = fields[4]
                        qual = fields[5]
                        filt = fields[6]
                        info = fields[7]
                        
                        format_samples = '\t'.join(fields[8:]) if len(fields) > 8 else ''
                        variant_key = f"{chrom}_{pos}_{ref}_{alt}"
                        
                        # 收集需要查询的位点
                        positions_to_query.add((chrom, pos, ref, alt))
                        
                        if variant_key not in all_variants:
                            all_variants[variant_key] = {
                                'chrom': chrom,
                                'pos': pos,
                                'id': var_id,
                                'ref': ref,
                                'alt': alt,
                                'qual': qual,
                                'filter': filt,
                                'info': info,
                                'format_samples': format_samples
                            }
                        else:
                            # 合并INFO
                            existing_info = all_variants[variant_key]['info']
                            existing_fields = {}
                            for part in existing_info.split(';'):
                                if '=' in part:
                                    k, v = part.split('=', 1)
                                    existing_fields[k] = v
                            
                            for part in info.split(';'):
                                if '=' in part:
                                    k, v = part.split('=', 1)
                                    if k not in existing_fields or existing_fields[k] == '.':
                                        existing_fields[k] = v
                            
                            new_info = ';'.join(f"{k}={v}" for k, v in existing_fields.items())
                            all_variants[variant_key]['info'] = new_info

    # 按染色体和位置排序
    def sort_key(item):
        variant_key, data = item
        chrom = data['chrom']
        pos = int(data['pos'])
        if chrom.startswith('chr'):
            chrom_num = chrom[3:]
            if chrom_num == 'X':
                chrom_num = 23
            elif chrom_num == 'Y':
                chrom_num = 24
            elif chrom_num == 'M':
                chrom_num = 25
            else:
                try:
                    chrom_num = int(chrom_num)
                except ValueError:
                    chrom_num = 26
        else:
            try:
                chrom_num = int(chrom)
            except ValueError:
                chrom_num = 26
        return (chrom_num, pos)

    sorted_variants = sorted(all_variants.items(), key=sort_key)

    # 写入合并后的文件
    with open(MERGED_VCF, 'w') as f:
        for line in header_lines:
            f.write(line)
        
        for variant_key, data in sorted_variants:
            line_parts = [
                data['chrom'], data['pos'], data['id'], data['ref'], data['alt'],
                data['qual'], data['filter'], data['info']
            ]
            if data['format_samples']:
                line_parts.append(data['format_samples'])
            f.write('\t'.join(line_parts) + '\n')

    logger.info(f"成功合并VCF文件: {len(sorted_variants)} 个唯一变异")
    return positions_to_query


def annotate_variants(positions_to_query):
    """对合并后的VCF文件进行完整注释（使用tabix按需查询）"""
    logger.info("开始对变异进行注释...")
    
    # 先加载dbscSNV（只查询需要的位点）
    dbscsnv_lookup = load_dbscsnv_for_positions(positions_to_query)
    
    dbnsfp_matched = 0
    dbscsnv_matched = 0
    total_count = 0
    
    with open(MERGED_VCF, 'r') as infile, open(ANNOTATED_VCF, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            total_count += 1
            fields = line.strip().split('\t')
            if len(fields) < 8:
                outfile.write(line)
                continue
            
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # 解析现有INFO字段
            info_dict = {}
            for part in info.split(';'):
                if '=' in part:
                    k, v = part.split('=', 1)
                    info_dict[k] = clean_info_value(k, v)
            
            # 使用tabix查询dbNSFP
            dbnsfp_data = query_dbnsfp_tabix(chrom, pos, ref, alt)
            if dbnsfp_data:
                dbnsfp_matched += 1
                for field_name, value in dbnsfp_data.items():
                    if field_name not in info_dict or info_dict[field_name] == '.':
                        info_dict[field_name] = clean_info_value(field_name, value)
            
            # 查询dbscSNV
            chrom_db = chrom[3:] if chrom.startswith('chr') else chrom
            dbscsnv_key = (chrom_db, pos, ref, alt)
            if dbscsnv_key in dbscsnv_lookup:
                dbscsnv_matched += 1
                dbscsnv_data = dbscsnv_lookup[dbscsnv_key]
                if dbscsnv_data.get('ada_score'):
                    if 'DBSCSNV_ADA' not in info_dict or info_dict['DBSCSNV_ADA'] == '.':
                        info_dict['DBSCSNV_ADA'] = dbscsnv_data['ada_score']
                if dbscsnv_data.get('rf_score'):
                    if 'DBSCSNV_RF' not in info_dict or info_dict['DBSCSNV_RF'] == '.':
                        info_dict['DBSCSNV_RF'] = dbscsnv_data['rf_score']
            
            # 重建INFO字段
            new_info_parts = []
            for field in FIELDS_TO_KEEP:
                if field in info_dict:
                    clean_val = clean_info_value(field, info_dict[field])
                    new_info_parts.append(f"{field}={clean_val}")
                else:
                    new_info_parts.append(f"{field}=.")
            
            new_info = ';'.join(new_info_parts)
            fields[7] = new_info
            
            outfile.write('\t'.join(fields) + '\n')
    
    logger.info(f"注释完成: 共处理 {total_count} 个变异")
    logger.info(f"  - dbNSFP匹配: {dbnsfp_matched} 个")
    logger.info(f"  - dbscSNV匹配: {dbscsnv_matched} 个")
    return True


def finalize_vcf():
    """最终清理并输出结果"""
    logger.info("生成最终结果文件...")
    
    # 直接复制，ANNOTATED_VCF已经是最终格式
    with open(ANNOTATED_VCF, 'r') as infile, open(FINAL_VCF, 'w') as outfile:
        for line in infile:
            outfile.write(line)
    
    # 统计最终结果
    variant_count = 0
    with open(FINAL_VCF, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                variant_count += 1
    
    logger.info(f"最终结果文件: {FINAL_VCF}")
    logger.info(f"最终变异数量: {variant_count}")
    return True


def main():
    """主函数"""
    import time
    start_time = time.time()
    
    logger.info("=" * 60)
    logger.info("步骤8: 整合三个决赛圈文件并补全注释 (优化版)")
    logger.info("=" * 60)
    
    # 检查输入文件
    if not check_input_files():
        logger.error("输入文件检查失败")
        sys.exit(1)
    
    # 检查tabix是否可用
    try:
        subprocess.run(['tabix', '--version'], capture_output=True, check=True)
    except Exception:
        logger.error("tabix未安装或不可用，请安装htslib: brew install htslib")
        sys.exit(1)
    
    # 合并VCF文件并收集需要查询的位点
    positions_to_query = merge_vcf_files()
    
    # 注释变异（使用按需查询）
    if not annotate_variants(positions_to_query):
        logger.error("注释失败")
        sys.exit(1)
    
    # 生成最终结果
    if not finalize_vcf():
        logger.error("最终处理失败")
        sys.exit(1)
    
    elapsed = time.time() - start_time
    logger.info("=" * 60)
    logger.info(f"步骤8完成! 总耗时: {elapsed:.1f} 秒")
    logger.info("=" * 60)


if __name__ == '__main__':
    main()
