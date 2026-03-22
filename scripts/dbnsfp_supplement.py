#!/usr/bin/env python3
"""
dbNSFP Supplement - 为决赛圈中缺失dbNSFP注释的变异补充注释

检查每个变异是否有dbNSFP核心字段，如果没有，则从dbNSFP数据库中查找并补充
"""

import sys
import os
import gzip
import argparse
from collections import defaultdict

# dbNSFP 核心字段 - 检查多种可能的字段名
DBNSFP_CORE_FIELDS = [
    # 直接字段名
    'REVEL_score', 'CADD_phred', 'MetaRNN_score', 'AlphaMissense_score',
    'ClinPred_score', 'phyloP100way_vertebrate', 'genename',
    # VEP 前缀
    'VEP_REVEL', 'VEP_CADD', 'VEP_AlphaMissense', 'VEP_PhyloP',
    # dbNSFP 前缀
    'dbNSFP_REVEL', 'dbNSFP_MetaRNN', 'dbNSFP_AlphaMissense',
]

# dbNSFP 需要补充的所有字段 (根据 dbNSFP lite 文件的 header)
DBNSFP_FIELDS_TO_EXTRACT = [
    'aaref', 'aaalt', 'rs_dbSNP', 'aapos', 'genename',
    'clinvar_id', 'clinvar_clnsig', 'clinvar_trait', 'clinvar_review',
    'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id',
    'clinvar_OMIM_id', 'clinvar_Orphanet_id',
    'MetaRNN_score', 'REVEL_score', 'PrimateAI_score', 'ClinPred_score',
    'ESM1b_score', 'AlphaMissense_score', 'CADD_phred',
    'gnomAD4.1_joint_flag', 'gnomAD4.1_joint_POPMAX_AF',
]

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

def has_dbnsfp_annotation(info_dict):
    """检查是否有dbNSFP核心注释"""
    # 检查是否有至少3个核心字段
    count = 0
    for field in DBNSFP_CORE_FIELDS:
        val = info_dict.get(field, '.')
        if val not in ['.', '', 'NA', 'None']:
            count += 1
    return count >= 3

def load_dbnsfp_index(dbnsfp_file):
    """
    为dbNSFP文件创建位置索引
    返回: {chrom: {pos: offset}} 用于快速查找
    """
    # 由于dbNSFP文件太大，我们使用tabix进行查找
    # 这里我们先返回文件路径，在查询时使用tabix
    return dbnsfp_file

def query_dbnsfp(dbnsfp_file, chrom, pos, ref, alt):
    """使用tabix查询dbNSFP"""
    import subprocess
    
    # 标准化染色体名 - dbNSFP使用不带chr前缀的染色体名
    query_chrom = chrom.replace('chr', '')
    
    try:
        # 使用tabix查询
        cmd = ['tabix', dbnsfp_file, f'{query_chrom}:{pos}-{pos}']
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        
        if result.returncode != 0:
            return None
        
        # 解析结果
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            
            # dbNSFP格式: chr, pos, ref, alt, ...
            db_chr, db_pos, db_ref, db_alt = fields[0], fields[1], fields[2], fields[3]
            
            # 匹配变异 (ref和alt都需要匹配)
            if str(db_pos) == str(pos) and db_ref == ref and db_alt == alt:
                return fields
            
            # 也尝试匹配不同的ref (有些indel可能ref不完全匹配)
            if str(db_pos) == str(pos) and db_alt == alt:
                return fields
        
        return None
    except subprocess.TimeoutExpired:
        return None
    except Exception as e:
        return None

def get_dbnsfp_header(dbnsfp_file):
    """获取dbNSFP的列头"""
    import subprocess
    
    try:
        # 读取第一行（header）
        if dbnsfp_file.endswith('.gz'):
            # 使用 gzip 模块读取
            import gzip
            with gzip.open(dbnsfp_file, 'rt') as f:
                header_line = f.readline().strip()
        else:
            with open(dbnsfp_file, 'r') as f:
                header_line = f.readline().strip()
        
        # 去掉 # 开头
        if header_line.startswith('#'):
            header_line = header_line[1:]
        
        return header_line.split('\t')
    except Exception as e:
        print(f"读取dbNSFP header失败: {e}", file=sys.stderr)
        return None

# 需要取最大值的数值字段
NUMERIC_FIELDS = {
    'REVEL_score', 'CADD_phred', 'MetaRNN_score', 'AlphaMissense_score',
    'ClinPred_score', 'PrimateAI_score', 'ESM1b_score',
    'gnomAD4.1_joint_POPMAX_AF', 'aapos',
}

# 需要取第一个有效值的字符串字段
STRING_FIELDS = {
    'aaref', 'aaalt', 'rs_dbSNP', 'genename',
    'clinvar_id', 'clinvar_clnsig', 'clinvar_trait', 'clinvar_review',
    'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id',
    'clinvar_OMIM_id', 'clinvar_Orphanet_id', 'gnomAD4.1_joint_flag',
}

def normalize_multi_value(val, field_name):
    """
    规范化多值字段:
    1. 如果全是空值(.;.;.;.)，返回单个 '.'
    2. 对于数值字段，取最大值
    3. 对于字符串字段，取第一个有效值
    """
    if not val or val in ['.', '', 'NA', 'None']:
        return '.'
    
    # 分割多值 (支持分号和逗号)
    if ';' in val:
        parts = val.split(';')
    elif ',' in val:
        parts = val.split(',')
    else:
        return val  # 单值，直接返回
    
    # 过滤空值
    valid_parts = [p.strip() for p in parts if p.strip() and p.strip() not in ['.', '', 'NA', 'None', '-']]
    
    # 如果没有有效值，返回 '.'
    if not valid_parts:
        return '.'
    
    # 判断是否是数值字段
    if field_name in NUMERIC_FIELDS:
        # 尝试转换为数值并取最大值
        numeric_values = []
        for p in valid_parts:
            try:
                numeric_values.append(float(p))
            except ValueError:
                continue
        
        if numeric_values:
            max_val = max(numeric_values)
            # 格式化输出
            if max_val == int(max_val):
                return str(int(max_val))
            elif abs(max_val) < 0.0001:
                return f"{max_val:.6e}"
            else:
                return f"{max_val:.6f}".rstrip('0').rstrip('.')
        else:
            return '.'
    else:
        # 字符串字段，取第一个有效值
        return valid_parts[0]

def extract_dbnsfp_annotations(fields, header, wanted_fields=None):
    """从dbNSFP记录中提取注释"""
    if wanted_fields is None:
        wanted_fields = DBNSFP_FIELDS_TO_EXTRACT
    
    annotations = {}
    
    # 创建列名到索引的映射
    col_map = {name: i for i, name in enumerate(header)}
    
    for field in wanted_fields:
        if field in col_map and col_map[field] < len(fields):
            val = fields[col_map[field]]
            # 规范化多值字段
            normalized_val = normalize_multi_value(val, field)
            if normalized_val != '.':
                annotations[field] = normalized_val
    
    return annotations

def process_vcf(input_vcf, output_vcf, dbnsfp_file, verbose=False):
    """处理VCF文件，补充缺失的dbNSFP注释"""
    
    print("=" * 60, file=sys.stderr)
    print("dbNSFP Supplement - 补充缺失的dbNSFP注释", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"输入文件: {input_vcf}", file=sys.stderr)
    print(f"输出文件: {output_vcf}", file=sys.stderr)
    print(f"dbNSFP文件: {dbnsfp_file}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    # 获取dbNSFP header
    print("\n加载dbNSFP header...", file=sys.stderr)
    dbnsfp_header = get_dbnsfp_header(dbnsfp_file)
    if not dbnsfp_header:
        print("警告: 无法读取dbNSFP header，将跳过补充注释", file=sys.stderr)
        # 直接复制文件
        with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
            fout.write(fin.read())
        return {'total': 0, 'missing': 0, 'supplemented': 0, 'failed': 0}
    
    print(f"  dbNSFP共有 {len(dbnsfp_header)} 列", file=sys.stderr)
    
    # 统计
    stats = {'total': 0, 'has_annotation': 0, 'missing': 0, 'supplemented': 0, 'failed': 0}
    
    print("\n处理VCF文件...", file=sys.stderr)
    
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            
            stats['total'] += 1
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                fout.write(line)
                continue
            
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info_str = fields[7]
            
            info_dict = parse_info(info_str)
            
            # 检查是否需要补充注释
            if has_dbnsfp_annotation(info_dict):
                stats['has_annotation'] += 1
                fout.write(line)
                continue
            
            stats['missing'] += 1
            
            # 查询dbNSFP
            if verbose:
                print(f"  补充: {chrom}:{pos} {ref}>{alt}", file=sys.stderr)
            
            db_record = query_dbnsfp(dbnsfp_file, chrom, pos, ref, alt)
            
            if db_record:
                # 提取注释
                annotations = extract_dbnsfp_annotations(db_record, dbnsfp_header)
                
                if annotations:
                    stats['supplemented'] += 1
                    
                    # 添加到INFO
                    new_info_parts = [info_str]
                    for key, val in annotations.items():
                        # 只添加不存在的字段
                        if key not in info_dict:
                            new_info_parts.append(f"{key}={val}")
                    
                    fields[7] = ';'.join(new_info_parts)
                    fout.write('\t'.join(fields) + '\n')
                else:
                    stats['failed'] += 1
                    fout.write(line)
            else:
                stats['failed'] += 1
                fout.write(line)
    
    # 打印统计
    print("\n" + "=" * 60, file=sys.stderr)
    print("统计结果:", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"  总变异数: {stats['total']}", file=sys.stderr)
    print(f"  已有注释: {stats['has_annotation']}", file=sys.stderr)
    print(f"  缺失注释: {stats['missing']}", file=sys.stderr)
    print(f"  成功补充: {stats['supplemented']}", file=sys.stderr)
    print(f"  补充失败: {stats['failed']}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='dbNSFP Supplement - 补充缺失的dbNSFP注释')
    parser.add_argument('-i', '--input', required=True, help='输入VCF文件')
    parser.add_argument('-o', '--output', required=True, help='输出VCF文件')
    parser.add_argument('-d', '--dbnsfp', required=True, help='dbNSFP数据库文件')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细输出')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"错误: 输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.dbnsfp):
        print(f"错误: dbNSFP文件不存在: {args.dbnsfp}", file=sys.stderr)
        sys.exit(1)
    
    process_vcf(args.input, args.output, args.dbnsfp, args.verbose)
    print("\n完成!", file=sys.stderr)

if __name__ == '__main__':
    main()
