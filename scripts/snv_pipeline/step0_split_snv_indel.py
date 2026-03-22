#!/usr/bin/env python3
"""
Step 0: 预处理 - 分离 SNV 和 Indel 变异
将输入 VCF 文件分成两个文件：
  - SNV (单核苷酸变异): REF 和 ALT 长度都为 1
  - Indel (插入/缺失): REF 和 ALT 长度不同
"""

import os
import sys
from datetime import datetime

# 项目路径
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def classify_variant(ref, alt):
    """
    根据 REF 和 ALT 判断变异类型
    返回: 'snv', 'insertion', 'deletion', 'complex', 'mnv'
    """
    # 处理多等位基因情况 (如 A,G,T)
    alts = alt.split(',')
    
    # 检查每个 ALT 等位基因
    types = set()
    for a in alts:
        ref_len = len(ref)
        alt_len = len(a)
        
        if ref_len == 1 and alt_len == 1:
            types.add('snv')
        elif ref_len == alt_len:
            types.add('mnv')  # 多核苷酸变异
        elif ref_len > alt_len:
            types.add('deletion')
        else:
            types.add('insertion')
    
    # 如果只有一种类型
    if len(types) == 1:
        return types.pop()
    
    # 混合类型（如一个等位基因是 SNV，另一个是 Indel）
    if 'snv' in types and ('insertion' in types or 'deletion' in types):
        return 'mixed'
    
    return 'complex'


def split_vcf(input_vcf, output_snv, output_indel):
    """
    将 VCF 文件分成 SNV 和 Indel 两个文件
    """
    print("=" * 60)
    print("Step 0: 分离 SNV 和 Indel 变异")
    print("=" * 60)
    print(f"输入文件: {input_vcf}")
    print(f"SNV 输出: {output_snv}")
    print(f"Indel 输出: {output_indel}")
    print("=" * 60)
    
    if not os.path.exists(input_vcf):
        print(f"❌ 错误: 输入文件不存在: {input_vcf}")
        return False
    
    start_time = datetime.now()
    
    # 统计计数
    stats = {
        'snv': 0,
        'insertion': 0,
        'deletion': 0,
        'mnv': 0,
        'complex': 0,
        'mixed': 0,
        'header_lines': 0
    }
    
    header_lines = []
    
    with open(input_vcf, 'r') as f_in, \
         open(output_snv, 'w') as f_snv, \
         open(output_indel, 'w') as f_indel:
        
        for line in f_in:
            # 头部行写入两个文件
            if line.startswith('#'):
                header_lines.append(line)
                f_snv.write(line)
                f_indel.write(line)
                stats['header_lines'] += 1
                continue
            
            # 解析变异行
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            
            ref = fields[3]
            alt = fields[4]
            
            var_type = classify_variant(ref, alt)
            stats[var_type] += 1
            
            # 分流写入
            if var_type == 'snv':
                f_snv.write(line)
            elif var_type in ['insertion', 'deletion', 'mnv', 'complex']:
                f_indel.write(line)
            elif var_type == 'mixed':
                # 混合类型同时写入两个文件
                f_snv.write(line)
                f_indel.write(line)
    
    elapsed = datetime.now() - start_time
    
    # 打印统计
    print("\n📊 变异类型统计:")
    print("-" * 40)
    
    total_snv = stats['snv']
    total_indel = stats['insertion'] + stats['deletion'] + stats['mnv'] + stats['complex']
    total = total_snv + total_indel + stats['mixed']
    
    print(f"SNV (单核苷酸变异):     {stats['snv']:>10,} ({stats['snv']*100/total:.2f}%)")
    print(f"Insertion (插入):       {stats['insertion']:>10,} ({stats['insertion']*100/total:.2f}%)")
    print(f"Deletion (缺失):        {stats['deletion']:>10,} ({stats['deletion']*100/total:.2f}%)")
    print(f"MNV (多核苷酸变异):     {stats['mnv']:>10,} ({stats['mnv']*100/total:.2f}%)")
    print(f"Complex:                {stats['complex']:>10,} ({stats['complex']*100/total:.2f}%)")
    print(f"Mixed (混合):           {stats['mixed']:>10,} ({stats['mixed']*100/total:.2f}%)")
    print("-" * 40)
    print(f"Total SNV:              {total_snv:>10,} ({total_snv*100/total:.2f}%)")
    print(f"Total Indel:            {total_indel:>10,} ({total_indel*100/total:.2f}%)")
    print(f"总计:                   {total:>10,}")
    
    print(f"\n⏱️ 处理时间: {elapsed}")
    
    # 验证输出文件
    snv_count = sum(1 for line in open(output_snv) if not line.startswith('#'))
    indel_count = sum(1 for line in open(output_indel) if not line.startswith('#'))
    
    print(f"\n✅ 输出文件验证:")
    print(f"   SNV 文件变异数: {snv_count:,}")
    print(f"   Indel 文件变异数: {indel_count:,}")
    
    return True


def main():
    # 默认路径
    input_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001.vcf")
    output_snv = os.path.join(PROJECT_DIR, "input_data", "HG001_snv.vcf")
    output_indel = os.path.join(PROJECT_DIR, "input_data", "HG001_indel.vcf")
    
    # 命令行参数
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_snv = sys.argv[2]
    if len(sys.argv) > 3:
        output_indel = sys.argv[3]
    
    success = split_vcf(input_vcf, output_snv, output_indel)
    
    if success:
        print("\n✅ Step 0 完成!")
        print("\n📌 下一步:")
        print(f"   - SNV 流水线: python scripts/run_pipeline.py {output_snv}")
        print(f"   - Indel 流水线: 待实现")
    else:
        print("\n❌ Step 0 失败!")
        sys.exit(1)


if __name__ == "__main__":
    main()
