#!/usr/bin/env python3
"""
Indel Pipeline - Step 0: 从 VCF 文件中提取 Indel 变异
Indel 定义: REF 和 ALT 的长度不同
"""

import os
import sys
from datetime import datetime

# 项目路径
PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def is_indel(ref, alt):
    """
    判断是否是 Indel 变异
    Indel: REF 和 ALT 长度不同
    """
    # 处理多等位基因情况
    alts = alt.split(',')
    
    for a in alts:
        if len(ref) != len(a):
            return True
    return False


def extract_indel(input_vcf, output_indel):
    """
    从 VCF 文件中提取 Indel 变异
    """
    print("=" * 60)
    print("Indel Pipeline - Step 0: 提取 Indel 变异")
    print("=" * 60)
    print(f"输入文件: {input_vcf}")
    print(f"输出文件: {output_indel}")
    print("=" * 60)
    
    if not os.path.exists(input_vcf):
        print(f"❌ 错误: 输入文件不存在: {input_vcf}")
        return False
    
    start_time = datetime.now()
    
    snv_count = 0
    indel_count = 0
    insertion_count = 0
    deletion_count = 0
    
    with open(input_vcf, 'r') as f_in, open(output_indel, 'w') as f_out:
        for line in f_in:
            # 头部行直接写入
            if line.startswith('#'):
                f_out.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            
            ref = fields[3]
            alt = fields[4]
            
            if is_indel(ref, alt):
                f_out.write(line)
                indel_count += 1
                
                # 统计插入/缺失
                if len(ref) < len(alt.split(',')[0]):
                    insertion_count += 1
                else:
                    deletion_count += 1
            else:
                snv_count += 1
    
    elapsed = datetime.now() - start_time
    
    print(f"\n📊 变异统计:")
    print(f"   SNV (已跳过): {snv_count:,}")
    print(f"   Indel (已提取): {indel_count:,}")
    print(f"     - Insertion: {insertion_count:,}")
    print(f"     - Deletion: {deletion_count:,}")
    print(f"\n⏱️ 处理时间: {elapsed}")
    print(f"\n✅ Indel 文件已生成: {output_indel}")
    
    return True


def main():
    input_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001.vcf")
    output_indel = os.path.join(PROJECT_DIR, "results", "indel", "indel_raw.vcf")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(output_indel), exist_ok=True)
    
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_indel = sys.argv[2]
    
    success = extract_indel(input_vcf, output_indel)
    
    if success:
        print("\n✅ Step 0 完成!")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
