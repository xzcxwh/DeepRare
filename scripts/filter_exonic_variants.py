#!/usr/bin/env python3
"""
筛选外显子区域变异
使用 UCSC hg38 外显子坐标或 GENCODE 注释

用法:
    python filter_exonic_variants.py input.vcf output.vcf

依赖:
    - bcftools (已安装)
    - 外显子 BED 文件 (脚本会自动下载)
"""

import os
import sys
import subprocess
import gzip
from datetime import datetime

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANNOTATION_DIR = os.path.join(PROJECT_DIR, "annotation_data")


def download_exon_bed():
    """下载 GENCODE hg38 外显子坐标 BED 文件"""
    
    os.makedirs(ANNOTATION_DIR, exist_ok=True)
    
    exon_bed = os.path.join(ANNOTATION_DIR, "hg38_exons.bed.gz")
    exon_bed_tbi = exon_bed + ".tbi"
    
    if os.path.exists(exon_bed) and os.path.exists(exon_bed_tbi):
        print(f"   ✓ 外显子 BED 文件已存在: {exon_bed}")
        return exon_bed
    
    print("   下载 GENCODE v44 基因注释...")
    
    # GENCODE GTF 下载地址
    gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz"
    gtf_file = os.path.join(ANNOTATION_DIR, "gencode.v44.basic.annotation.gtf.gz")
    
    if not os.path.exists(gtf_file):
        print(f"   正在下载 GENCODE GTF (约 50MB)...")
        cmd = f"curl -L -o {gtf_file} {gencode_url}"
        result = subprocess.run(cmd, shell=True, capture_output=True)
        if result.returncode != 0:
            print(f"   ❌ 下载失败，尝试备用方法...")
            return create_exon_bed_from_ucsc()
    
    # 从 GTF 提取外显子区域
    print("   从 GTF 提取外显子区域...")
    temp_bed = os.path.join(ANNOTATION_DIR, "hg38_exons_temp.bed")
    
    with gzip.open(gtf_file, 'rt') as f_in, open(temp_bed, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'exon':
                continue
            
            chrom = fields[0]
            start = int(fields[3]) - 1  # BED is 0-based
            end = int(fields[4])
            
            # 只保留主要染色体
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            if chrom not in [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']:
                continue
            
            f_out.write(f"{chrom}\t{start}\t{end}\n")
    
    # 排序、合并、压缩
    print("   排序和合并重叠区域...")
    sorted_bed = os.path.join(ANNOTATION_DIR, "hg38_exons_sorted.bed")
    
    # 排序
    cmd_sort = f"sort -k1,1 -k2,2n {temp_bed} > {sorted_bed}"
    subprocess.run(cmd_sort, shell=True, check=True)
    
    # 合并重叠区域 (使用简单的 Python 实现，因为没有 bedtools)
    print("   合并重叠区域...")
    merged_bed = os.path.join(ANNOTATION_DIR, "hg38_exons_merged.bed")
    merge_overlapping_regions(sorted_bed, merged_bed)
    
    # 压缩和索引
    print("   压缩和建立索引...")
    cmd_bgzip = f"bgzip -c {merged_bed} > {exon_bed}"
    subprocess.run(cmd_bgzip, shell=True, check=True)
    
    cmd_tabix = f"tabix -p bed {exon_bed}"
    subprocess.run(cmd_tabix, shell=True, check=True)
    
    # 清理临时文件
    for f in [temp_bed, sorted_bed, merged_bed]:
        if os.path.exists(f):
            os.remove(f)
    
    print(f"   ✓ 外显子 BED 文件已创建: {exon_bed}")
    return exon_bed


def merge_overlapping_regions(input_bed, output_bed):
    """合并重叠的基因组区域"""
    
    regions = []
    current_chrom = None
    current_start = None
    current_end = None
    
    with open(input_bed, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            
            if current_chrom is None:
                current_chrom = chrom
                current_start = start
                current_end = end
            elif chrom == current_chrom and start <= current_end:
                # 重叠或相邻，扩展区域
                current_end = max(current_end, end)
            else:
                # 新区域，保存之前的
                regions.append((current_chrom, current_start, current_end))
                current_chrom = chrom
                current_start = start
                current_end = end
    
    # 保存最后一个区域
    if current_chrom is not None:
        regions.append((current_chrom, current_start, current_end))
    
    # 写入输出
    with open(output_bed, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")
    
    print(f"      合并后: {len(regions):,} 个外显子区域")


def create_exon_bed_from_ucsc():
    """备用方案：从 UCSC Table Browser 下载"""
    
    print("   尝试从 UCSC 下载外显子坐标...")
    
    exon_bed = os.path.join(ANNOTATION_DIR, "hg38_exons.bed.gz")
    
    # UCSC Table Browser API
    ucsc_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
    refgene_file = os.path.join(ANNOTATION_DIR, "refGene.txt.gz")
    
    cmd = f"curl -L -o {refgene_file} {ucsc_url}"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    
    if result.returncode != 0 or not os.path.exists(refgene_file):
        print("   ❌ 无法下载注释文件")
        return None
    
    # 解析 refGene 格式提取外显子
    temp_bed = os.path.join(ANNOTATION_DIR, "hg38_exons_temp.bed")
    
    with gzip.open(refgene_file, 'rt') as f_in, open(temp_bed, 'w') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            
            chrom = fields[2]
            exon_starts = fields[9].rstrip(',').split(',')
            exon_ends = fields[10].rstrip(',').split(',')
            
            # 只保留主要染色体
            if chrom not in [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']:
                continue
            
            for start, end in zip(exon_starts, exon_ends):
                if start and end:
                    f_out.write(f"{chrom}\t{start}\t{end}\n")
    
    # 排序、合并、压缩
    sorted_bed = os.path.join(ANNOTATION_DIR, "hg38_exons_sorted.bed")
    merged_bed = os.path.join(ANNOTATION_DIR, "hg38_exons_merged.bed")
    
    subprocess.run(f"sort -k1,1 -k2,2n {temp_bed} > {sorted_bed}", shell=True, check=True)
    merge_overlapping_regions(sorted_bed, merged_bed)
    
    subprocess.run(f"bgzip -c {merged_bed} > {exon_bed}", shell=True, check=True)
    subprocess.run(f"tabix -p bed {exon_bed}", shell=True, check=True)
    
    # 清理
    for f in [temp_bed, sorted_bed, merged_bed, refgene_file]:
        if os.path.exists(f):
            os.remove(f)
    
    return exon_bed


def filter_exonic_variants(input_vcf, output_vcf, exon_bed):
    """使用 Python 直接筛选外显子区域变异（避免排序要求）"""
    
    print(f"\n[2/2] 筛选外显子区域变异...")
    print(f"   输入: {input_vcf}")
    print(f"   外显子BED: {exon_bed}")
    
    # 加载外显子区域到内存（使用区间树结构）
    print("   加载外显子区域...")
    exon_regions = {}  # {chrom: [(start, end), ...]}
    
    with gzip.open(exon_bed, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            
            if chrom not in exon_regions:
                exon_regions[chrom] = []
            exon_regions[chrom].append((start, end))
    
    # 对每个染色体的区域排序
    for chrom in exon_regions:
        exon_regions[chrom].sort()
    
    total_regions = sum(len(regions) for regions in exon_regions.values())
    print(f"   ✓ 加载 {total_regions:,} 个外显子区域")
    
    # 使用二分查找检查位置是否在外显子区域内
    def is_in_exon(chrom, pos):
        if chrom not in exon_regions:
            return False
        
        regions = exon_regions[chrom]
        # 二分查找
        left, right = 0, len(regions) - 1
        
        while left <= right:
            mid = (left + right) // 2
            start, end = regions[mid]
            
            if start <= pos <= end:
                return True
            elif pos < start:
                right = mid - 1
            else:
                left = mid + 1
        
        return False
    
    # 筛选变异
    print("   筛选外显子区域变异...")
    input_count = 0
    output_count = 0
    
    with open(input_vcf, 'r') as f_in, open(output_vcf, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
                continue
            
            input_count += 1
            
            fields = line.split('\t', 3)
            chrom = fields[0]
            pos = int(fields[1])
            
            if is_in_exon(chrom, pos):
                f_out.write(line)
                output_count += 1
            
            # 进度显示
            if input_count % 500000 == 0:
                print(f"      已处理 {input_count:,} 变异, 保留 {output_count:,}")
    
    print(f"\n   结果统计:")
    print(f"   ├─ 输入变异数: {input_count:,}")
    print(f"   ├─ 外显子变异数: {output_count:,}")
    print(f"   └─ 保留比例: {output_count*100/input_count:.2f}%")
    
    return True


def main():
    print("=" * 60)
    print("筛选外显子区域变异")
    print("=" * 60)
    
    # 默认输入/输出
    input_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001.vcf.backup")
    output_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001_exonic.vcf")
    
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_vcf = sys.argv[2]
    
    if not os.path.exists(input_vcf):
        print(f"❌ 输入文件不存在: {input_vcf}")
        sys.exit(1)
    
    print(f"输入文件: {input_vcf}")
    print(f"输出文件: {output_vcf}")
    print("=" * 60)
    
    start_time = datetime.now()
    
    # Step 1: 准备外显子 BED 文件
    print("\n[1/2] 准备外显子坐标文件...")
    exon_bed = download_exon_bed()
    
    if not exon_bed:
        print("❌ 无法获取外显子坐标文件")
        sys.exit(1)
    
    # Step 2: 筛选
    success = filter_exonic_variants(input_vcf, output_vcf, exon_bed)
    
    elapsed = datetime.now() - start_time
    
    if success:
        print(f"\n✅ 完成! 耗时: {elapsed}")
        print(f"   输出文件: {output_vcf}")
    else:
        print(f"\n❌ 筛选失败")
        sys.exit(1)


if __name__ == "__main__":
    main()
