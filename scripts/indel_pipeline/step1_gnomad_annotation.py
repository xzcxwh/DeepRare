#!/usr/bin/env python3
"""
Indel Pipeline - Step 1: gnomAD 注释与频率筛选 (优化版 v8)
核心优化: 使用 bytes 模式读取，避免 text 解码开销
"""

import os
import sys
import subprocess
import re
from datetime import datetime
from collections import defaultdict

# 强制刷新所有输出
sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GNOMAD_DIR = os.path.join(PROJECT_DIR, "gnomAD_dataset")
WINDOW_SIZE = 100000  # 100kb

# 预编译正则表达式 (bytes 模式)
AF_PATTERN = re.compile(rb';AF=([0-9.eE+-]+)')

def log(msg):
    """强制刷新的打印函数"""
    print(msg, flush=True)


def get_gnomad_file(chrom):
    chrom_clean = chrom.replace('chr', '')
    gnomad_file = os.path.join(GNOMAD_DIR, f"gnomad.exomes.v4.1.sites.chr{chrom_clean}.vcf.bgz")
    if os.path.exists(gnomad_file):
        return gnomad_file
    return None


def query_gnomad_region_fast(gnomad_file, chrom, start, end):
    """
    高效查询 gnomAD 区域
    核心优化: 使用 bytes 模式，避免 text 解码开销
    单个 100kb 窗口查询时间: ~0.3秒
    """
    chrom_query = chrom if chrom.startswith('chr') else f"chr{chrom}"
    
    af_map = {}
    
    try:
        proc = subprocess.Popen(
            ['tabix', gnomad_file, f"{chrom_query}:{start}-{end}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL
        )
        
        for line in proc.stdout:
            # 使用 bytes.split 只分割前8个字段
            fields = line.split(b'\t', 8)
            if len(fields) < 8:
                continue
            
            pos = int(fields[1])
            ref = fields[3].decode('utf-8', errors='ignore')
            alt_field = fields[4].decode('utf-8', errors='ignore')
            info = fields[7][:500]  # 只取 INFO 字段前500字节
            
            # 快速提取 AF
            match = AF_PATTERN.search(b';' + info)
            if match:
                af_str = match.group(1).decode('utf-8', errors='ignore')
                afs = af_str.split(',')
                alts = alt_field.split(',')
                for i, alt in enumerate(alts):
                    if i < len(afs):
                        try:
                            af_map[(pos, ref, alt)] = float(afs[i])
                        except ValueError:
                            pass
        
        proc.wait()
                        
    except Exception as e:
        log(f"      ⚠️ 查询错误: {e}")
    
    return af_map


def annotate_and_filter_indel(input_vcf, output_vcf, max_af=0.01):
    log("=" * 60)
    log("Indel Pipeline - Step 1: gnomAD 注释与频率筛选 (v8)")
    log("=" * 60)
    log(f"输入文件: {input_vcf}")
    log(f"输出文件: {output_vcf}")
    log(f"筛选条件: gnomAD AF 有值且 < {max_af}")
    log(f"窗口大小: {WINDOW_SIZE//1000}kb")
    log("=" * 60)
    
    if not os.path.exists(input_vcf):
        log(f"❌ 错误: 输入文件不存在: {input_vcf}")
        return False
    
    start_time = datetime.now()
    
    # 读取变异
    log("\n[1/3] 读取 Indel 变异...")
    variants_by_chrom = defaultdict(list)
    header_lines = []
    total_variants = 0
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                if line.startswith('#CHROM'):
                    header_lines.insert(-1, '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomAD v4.1 exomes allele frequency">\n')
            else:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chrom = fields[0]
                    pos = int(fields[1])
                    variants_by_chrom[chrom].append({
                        'pos': pos,
                        'ref': fields[3],
                        'alt': fields[4],
                        'fields': fields
                    })
                    total_variants += 1
    
    for chrom in variants_by_chrom:
        variants_by_chrom[chrom].sort(key=lambda x: x['pos'])
    
    log(f"   读取 {total_variants:,} 个 Indel 变异")
    log(f"   分布在 {len(variants_by_chrom)} 个染色体上")
    
    os.makedirs(os.path.dirname(output_vcf), exist_ok=True)
    with open(output_vcf, 'w') as f:
        for line in header_lines:
            f.write(line)
    
    log("\n[2/3] 查询 gnomAD 并筛选...")
    
    has_af = 0
    passed_filter = 0
    no_af = 0
    high_af = 0
    total_windows = 0
    
    def chrom_sort_key(x):
        c = x.replace('chr', '')
        if c.isdigit():
            return (0, int(c))
        elif c == 'X':
            return (1, 23)
        elif c == 'Y':
            return (1, 24)
        else:
            return (2, ord(c[0]) if c else 0)
    
    for chrom in sorted(variants_by_chrom.keys(), key=chrom_sort_key):
        chrom_variants = variants_by_chrom[chrom]
        chrom_start_time = datetime.now()
        
        gnomad_file = get_gnomad_file(chrom)
        if not gnomad_file:
            log(f"   ⚠️ {chrom}: 无 gnomAD 数据，跳过 {len(chrom_variants):,} 个变异")
            no_af += len(chrom_variants)
            continue
        
        chrom_has_af = 0
        chrom_passed = 0
        chrom_no_af = 0
        passed_lines = []
        
        idx = 0
        chrom_windows = 0
        
        while idx < len(chrom_variants):
            window_start = chrom_variants[idx]['pos']
            window_end = window_start + WINDOW_SIZE
            
            window_variants = []
            while idx < len(chrom_variants) and chrom_variants[idx]['pos'] <= window_end:
                window_variants.append(chrom_variants[idx])
                idx += 1
            
            af_map = query_gnomad_region_fast(gnomad_file, chrom, window_start, window_end)
            chrom_windows += 1
            total_windows += 1
            
            for v in window_variants:
                key = (v['pos'], v['ref'], v['alt'])
                af = af_map.get(key)
                
                if af is not None:
                    chrom_has_af += 1
                    has_af += 1
                    
                    if af < max_af:
                        fields = v['fields'].copy()
                        info = fields[7] if len(fields) > 7 else "."
                        new_info = f"gnomAD_AF={af:.6e}" if info == "." else f"{info};gnomAD_AF={af:.6e}"
                        fields[7] = new_info
                        passed_lines.append('\t'.join(fields) + '\n')
                        passed_filter += 1
                        chrom_passed += 1
                    else:
                        high_af += 1
                else:
                    chrom_no_af += 1
                    no_af += 1
        
        with open(output_vcf, 'a') as f:
            for line in passed_lines:
                f.write(line)
        
        chrom_elapsed = (datetime.now() - chrom_start_time).total_seconds()
        log(f"   {chrom}: {len(chrom_variants):,} 变异, {chrom_windows} 窗口, {chrom_has_af:,} 有AF, {chrom_passed:,} 通过 ({chrom_elapsed:.1f}s)")
    
    elapsed = datetime.now() - start_time
    
    log("\n[3/3] 结果统计:")
    log("=" * 60)
    log(f"   总 Indel 变异数: {total_variants:,}")
    log(f"   有 gnomAD AF 的变异: {has_af:,} ({has_af*100/total_variants:.2f}%)")
    log(f"   无 gnomAD AF 的变异: {no_af:,} ({no_af*100/total_variants:.2f}%)")
    log(f"   AF >= {max_af} (过滤): {high_af:,}")
    log(f"   AF < {max_af} (通过): {passed_filter:,}")
    log("=" * 60)
    log(f"   ✅ 最终保留: {passed_filter:,} 个变异")
    log(f"   ⏱️ 处理时间: {elapsed}")
    
    return True


def main():
    input_vcf = os.path.join(PROJECT_DIR, "results", "indel", "indel_raw.vcf")
    output_vcf = os.path.join(PROJECT_DIR, "results", "indel", "indel_gnomad_filtered.vcf")
    
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_vcf = sys.argv[2]
    
    success = annotate_and_filter_indel(input_vcf, output_vcf, max_af=0.01)
    
    if success:
        log("\n✅ Step 1 完成!")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
