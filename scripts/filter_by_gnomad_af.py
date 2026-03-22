#!/usr/bin/env python3
"""
使用 gnomAD AF 筛选变异
- 如果 gnomAD_AF < 0.01 或为空(无AF): 保留
- 如果 gnomAD_AF >= 0.01: 去除

用法:
    python filter_by_gnomad_af.py input.vcf output.vcf
"""

import os
import sys
import subprocess
import re
from datetime import datetime
from collections import defaultdict

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GNOMAD_DIR = os.path.join(PROJECT_DIR, "gnomAD_dataset")

# 确保输出立即刷新
sys.stdout.reconfigure(line_buffering=True)


def log(msg):
    print(msg, flush=True)


def get_gnomad_file(chrom):
    """获取对应染色体的 gnomAD 文件"""
    chrom_clean = chrom.replace('chr', '')
    gnomad_file = os.path.join(GNOMAD_DIR, f"gnomad.exomes.v4.1.sites.chr{chrom_clean}.vcf.bgz")
    if os.path.exists(gnomad_file):
        return gnomad_file
    return None


def query_gnomad_af(gnomad_file, chrom, positions):
    """
    批量查询 gnomAD AF
    positions: set of (pos, ref, alt)
    返回: {(pos, ref, alt): af, ...}
    """
    if not positions:
        return {}
    
    # 找出位置范围
    pos_list = [p[0] for p in positions]
    min_pos = min(pos_list)
    max_pos = max(pos_list)
    
    chrom_query = chrom if chrom.startswith('chr') else f"chr{chrom}"
    
    af_map = {}
    af_pattern = re.compile(r';AF=([^;]+)')
    
    try:
        # 使用 tabix 查询区间
        proc = subprocess.Popen(
            ['tabix', gnomad_file, f"{chrom_query}:{min_pos}-{max_pos}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True
        )
        
        for line in proc.stdout:
            fields = line.split('\t', 9)
            if len(fields) < 8:
                continue
            
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(',')
            info = fields[7]
            
            # 快速提取 AF
            match = af_pattern.search(';' + info)
            if match:
                afs = match.group(1).split(',')
                for i, alt in enumerate(alts):
                    key = (pos, ref, alt)
                    if key in positions and i < len(afs):
                        try:
                            af_map[key] = float(afs[i])
                        except ValueError:
                            pass
        
        proc.wait()
        
    except Exception as e:
        log(f"   ⚠️ 查询错误: {e}")
    
    return af_map


def filter_by_gnomad_af(input_vcf, output_vcf, max_af=0.01):
    """
    使用 gnomAD AF 筛选变异
    保留: AF < max_af 或无 AF
    去除: AF >= max_af
    """
    log("=" * 60)
    log("gnomAD AF 筛选")
    log("=" * 60)
    log(f"输入文件: {input_vcf}")
    log(f"输出文件: {output_vcf}")
    log(f"筛选条件: 保留 gnomAD_AF < {max_af} 或无AF")
    log("=" * 60)
    
    if not os.path.exists(input_vcf):
        log(f"❌ 错误: 输入文件不存在: {input_vcf}")
        return False
    
    start_time = datetime.now()
    
    # 读取变异并按染色体分组
    log("\n[1/3] 读取变异...")
    variants_by_chrom = defaultdict(list)
    header_lines = []
    total_variants = 0
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                # 添加 gnomAD_AF INFO 字段定义
                if line.startswith('#CHROM'):
                    # 检查是否已有 gnomAD_AF 定义
                    has_gnomad_af = any('gnomAD_AF' in h for h in header_lines)
                    if not has_gnomad_af:
                        header_lines.insert(-1, '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomAD v4.1 exomes allele frequency">\n')
            else:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    variants_by_chrom[chrom].append({
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'fields': fields,
                        'line': line
                    })
                    total_variants += 1
    
    log(f"   读取 {total_variants:,} 个变异")
    log(f"   分布在 {len(variants_by_chrom)} 个染色体上")
    
    # 准备输出
    os.makedirs(os.path.dirname(output_vcf) if os.path.dirname(output_vcf) else '.', exist_ok=True)
    
    with open(output_vcf, 'w') as f:
        for line in header_lines:
            f.write(line)
    
    log("\n[2/3] 查询 gnomAD 并筛选...")
    
    # 统计
    has_af = 0
    no_af = 0
    passed_low_af = 0
    passed_no_af = 0
    filtered_high_af = 0
    
    # 染色体排序
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
        chrom_start = datetime.now()
        
        # 构建查询集合
        positions = set()
        for v in chrom_variants:
            positions.add((v['pos'], v['ref'], v['alt']))
        
        # 查询 gnomAD
        gnomad_file = get_gnomad_file(chrom)
        af_map = {}
        
        if gnomad_file:
            af_map = query_gnomad_af(gnomad_file, chrom, positions)
        
        # 处理变异
        passed_lines = []
        chrom_has_af = 0
        chrom_no_af = 0
        chrom_passed = 0
        chrom_filtered = 0
        
        for v in chrom_variants:
            key = (v['pos'], v['ref'], v['alt'])
            af = af_map.get(key)
            
            if af is not None:
                chrom_has_af += 1
                has_af += 1
                
                if af < max_af:
                    # AF < 0.01, 保留并添加注释
                    fields = v['fields'].copy()
                    info = fields[7] if len(fields) > 7 else "."
                    new_info = f"gnomAD_AF={af:.6e}" if info == "." else f"{info};gnomAD_AF={af:.6e}"
                    fields[7] = new_info
                    passed_lines.append('\t'.join(fields) + '\n')
                    passed_low_af += 1
                    chrom_passed += 1
                else:
                    # AF >= 0.01, 过滤掉
                    filtered_high_af += 1
                    chrom_filtered += 1
            else:
                # 无 AF, 保留原样
                chrom_no_af += 1
                no_af += 1
                passed_lines.append(v['line'])
                passed_no_af += 1
                chrom_passed += 1
        
        # 写入结果
        with open(output_vcf, 'a') as f:
            for line in passed_lines:
                f.write(line)
        
        chrom_elapsed = (datetime.now() - chrom_start).total_seconds()
        log(f"   {chrom}: {len(chrom_variants):,} 变异, 有AF {chrom_has_af:,}, 无AF {chrom_no_af:,}, 通过 {chrom_passed:,}, 过滤 {chrom_filtered:,} ({chrom_elapsed:.1f}s)")
    
    elapsed = datetime.now() - start_time
    
    total_passed = passed_low_af + passed_no_af
    
    log("\n[3/3] 结果统计:")
    log("=" * 60)
    log(f"   输入变异数: {total_variants:,}")
    log(f"   有 gnomAD AF: {has_af:,} ({has_af*100/total_variants:.2f}%)")
    log(f"   无 gnomAD AF: {no_af:,} ({no_af*100/total_variants:.2f}%)")
    log("-" * 60)
    log(f"   AF < {max_af} (保留): {passed_low_af:,}")
    log(f"   无 AF (保留): {passed_no_af:,}")
    log(f"   AF >= {max_af} (过滤): {filtered_high_af:,}")
    log("=" * 60)
    log(f"   ✅ 最终保留: {total_passed:,} 个变异 ({total_passed*100/total_variants:.2f}%)")
    log(f"   ⏱️ 处理时间: {elapsed}")
    
    return True


def main():
    # 默认输入/输出
    input_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001_exonic.vcf")
    output_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001_exonic_rare.vcf")
    
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_vcf = sys.argv[2]
    
    success = filter_by_gnomad_af(input_vcf, output_vcf, max_af=0.01)
    
    if success:
        log("\n✅ 筛选完成!")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
