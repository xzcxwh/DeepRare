#!/usr/bin/env python3
"""
快速 gnomAD AF 注释与筛选
使用预构建的位置索引进行快速查找
筛选条件: gnomAD_AF < 0.01 或 AF 为空 → 保留
"""

import os
import sys
import subprocess
import gzip
from datetime import datetime
from collections import defaultdict

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GNOMAD_DIR = os.path.join(PROJECT_DIR, "gnomAD_dataset")


def get_gnomad_file(chrom):
    """获取染色体对应的 gnomAD 文件"""
    chrom_clean = chrom.replace('chr', '')
    gnomad_file = os.path.join(GNOMAD_DIR, f"gnomad.exomes.v4.1.sites.chr{chrom_clean}.vcf.bgz")
    if os.path.exists(gnomad_file):
        return gnomad_file
    return None


def build_gnomad_af_index(gnomad_file, positions_set, chrom):
    """
    为指定位置集合构建 gnomAD AF 索引
    positions_set: set of (pos, ref, alt)
    返回: {(pos, ref, alt): af}
    """
    if not positions_set:
        return {}
    
    # 获取位置范围
    all_pos = [p[0] for p in positions_set]
    min_pos = min(all_pos)
    max_pos = max(all_pos)
    
    chrom_query = chrom if chrom.startswith('chr') else f"chr{chrom}"
    
    af_map = {}
    
    # 使用 tabix 查询整个范围，用 awk 只提取需要的字段
    cmd = f"tabix {gnomad_file} {chrom_query}:{min_pos}-{max_pos} | cut -f2,4,5,8"
    
    try:
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        
        for line in proc.stdout:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            pos = int(parts[0])
            ref = parts[1]
            alts = parts[2].split(',')
            info = parts[3]
            
            # 快速提取 AF
            af_start = info.find('AF=')
            if af_start == -1:
                continue
            
            af_end = info.find(';', af_start)
            if af_end == -1:
                af_str = info[af_start+3:]
            else:
                af_str = info[af_start+3:af_end]
            
            afs = af_str.split(',')
            
            for i, alt in enumerate(alts):
                key = (pos, ref, alt)
                if key in positions_set and i < len(afs):
                    try:
                        af_map[key] = float(afs[i])
                    except ValueError:
                        pass
        
        proc.wait()
        
    except Exception as e:
        print(f"   ⚠️ 查询错误: {e}", flush=True)
    
    return af_map


def annotate_and_filter(input_vcf, output_vcf, max_af=0.01):
    """注释并筛选变异"""
    
    print("=" * 60, flush=True)
    print("gnomAD AF 注释与筛选 (优化版)", flush=True)
    print("=" * 60, flush=True)
    print(f"输入: {input_vcf}", flush=True)
    print(f"输出: {output_vcf}", flush=True)
    print(f"筛选: 保留 AF < {max_af} 或 AF 为空", flush=True)
    print("=" * 60, flush=True)
    
    start_time = datetime.now()
    
    # 读取所有变异，按染色体分组
    print("\n[1/3] 读取变异并按染色体分组...", flush=True)
    
    variants_by_chrom = defaultdict(list)
    header_lines = []
    total = 0
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                if line.startswith('#CHROM'):
                    # 插入 gnomAD_AF 定义
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
                        'fields': fields
                    })
                    total += 1
    
    print(f"   读取 {total:,} 个变异，分布在 {len(variants_by_chrom)} 个染色体", flush=True)
    
    # 写入头部
    os.makedirs(os.path.dirname(output_vcf) if os.path.dirname(output_vcf) else '.', exist_ok=True)
    with open(output_vcf, 'w') as f:
        for line in header_lines:
            f.write(line)
    
    # 处理每个染色体
    print("\n[2/3] 查询 gnomAD 并筛选...", flush=True)
    
    stats = {'has_af': 0, 'no_af': 0, 'passed': 0, 'filtered': 0}
    
    # 排序染色体
    def chrom_key(c):
        c = c.replace('chr', '')
        if c.isdigit():
            return (0, int(c))
        elif c == 'X':
            return (1, 23)
        elif c == 'Y':
            return (1, 24)
        else:
            return (2, ord(c[0]) if c else 0)
    
    for chrom in sorted(variants_by_chrom.keys(), key=chrom_key):
        chrom_variants = variants_by_chrom[chrom]
        chrom_start = datetime.now()
        
        gnomad_file = get_gnomad_file(chrom)
        
        if not gnomad_file:
            # 无 gnomAD 数据，保留所有变异（AF 为空）
            print(f"   {chrom}: {len(chrom_variants):,} 变异 (无gnomAD数据，全部保留)", flush=True)
            with open(output_vcf, 'a') as f:
                for v in chrom_variants:
                    f.write('\t'.join(v['fields']) + '\n')
            stats['no_af'] += len(chrom_variants)
            stats['passed'] += len(chrom_variants)
            continue
        
        # 构建查找集合
        positions_set = set((v['pos'], v['ref'], v['alt']) for v in chrom_variants)
        
        # 批量查询 gnomAD
        af_map = build_gnomad_af_index(gnomad_file, positions_set, chrom)
        
        # 筛选并写入
        chrom_passed = 0
        chrom_has_af = 0
        passed_lines = []
        
        for v in chrom_variants:
            key = (v['pos'], v['ref'], v['alt'])
            af = af_map.get(key)
            
            if af is not None:
                chrom_has_af += 1
                stats['has_af'] += 1
                
                if af < max_af:
                    # 添加 gnomAD_AF 注释
                    fields = v['fields'].copy()
                    info = fields[7] if len(fields) > 7 else "."
                    new_info = f"gnomAD_AF={af:.6e}" if info == "." else f"{info};gnomAD_AF={af:.6e}"
                    fields[7] = new_info
                    passed_lines.append('\t'.join(fields) + '\n')
                    chrom_passed += 1
                    stats['passed'] += 1
                else:
                    stats['filtered'] += 1
            else:
                # AF 为空，保留
                stats['no_af'] += 1
                stats['passed'] += 1
                passed_lines.append('\t'.join(v['fields']) + '\n')
                chrom_passed += 1
        
        # 写入结果
        with open(output_vcf, 'a') as f:
            for line in passed_lines:
                f.write(line)
        
        elapsed = (datetime.now() - chrom_start).total_seconds()
        print(f"   {chrom}: {len(chrom_variants):,} → {chrom_passed:,} 保留 ({chrom_has_af} 有AF) [{elapsed:.1f}s]", flush=True)
    
    # 统计
    elapsed = datetime.now() - start_time
    
    print("\n[3/3] 结果统计:", flush=True)
    print("=" * 60, flush=True)
    print(f"   输入变异: {total:,}", flush=True)
    print(f"   有 gnomAD AF: {stats['has_af']:,} ({stats['has_af']*100/total:.1f}%)", flush=True)
    print(f"   无 gnomAD AF: {stats['no_af']:,} ({stats['no_af']*100/total:.1f}%)", flush=True)
    print(f"   过滤 (AF >= {max_af}): {stats['filtered']:,}", flush=True)
    print(f"   保留: {stats['passed']:,} ({stats['passed']*100/total:.1f}%)", flush=True)
    print("=" * 60, flush=True)
    print(f"✅ 完成! 耗时: {elapsed}", flush=True)
    print(f"   输出: {output_vcf}", flush=True)
    
    return True


def main():
    input_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001_exonic.vcf")
    output_vcf = os.path.join(PROJECT_DIR, "input_data", "HG001_exonic_rare.vcf")
    
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_vcf = sys.argv[2]
    
    if not os.path.exists(input_vcf):
        print(f"❌ 文件不存在: {input_vcf}")
        sys.exit(1)
    
    annotate_and_filter(input_vcf, output_vcf, max_af=0.01)


if __name__ == "__main__":
    main()
