#!/usr/bin/env python3
"""
VCF dbscSNV注释处理脚本 - 优化版
按需查询dbscSNV数据库，而不是加载整个文件

优化说明：
- 原版本加载整个dbscSNV(469MB, 12M+记录)到内存，耗时约7秒
- 优化版只查询vcf4中的变异位点（约271个）
- 预计耗时：2-3秒
"""

import os
import time

# 文件路径配置
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

DBSCSNV_FILE = os.path.join(PROJECT_ROOT, "dbscSNV1", "dbscSNV1.1_hg38_sorted.tsv")
VCF4_PATH = os.path.join(PROJECT_ROOT, "results", "vcf4.vcf")
FINAL_DBSCSNV_PATH = os.path.join(PROJECT_ROOT, "results", "final_dbscSNV.vcf")
VCF5_PATH = os.path.join(PROJECT_ROOT, "results", "vcf5.vcf")


def collect_variants_from_vcf4():
    """
    从vcf4收集所有需要查询的变异位点
    返回: set of (chrom_db, pos, ref, alt)
    """
    positions = set()
    with open(VCF4_PATH, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            # 转换为dbscSNV格式（不带chr前缀）
            chrom_db = chrom[3:] if chrom.startswith('chr') else chrom
            for single_alt in alt.split(','):
                positions.add((chrom_db, pos, ref, single_alt))
    
    return positions


def load_dbscsnv_for_positions(positions):
    """
    只加载需要的dbscSNV位点
    positions: set of (chrom, pos, ref, alt)
    返回: {(chrom, pos, ref, alt): (ada_score, rf_score)}
    """
    print(f"从dbscSNV查询 {len(positions)} 个位点...")
    
    lookup = {}
    positions_set = set(positions)  # 复制一份用于删除
    
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
                ada_score = fields[ada_idx] if fields[ada_idx] != '.' and fields[ada_idx] != '' else '.'
                rf_score = fields[rf_idx] if fields[rf_idx] != '.' and fields[rf_idx] != '' else '.'
                
                # 转换为数值
                try:
                    ada_val = float(ada_score) if ada_score != '.' else '.'
                except ValueError:
                    ada_val = '.'
                
                try:
                    rf_val = float(rf_score) if rf_score != '.' else '.'
                except ValueError:
                    rf_val = '.'
                
                lookup[key] = (ada_val, rf_val)
                
                # 移除已找到的位点
                positions_set.discard(key)
                
                # 如果所有位点都找到了，提前退出
                if not positions_set:
                    break
    
    print(f"dbscSNV匹配到 {len(lookup)} 个位点")
    return lookup


def should_enter_final_dbscSNV(ada_score, rf_score):
    """
    判断是否应该进入final_dbscSNV.vcf
    ada_score > 0.6 且 rf_score > 0.6
    """
    if ada_score == '.' or rf_score == '.':
        return False

    try:
        ada_val = float(ada_score)
        rf_val = float(rf_score)
        return ada_val > 0.6 and rf_val > 0.6
    except (ValueError, TypeError):
        return False


def process_vcf4():
    """处理vcf4文件（优化版 - 按需查询）"""
    start_time = time.time()
    
    print("🧬 dbscSNV注释 - 优化版 (按需查询)")
    print("=" * 60)
    
    # 检查输入文件
    if not os.path.exists(VCF4_PATH):
        print(f"❌ 输入文件不存在: {VCF4_PATH}")
        return
    
    # 收集需要查询的位点
    positions = collect_variants_from_vcf4()
    print(f"📊 vcf4变异数量: {len(positions)}")
    
    # 加载需要的dbscSNV数据
    dbscSNV_lookup = load_dbscsnv_for_positions(positions)
    
    print("\n开始处理VCF文件...")

    final_dbscSNV_count = 0
    vcf5_count = 0
    total_variants = 0
    matched_variants = 0
    high_score_variants = 0

    with open(VCF4_PATH, 'r') as infile, \
         open(FINAL_DBSCSNV_PATH, 'w') as final_out, \
         open(VCF5_PATH, 'w') as vcf5_out:

        # 写入VCF头部
        header_lines = []
        for line in infile:
            if line.startswith('#'):
                header_lines.append(line)
                final_out.write(line)
                vcf5_out.write(line)
            else:
                infile.seek(0)
                for _ in header_lines:
                    next(infile)
                break

        # 处理变异行
        for line in infile:
            if line.startswith('#') or line.strip() == '':
                continue

            total_variants += 1
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
                chrom_db = chrom[3:] if chrom.startswith('chr') else chrom
                key = (chrom_db, pos, ref, single_alt)

                # 从查找表获取注释
                ada_score = '.'
                rf_score = '.'
                if key in dbscSNV_lookup:
                    ada_score, rf_score = dbscSNV_lookup[key]
                    matched_variants += 1

                    if should_enter_final_dbscSNV(ada_score, rf_score):
                        high_score_variants += 1

                # 添加dbscSNV注释到INFO字段
                updated_info = info
                if updated_info == '.':
                    updated_info = f'DBSCSNV_ADA={ada_score};DBSCSNV_RF={rf_score}'
                else:
                    updated_info += f';DBSCSNV_ADA={ada_score};DBSCSNV_RF={rf_score}'

                # 更新字段
                updated_fields = fields.copy()
                updated_fields[4] = single_alt
                updated_fields[7] = updated_info
                updated_line = '\t'.join(updated_fields) + '\n'

                # 判断分类
                if should_enter_final_dbscSNV(ada_score, rf_score):
                    final_out.write(updated_line)
                    final_dbscSNV_count += 1
                    print(f"变异 {chrom}:{pos} {ref}>{single_alt} 进入dbscSNV决赛圈: ada={ada_score}, rf={rf_score}")
                else:
                    vcf5_out.write(updated_line)
                    vcf5_count += 1

    elapsed = time.time() - start_time
    
    print(f"\n处理完成！")
    print(f"总变异行数: {total_variants}")
    print(f"dbscSNV匹配成功: {matched_variants} ({matched_variants/total_variants*100:.1f}%)")
    print(f"高分变异(ada>0.6&rf>0.6): {high_score_variants}")
    print(f"进入dbscSNV决赛圈 (final_dbscSNV.vcf): {final_dbscSNV_count}")
    print(f"保留在vcf5中: {vcf5_count}")
    print(f"验证: {final_dbscSNV_count + vcf5_count} = {total_variants}")

    if matched_variants > 0:
        print(f"\n高分变异比例: {high_score_variants/matched_variants*100:.1f}% (相对于有注释的变异)")
    
    print(f"\n⏱️ 总耗时: {elapsed:.2f} 秒")


if __name__ == "__main__":
    process_vcf4()
