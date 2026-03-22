#!/usr/bin/env python3
"""
使用 bcftools 对 VCF 进行快速 gnomAD 注释和筛选
"""

import os
import sys
import subprocess
import glob

def filter_gnomad_bcftools(input_vcf, output_vcf):
    print("=" * 60)
    print("gnomAD 注释与筛选 (bcftools 高速版)")
    print("=" * 60)
    print(f"输入文件: {input_vcf}")
    print(f"输出文件: {output_vcf}")
    print("筛选条件: gnomAD_AF < 0.01 或为空")
    print("=" * 60)

    # 1. 准备工作
    # 压缩并索引输入文件 (这对性能至关重要，也防止 bcftools 报错)
    input_gz = input_vcf + ".tmp.gz"
    print("   压缩并索引输入文件...")
    # -f force overwrite
    subprocess.run(f"bcftools view {input_vcf} -O z -o {input_gz}", shell=True, check=True)
    subprocess.run(f"bcftools index -f {input_gz}", shell=True, check=True)

    # 获取所有染色体
    chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    
    temp_files = []
    
    print("\n[1/2] 按染色体并行注释...")
    for chrom in chroms:
        # 对应的 gnomAD 文件
        chrom_num = chrom.replace("chr", "")
        gnomad_file = f"/Volumes/T9/herita-project/gnomAD_dataset/gnomad.exomes.v4.1.sites.chr{chrom_num}.vcf.bgz"
        
        if not os.path.exists(gnomad_file):
            # print(f"   ⚠️ 跳过 {chrom} (无 gnomAD 数据)")
            continue
            
        # print(f"   处理 {chrom} ...")
        chrom_out = f"{output_vcf}.{chrom}.tmp.vcf.gz"
        
        # 构造 bcftools 命令
        # 关键修正: -c INFO/gnomAD_AF:=INFO/AF (将注释文件中的 AF 重命名为 gnomAD_AF)
        # 使用 -r chrom 限制处理区域
        cmd = (
            f"bcftools annotate -a {gnomad_file} "
            f"-c INFO/gnomAD_AF:=INFO/AF "
            f"-r {chrom} "
            f"{input_gz} | "
            f"bcftools filter -i 'gnomAD_AF<0.01 || gnomAD_AF=\".\"' "
            f"-O z -o {chrom_out}"
        )
        
        try:
            # capture_output=True to hide verbose output unless error
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
            temp_files.append(chrom_out)
            print(f"   ✅ {chrom} 完成")
        except subprocess.CalledProcessError as e:
            # 如果该染色体没有变异，bcftools 可能会报错或产生空文件
            # 检查 stderr 是否包含 "no matching records" 之类的
            err_msg = e.stderr.decode()
            if "index" in err_msg:
                 print(f"   ❌ {chrom} 索引错误: {err_msg}")
            else:
                 # 很多时候是因为该染色体在 input 中没有变异，bcftools -r 会返回空，导致后续 filter 报错
                 # print(f"   (跳过 {chrom}: 可能无变异)")
                 pass
    
    # 合并结果
    print("\n[2/2] 合并结果...")
    if temp_files:
        cmd_concat = f"bcftools concat -a {' '.join(temp_files)} -O v -o {output_vcf}"
        subprocess.run(cmd_concat, shell=True, check=True)
        print(f"✅ 完成! 输出文件: {output_vcf}")
        
        # 清理临时文件
        if os.path.exists(input_gz): os.remove(input_gz)
        if os.path.exists(input_gz + ".csi"): os.remove(input_gz + ".csi")
        for f in temp_files:
            if os.path.exists(f): os.remove(f)
    else:
        print("❌ 错误: 没有生成任何结果文件")

if __name__ == "__main__":
    input_vcf = "input_data/HG001_cds.vcf"
    output_vcf = "input_data/HG001_cds_filtered.vcf"
    
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_vcf = sys.argv[2]
        
    filter_gnomad_bcftools(input_vcf, output_vcf)
