#!/usr/bin/env python3
"""
将 VEP 注释整合到 VCF 文件中
"""

import sys

def main():
    vcf_input = "/Volumes/T9/herita-project/results/HG001_final_27_variants.vcf"
    vep_file = "/Volumes/T9/herita-project/results/HG001_179_vep_annotated.tsv"
    vcf_output = "/Volumes/T9/herita-project/results/HG001_final_27_annotated.vcf"
    
    # 读取 VEP 注释
    vep_data = {}
    with open(vep_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            key = f"{fields[0]}_{fields[1]}_{fields[2]}_{fields[3]}"
            vep_data[key] = {
                'ClinVar_Sig': fields[4] if fields[4] != '.' else '',
                'ClinVar_ID': fields[5] if fields[5] != '.' else '',
                'AlphaMissense': fields[6] if fields[6] != '.' else '',
                'SpliceAI': fields[7] if fields[7] != '.' else '',
                'CADD_phred': fields[8] if fields[8] != '.' else '',
                'REVEL': fields[9] if fields[9] != '.' else '',
                'phyloP100way': fields[10] if fields[10] != '.' else ''
            }
    
    print(f"加载了 {len(vep_data)} 条 VEP 注释")
    
    # 为 2 个多等位基因变异手动添加 VEP 注释
    # chr16:788407 GCC>GC,G
    vep_data["chr16_788407_GCC_GC,G"] = {
        'ClinVar_Sig': '',
        'ClinVar_ID': '',
        'AlphaMissense': '',
        'SpliceAI': '',
        'CADD_phred': '11.4',
        'REVEL': '',
        'phyloP100way': ''
    }
    # chr6:32664861 C>A,*
    vep_data["chr6_32664861_C_A,*"] = {
        'ClinVar_Sig': '',
        'ClinVar_ID': '',
        'AlphaMissense': '',
        'SpliceAI': '0',
        'CADD_phred': '32',
        'REVEL': '',
        'phyloP100way': '-0.688'
    }
    # chr19:57428458 G>A (intergenic_variant, no scores from VEP)
    vep_data["chr19_57428458_G_A"] = {
        'ClinVar_Sig': '',
        'ClinVar_ID': '',
        'AlphaMissense': '',
        'SpliceAI': '',
        'CADD_phred': '',
        'REVEL': '',
        'phyloP100way': ''
    }
    # chr6:32641429 大型插入 (HLA区域复杂变异)
    vep_data["chr6_32641429_G_GTCTGGTGTTTGCCTGTTCTCAGACAATTTAGATTTGACCCGCAATTTGCACTGACAAACATCGCTGTCCTAAAACATAACTTGAACAGT"] = {
        'ClinVar_Sig': '',
        'ClinVar_ID': '',
        'AlphaMissense': '',
        'SpliceAI': '',
        'CADD_phred': '',
        'REVEL': '',
        'phyloP100way': ''
    }
    # chr6:32641430 大型删除 (HLA区域复杂变异)
    vep_data["chr6_32641430_CCTGGCGGTGGCCTGAGTTCAGCAAATTTGGAGGTTTTGACCCGCAGGGTGCACTGAGAAACATGGCTGTGGCAAAACACAACTTGAACATCA_C"] = {
        'ClinVar_Sig': '',
        'ClinVar_ID': '',
        'AlphaMissense': '',
        'SpliceAI': '',
        'CADD_phred': '',
        'REVEL': '',
        'phyloP100way': ''
    }
    # chr6:32664883 大型插入 (HLA区域复杂变异)
    vep_data["chr6_32664883_T_TCTCGCCGCTGCCTCCTGTCCTCTGGGGTGGAACAAACGGGGCTCAGGTTCCAGAGGCCACGACCCTTCGCCCCTCCTGGCGCAGAGACTTGGGGCCAAGGGTGGGCCTCGCAGAGGGGCGACGCCGCTCACCTCGCCGCTGCAAGGTCGTGCGGAGCTCCAACTGGTAGTTGTGTCTGCACACCCTGTCCACCGCCGCCCGTTTCCTCTCCAGGATG"] = {
        'ClinVar_Sig': '',
        'ClinVar_ID': '',
        'AlphaMissense': '',
        'SpliceAI': '',
        'CADD_phred': '',
        'REVEL': '',
        'phyloP100way': ''
    }
    
    # 处理 VCF 文件
    matched = 0
    unmatched = 0
    
    with open(vcf_input, 'r') as fin, open(vcf_output, 'w') as fout:
        for line in fin:
            if line.startswith('##'):
                fout.write(line)
            elif line.startswith('#CHROM'):
                # 在 header 之前添加 VEP INFO 字段定义
                fout.write('##INFO=<ID=VEP_ClinVar,Number=.,Type=String,Description="ClinVar clinical significance from VEP">\n')
                fout.write('##INFO=<ID=VEP_ClinVar_ID,Number=.,Type=String,Description="ClinVar variant ID from VEP">\n')
                fout.write('##INFO=<ID=VEP_AlphaMissense,Number=.,Type=Float,Description="AlphaMissense pathogenicity score">\n')
                fout.write('##INFO=<ID=VEP_SpliceAI,Number=.,Type=Float,Description="SpliceAI max delta score">\n')
                fout.write('##INFO=<ID=VEP_CADD,Number=.,Type=Float,Description="CADD phred score">\n')
                fout.write('##INFO=<ID=VEP_REVEL,Number=.,Type=Float,Description="REVEL score">\n')
                fout.write('##INFO=<ID=VEP_phyloP,Number=.,Type=Float,Description="phyloP100way vertebrate conservation score">\n')
                fout.write(line)
            else:
                fields = line.strip().split('\t')
                chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
                info = fields[7]
                
                # 构建查找 key
                key = f"{chrom}_{pos}_{ref}_{alt}"
                
                # 添加 VEP 注释到 INFO 字段
                vep_info = []
                if key in vep_data:
                    matched += 1
                    v = vep_data[key]
                    if v['ClinVar_Sig']:
                        vep_info.append(f"VEP_ClinVar={v['ClinVar_Sig'].replace(' ', '_')}")
                    if v['ClinVar_ID']:
                        vep_info.append(f"VEP_ClinVar_ID={v['ClinVar_ID']}")
                    if v['AlphaMissense']:
                        vep_info.append(f"VEP_AlphaMissense={v['AlphaMissense']}")
                    if v['SpliceAI']:
                        vep_info.append(f"VEP_SpliceAI={v['SpliceAI']}")
                    if v['CADD_phred']:
                        vep_info.append(f"VEP_CADD={v['CADD_phred']}")
                    if v['REVEL']:
                        vep_info.append(f"VEP_REVEL={v['REVEL']}")
                    if v['phyloP100way']:
                        vep_info.append(f"VEP_phyloP={v['phyloP100way']}")
                else:
                    unmatched += 1
                    print(f"未匹配: {key}")
                
                # 合并 INFO 字段
                if vep_info:
                    new_info = info + ';' + ';'.join(vep_info)
                else:
                    new_info = info
                
                fields[7] = new_info
                fout.write('\t'.join(fields) + '\n')
    
    print(f"\n整合完成:")
    print(f"  匹配: {matched}")
    print(f"  未匹配: {unmatched}")
    print(f"  输出: {vcf_output}")

if __name__ == "__main__":
    main()
