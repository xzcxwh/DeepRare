#!/usr/bin/env python3
"""
VEP REST API Annotation Script
对 VCF 文件中的变异进行 VEP 注释，获取以下信息：
- ClinVar (通过 colocated_variants)
- AlphaMissense
- SpliceAI delta score
- CADD phred
- REVEL
- 保守性 (通过 dbNSFP 的 phyloP100way)
"""

import requests
import json
import time
import sys
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

# VEP REST API 配置
VEP_API_URL = "https://rest.ensembl.org/vep/human/region"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

# API 参数 - 根据 VEP REST API 文档配置
VEP_PARAMS = {
    "AlphaMissense": 1,           # AlphaMissense 病理性预测
    "CADD": 1,                    # CADD 分数 (snv_indels)
    "REVEL": 1,                   # REVEL 分数
    "SpliceAI": 1,                # SpliceAI 剪接预测
    "dbNSFP": "phyloP100way_vertebrate",  # dbNSFP 保守性分数
    "hgvs": 1,
    "canonical": 1,
    "mane": 1,
    "numbers": 1,
    "variant_class": 1
}

def parse_vcf(vcf_file):
    """解析 VCF 文件，提取变异信息"""
    variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            
            chrom = fields[0].replace('chr', '')  # VEP 使用不带 chr 前缀的染色体
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(',')[0]  # 只取第一个 ALT
            
            # 构建 VEP 格式的变异字符串
            # 格式: chrom start end allele_string strand
            if len(ref) == 1 and len(alt) == 1:
                # SNV
                vep_string = f"{chrom} {pos} {pos} {ref}/{alt} 1"
            elif len(ref) > len(alt):
                # Deletion
                deleted = ref[len(alt):]
                start = pos + len(alt)
                end = pos + len(ref) - 1
                vep_string = f"{chrom} {start} {end} {deleted}/- 1"
            else:
                # Insertion
                inserted = alt[len(ref):]
                start = pos + len(ref)
                end = pos + len(ref) - 1
                vep_string = f"{chrom} {start} {end} -/{inserted} 1"
            
            variants.append({
                'original': f"{fields[0]}\t{fields[1]}\t{fields[3]}\t{alt}",
                'vep_string': vep_string,
                'chrom': fields[0],
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'line': line.strip()
            })
    
    return variants

def call_vep_api(variants_batch):
    """调用 VEP REST API"""
    payload = {
        "variants": [v['vep_string'] for v in variants_batch]
    }
    payload.update(VEP_PARAMS)
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            response = requests.post(
                VEP_API_URL,
                headers=HEADERS,
                json=payload,
                timeout=120
            )
            
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 429:
                # Rate limited
                wait_time = int(response.headers.get('Retry-After', 5))
                print(f"  Rate limited, waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"  API error: {response.status_code} - {response.text[:200]}")
                time.sleep(2)
        except Exception as e:
            print(f"  Request error: {e}")
            time.sleep(2)
    
    return None

def extract_annotations(vep_result):
    """从 VEP 结果中提取注释信息"""
    annotations = {
        'clinvar_clnsig': '.',
        'clinvar_trait': '.',
        'alpha_missense': '.',
        'spliceai_max': '.',
        'cadd_phred': '.',
        'revel': '.',
        'phylop100way': '.'
    }
    
    if not vep_result:
        return annotations
    
    # 从 colocated_variants 获取 ClinVar
    if 'colocated_variants' in vep_result:
        for cv in vep_result['colocated_variants']:
            if 'clin_sig' in cv and annotations['clinvar_clnsig'] == '.':
                annotations['clinvar_clnsig'] = ','.join(cv['clin_sig']) if isinstance(cv['clin_sig'], list) else str(cv['clin_sig'])
            if 'var_synonyms' in cv and 'ClinVar' in cv['var_synonyms'] and annotations['clinvar_trait'] == '.':
                clinvar_ids = cv['var_synonyms']['ClinVar']
                annotations['clinvar_trait'] = clinvar_ids[0] if clinvar_ids else '.'
    
    # 从 transcript_consequences 获取其他注释
    if 'transcript_consequences' in vep_result:
        for tc in vep_result['transcript_consequences']:
            # AlphaMissense - 可能是嵌套对象 {'am_class': 'xxx', 'am_pathogenicity': 0.xxx}
            if 'alphamissense' in tc and annotations['alpha_missense'] == '.':
                am = tc['alphamissense']
                if isinstance(am, dict) and 'am_pathogenicity' in am:
                    annotations['alpha_missense'] = str(am['am_pathogenicity'])
                elif am is not None:
                    annotations['alpha_missense'] = str(am)
            # 也检查直接的 am_pathogenicity 字段
            if 'am_pathogenicity' in tc and annotations['alpha_missense'] == '.':
                annotations['alpha_missense'] = str(tc['am_pathogenicity'])
            
            # CADD phred
            if 'cadd_phred' in tc and annotations['cadd_phred'] == '.':
                val = tc['cadd_phred']
                if val is not None and str(val) != 'invalid_field':
                    annotations['cadd_phred'] = str(val)
            
            # REVEL - 优先使用简单的 revel 字段
            if 'revel' in tc and annotations['revel'] == '.':
                val = tc['revel']
                if val is not None and str(val) != '.':
                    annotations['revel'] = str(val)
            
            # SpliceAI - 字段名是大写 DS_AG 等
            if 'spliceai' in tc and annotations['spliceai_max'] == '.':
                spliceai = tc['spliceai']
                if isinstance(spliceai, dict):
                    scores = []
                    for key in ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']:
                        if key in spliceai and spliceai[key] is not None:
                            try:
                                scores.append(float(spliceai[key]))
                            except (ValueError, TypeError):
                                pass
                    if scores:
                        annotations['spliceai_max'] = f"{max(scores):.4f}"
            
            # phyloP100way
            if 'phylop100way_vertebrate' in tc and annotations['phylop100way'] == '.':
                val = tc['phylop100way_vertebrate']
                if val is not None:
                    annotations['phylop100way'] = str(val)
    
    # 从顶层获取 (有些注释在顶层)
    if 'cadd_phred' in vep_result and annotations['cadd_phred'] == '.':
        val = vep_result['cadd_phred']
        if val is not None and str(val) != 'invalid_field':
            annotations['cadd_phred'] = str(val)
    
    return annotations

def main():
    vcf_file = "/Volumes/T9/herita-project/results/HG001_543_high_moderate.vcf"
    output_file = "/Volumes/T9/herita-project/results/HG001_179_vep_annotated.tsv"
    
    print("=" * 60)
    print("VEP REST API Annotation Tool")
    print("=" * 60)
    print(f"Input: {vcf_file}")
    print(f"Output: {output_file}")
    print("=" * 60)
    
    # 解析 VCF
    print("Parsing VCF file...")
    variants = parse_vcf(vcf_file)
    print(f"Found {len(variants)} variants")
    
    # 批量调用 API (每批最多 200 个，API 限制)
    batch_size = 200
    all_results = []
    
    print(f"\nCalling VEP API (batch size: {batch_size})...")
    start_time = time.time()
    
    for i in range(0, len(variants), batch_size):
        batch = variants[i:i + batch_size]
        print(f"  Processing batch {i // batch_size + 1}/{(len(variants) - 1) // batch_size + 1} ({len(batch)} variants)...")
        
        results = call_vep_api(batch)
        
        if results:
            # 将结果与原始变异匹配
            result_map = {}
            for r in results:
                if 'input' in r:
                    result_map[r['input']] = r
            
            for v in batch:
                vep_result = result_map.get(v['vep_string'])
                annotations = extract_annotations(vep_result)
                all_results.append({
                    'variant': v,
                    'annotations': annotations,
                    'raw': vep_result
                })
        else:
            # API 调用失败，填充空注释
            for v in batch:
                all_results.append({
                    'variant': v,
                    'annotations': extract_annotations(None),
                    'raw': None
                })
        
        # 避免 API 限流
        time.sleep(1)
    
    elapsed = time.time() - start_time
    print(f"\nAPI calls completed in {elapsed:.1f}s")
    
    # 输出结果
    print(f"\nWriting results to {output_file}...")
    with open(output_file, 'w') as f:
        # 写入表头
        headers = [
            "CHROM", "POS", "REF", "ALT",
            "ClinVar_Significance", "ClinVar_Trait",
            "AlphaMissense", "SpliceAI_Max_Delta",
            "CADD_phred", "REVEL", "phyloP100way"
        ]
        f.write('\t'.join(headers) + '\n')
        
        # 写入数据
        for r in all_results:
            v = r['variant']
            a = r['annotations']
            row = [
                v['chrom'], str(v['pos']), v['ref'], v['alt'],
                a['clinvar_clnsig'], a['clinvar_trait'],
                a['alpha_missense'], a['spliceai_max'],
                a['cadd_phred'], a['revel'], a['phylop100way']
            ]
            f.write('\t'.join(row) + '\n')
    
    print("=" * 60)
    print("Annotation Complete!")
    print(f"Total variants: {len(all_results)}")
    print(f"Output file: {output_file}")
    print("=" * 60)
    
    # 统计注释覆盖率
    stats = {
        'clinvar': sum(1 for r in all_results if r['annotations']['clinvar_clnsig'] != '.'),
        'alpha_missense': sum(1 for r in all_results if r['annotations']['alpha_missense'] != '.'),
        'spliceai': sum(1 for r in all_results if r['annotations']['spliceai_max'] != '.'),
        'cadd': sum(1 for r in all_results if r['annotations']['cadd_phred'] != '.'),
        'revel': sum(1 for r in all_results if r['annotations']['revel'] != '.'),
        'phylop': sum(1 for r in all_results if r['annotations']['phylop100way'] != '.')
    }
    
    print("\nAnnotation Coverage:")
    for key, count in stats.items():
        print(f"  {key}: {count}/{len(all_results)} ({count/len(all_results)*100:.1f}%)")

if __name__ == "__main__":
    main()
