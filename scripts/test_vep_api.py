#!/usr/bin/env python3
"""测试 VEP API 支持的插件和返回的字段结构 - 使用正确的参数"""
import requests
import json

VEP_API_URL = "https://rest.ensembl.org/vep/human/region"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# 使用一个已知的 missense 变异来获取 AlphaMissense/REVEL 分数
# BRCA1 p.Ile1854Val (可能有 REVEL 分数)
payload = {
    "variants": [
        "17 43071077 43071077 T/C 1",  # BRCA1 rs1799966
    ],
    # VEP REST API 支持的参数 - 注意布尔值用 1
    "AlphaMissense": 1,       # AlphaMissense 预测
    "CADD": 1,                # CADD 分数 (snv_indels)
    "REVEL": 1,               # REVEL 分数
    "SpliceAI": 1,            # SpliceAI 预测 (1 = canonical transcripts)
    "dbNSFP": "phyloP100way_vertebrate,REVEL_score,CADD_phred",  # dbNSFP 字段
    "dbscSNV": 1,             # dbscSNV 剪接预测
    "canonical": 1,
    "hgvs": 1,
    "protein": 1
}

print("发送请求 (使用正确的参数格式)...")
response = requests.post(VEP_API_URL, headers=HEADERS, json=payload, timeout=60)
result = response.json()

print(f"返回 {len(result)} 个结果\n")

# 打印完整的 transcript_consequences
print("=== transcript_consequences 完整内容 ===")
if result and 'transcript_consequences' in result[0]:
    for i, tc in enumerate(result[0]['transcript_consequences'][:3]):
        print(f"\n--- Transcript {i+1}: {tc.get('transcript_id', 'N/A')} ---")
        for k, v in sorted(tc.items()):
            print(f"  {k}: {v}")

# 检查顶层是否有 cadd
print("\n=== 顶层字段 ===")
for k in result[0].keys():
    if k not in ['transcript_consequences', 'colocated_variants']:
        print(f"  {k}: {result[0][k]}")

# 递归查找所有包含特定关键字的字段
def find_fields_with_value(obj, keywords, prefix=""):
    """递归查找包含特定关键字的字段"""
    results = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            full_key = f"{prefix}.{k}" if prefix else k
            k_lower = k.lower()
            if any(kw in k_lower for kw in keywords):
                results.append((full_key, v))
            results.extend(find_fields_with_value(v, keywords, full_key))
    elif isinstance(obj, list) and obj:
        results.extend(find_fields_with_value(obj[0], keywords, prefix + "[0]"))
    return results

print("\n=== 查找所有注释相关字段 ===")
keywords = ['am_', 'alpha', 'revel', 'cadd', 'splice', 'phylo', 'conservation']
found = find_fields_with_value(result, keywords)
for key, value in found:
    print(f"  {key}: {value}")
