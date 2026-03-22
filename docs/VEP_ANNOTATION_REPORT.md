# VEP REST API 注释报告

## 概述

本报告总结了对 HG001 样本 179 个 HIGH/MODERATE 影响变异的 VEP REST API 注释结果。

## 输入输出

| 项目 | 值 |
|------|-----|
| 输入文件 | `results/HG001_543_high_moderate.vcf` |
| 输出文件 | `results/HG001_179_vep_annotated.tsv` |
| 变异数量 | 179 |
| API 调用耗时 | 40.3 秒 |

## 注释覆盖率

| 注释类型 | 有注释 | 总数 | 覆盖率 |
|----------|--------|------|--------|
| ClinVar | 45 | 179 | 25.1% |
| AlphaMissense | 104 | 179 | 58.1% |
| SpliceAI | 133 | 179 | 74.3% |
| CADD | 163 | 179 | 91.1% |
| REVEL | 104 | 179 | 58.1% |
| phyloP100way | 140 | 179 | 78.2% |

## ClinVar 临床意义分布

| 分类 | 数量 |
|------|------|
| benign | 19 |
| likely_benign | 7 |
| benign,likely_benign | 7 |
| uncertain_significance | 5 |
| 混合/其他 | 7 |
| 无注释 | 134 |

## AlphaMissense 分布 (104 个有分数)

AlphaMissense 使用以下阈值判断致病性：
- < 0.34: likely benign
- 0.34 - 0.564: ambiguous
- ≥ 0.564: likely pathogenic

| 分类 | 数量 |
|------|------|
| likely_benign (< 0.34) | 100 |
| ambiguous (0.34 - 0.564) | 2 |
| **likely_pathogenic (≥ 0.564)** | **2** |

## SpliceAI 分布 (133 个有分数)

SpliceAI delta score 表示剪接影响程度：
- 0: 无剪接影响
- 0.2 - 0.5: 中等影响
- ≥ 0.5: 高影响

| 分类 | 数量 |
|------|------|
| 0 (无影响) | 102 |
| 0 - 0.2 (低) | 29 |
| 0.2 - 0.5 (中等) | 1 |
| **≥ 0.5 (高)** | **1** |

## REVEL 分布 (104 个有分数)

REVEL 使用以下阈值：
- < 0.5: likely benign
- 0.5 - 0.7: uncertain
- ≥ 0.7: likely pathogenic

| 分类 | 数量 |
|------|------|
| likely_benign (< 0.5) | 102 |
| uncertain (0.5 - 0.7) | 1 |
| **likely_pathogenic (≥ 0.7)** | **1** |

## 潜在致病性变异

### 1. chr15:82151813 G>T (多重证据支持)

| 注释 | 值 | 解读 |
|------|-----|------|
| AlphaMissense | 0.6575 | **Likely pathogenic** |
| REVEL | 0.908 | **Likely pathogenic** |
| CADD | 25.0 | 高致病性分数 |
| phyloP100way | 9.495 | 高度保守 |
| ClinVar | uncertain/benign/likely_benign | 分类冲突 |

**注意**: 此变异在多个计算预测工具中显示可能致病性，但 ClinVar 有矛盾分类，值得进一步研究。

### 2. chr22:24728176 C>T (剪接影响)

| 注释 | 值 | 解读 |
|------|-----|------|
| SpliceAI | 0.86 | **高剪接影响** |
| CADD | 22.6 | 较高致病性分数 |
| ClinVar | . | 无 ClinVar 记录 |

**注意**: 此变异可能显著影响剪接，需要进一步功能验证。

### 3. chr1:13196958 G>C (AlphaMissense 预测)

| 注释 | 值 | 解读 |
|------|-----|------|
| AlphaMissense | 0.6913 | **Likely pathogenic** |
| CADD | 0.211 | 低 CADD 分数 |
| phyloP | -0.002 | 不保守 |
| ClinVar | . | 无 ClinVar 记录 |

**注意**: AlphaMissense 预测可能致病，但其他证据不支持，需谨慎解读。

## 输出文件格式

TSV 文件包含以下列：

1. **CHROM** - 染色体
2. **POS** - 位置
3. **REF** - 参考等位基因
4. **ALT** - 变异等位基因
5. **ClinVar_Significance** - ClinVar 临床意义
6. **ClinVar_Trait** - ClinVar 变异 ID
7. **AlphaMissense** - AlphaMissense 病理性分数
8. **SpliceAI_Max_Delta** - SpliceAI 最大 delta 分数 (max of DS_AG, DS_AL, DS_DG, DS_DL)
9. **CADD_phred** - CADD phred 分数
10. **REVEL** - REVEL 分数
11. **phyloP100way** - phyloP100way 保守性分数

## 结论

- 179 个 HIGH/MODERATE 影响变异中，大多数被预测为良性
- **2 个变异** 被 AlphaMissense 预测为可能致病
- **1 个变异** 有高 SpliceAI 分数 (>0.5)
- **1 个变异** (chr15:82151813) 有多重致病性证据，值得进一步研究
- 45 个变异 (25.1%) 有 ClinVar 记录，但无明确致病性分类

## 技术细节

- 使用 Ensembl VEP REST API (https://rest.ensembl.org/vep/human/region)
- 启用的插件: AlphaMissense, CADD, REVEL, SpliceAI, dbNSFP (phyloP100way_vertebrate)
- 参考基因组: GRCh38
- 注释时间: 2024年
