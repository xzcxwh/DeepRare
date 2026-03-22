# ACMG 致病性分类报告

## 概述

本报告基于 ACMG/AMP 2015 指南对变异进行自动化致病性评级。

**注意**: 自动化分类仅供参考，临床使用需专业人员复核。

## 分类结果统计

| 分类 | 数量 | 百分比 |
|------|------|--------|
| Pathogenic | 0 | 0.0% |
| Likely_Pathogenic | 1 | 1.4% |
| VUS | 71 | 98.6% |
| Likely_Benign | 0 | 0.0% |
| Benign | 0 | 0.0% |
| **Total** | **72** | **100%** |

## 使用的证据类型

### 致病性证据

| 代码 | 强度 | 说明 |
|------|------|------|
| PS1 | Strong | ClinVar 致病变异（专家审核） |
| PM2 | Moderate | 人群频率 < 0.01% 或缺失 |
| PP3 | Supporting | ≥3 个计算工具预测有害 |
| PP3_splicing | Supporting | 剪接预测支持有害 |
| PP5 | Supporting | ClinVar 报告致病 |
| PP5_strong | Moderate | ClinVar 致病（多提交者无冲突） |

### 良性证据

| 代码 | 强度 | 说明 |
|------|------|------|
| BA1 | Stand-alone | 人群频率 > 5% |
| BS1 | Strong | 人群频率 > 1% |
| BP4 | Supporting | ≥3 个计算工具预测良性 |
| BP4_splicing | Supporting | 剪接预测支持良性 |
| BP6 | Supporting | ClinVar 报告良性 |
| BP6_strong | Strong | ClinVar 良性（专家审核） |

## 计算预测工具阈值

| 工具 | 有害阈值 | 良性阈值 |
|------|----------|----------|
| AlphaMissense | ≥ 0.564 | < 0.34 |
| REVEL | ≥ 0.5 | < 0.25 |
| MetaRNN | ≥ 0.5 | < 0.25 |
| ClinPred | ≥ 0.5 | < 0.25 |
| CADD_phred | ≥ 20 | < 10 |
| PrimateAI | ≥ 0.7 | < 0.3 |
| ESM1b | ≤ -7 | > -3 |
| dbscSNV (ADA+RF) | ≥ 0.6 | < 0.3 |

## 输出文件

- 分类结果: `/Volumes/T9/herita-project/results/final_acmg_classified.vcf`
- 新增字段:
  - `ACMG_CLASS`: 分类结果
  - `ACMG_EVIDENCE`: 支持分类的证据代码

## 参考文献

Richards S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med. 2015;17(5):405-424.
