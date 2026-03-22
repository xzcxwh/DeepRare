# HERITA 流水线 ACMG 分类器集成报告

## 概述

已成功将 ACMG 分类步骤集成到 HERITA 变异注释与筛选流水线中，作为第 9 步运行。

## 流水线结构（9 步）

| 步骤 | 名称 | 输出文件 |
|------|------|----------|
| 1 | gnomAD 频率过滤 | `*_gnomad_filtered.vcf` |
| 2 | ClinVar 注释 | `*_clinvar_annotated.vcf` |
| 3 | ClinVar 致病性分流 | `*_final_round.vcf` / `*_next_level.vcf` |
| 4 | snpEff HIGH 过滤 | `*_after_snpeff.vcf` |
| 5 | VEP API 注释 | `*_after_vep.vcf` |
| 6 | dbNSFP 深度注释 | `*_dbnsfp_annotated.vcf` |
| 7 | 综合打分筛选 TOP 50 | `*_top50_finalists.vcf` |
| 8 | 合并决赛圈 | `*_finalists_merged.vcf` |
| **9** | **ACMG 分类** | `*_acmg_classified.vcf` |

## ACMG 分类器实现

### 脚本位置
`scripts/acmg_classifier.py`

### 知识库
使用 InterVar 项目的知识库 (`InterVar-master/intervardb/`)：

| 文件 | 用途 | 记录数 |
|------|------|--------|
| `PVS1.LOF.genes.hg38` | LOF 机制致病基因 | 6,213 |
| `PP2.genes.hg38` | 低良性错义变异率基因 | 931 |
| `BP1.genes.hg38` | 截断变异致病基因 | 1,449 |
| `PS1.AA.change.patho.hg38` | 已知致病 AA 变化 | 35,530 |
| `PM1_domains_with_benigns.hg38` | 功能域 | 1,724 |
| `PS4.variants.hg38` | GWAS 显著变异 | 1,159 |
| `BS2_hom_het.hg38` | 健康个体观察到的变异 | 5,665,742 |

### 实现的 ACMG 标准

#### 致病性证据
- **PVS1**: LOF 机制基因中的无义/移码/剪接变异
- **PS1**: 与已知致病变异相同的 AA 变化
- **PS3**: 多个预测器高分（AlphaMissense ≥0.9, REVEL ≥0.9）
- **PS4**: GWAS 显著富集变异
- **PM1**: 位于功能域
- **PM2**: 人群中缺失或极低频 (<0.0001)
- **PM4**: 框内缺失/插入或终止密码子丢失
- **PM5**: 已知致病位点的不同错义变化
- **PP2**: 低良性错义率基因中的错义变异
- **PP3**: 多个计算预测器支持致病（≥2 个）
- **PP5**: ClinVar 报告为致病

#### 良性证据
- **BA1**: 等位基因频率 > 5%
- **BS1**: 等位基因频率 > 1%
- **BS2**: 在健康个体中观察到
- **BP1**: 截断变异致病基因中的错义变异
- **BP4**: 多个计算预测器支持良性（≥2 个）
- **BP6**: ClinVar 报告为良性
- **BP7**: 同义变异且不影响剪接

### 分类规则

根据 ACMG/AMP 2015 指南：

**Pathogenic**:
- PVS1 + ≥1 PS
- PVS1 + ≥2 PM
- PVS1 + 1 PM + ≥1 PP
- PVS1 + ≥2 PP
- ≥2 PS
- 1 PS + ≥3 PM
- 1 PS + 2 PM + ≥2 PP
- 1 PS + 1 PM + ≥4 PP

**Likely_pathogenic**:
- PVS1 + 1 PM
- 1 PS + 1-2 PM
- 1 PS + ≥2 PP
- ≥3 PM
- 2 PM + ≥2 PP
- 1 PM + ≥4 PP

**Benign**:
- BA1
- ≥2 BS

**Likely_benign**:
- 1 BS + ≥1 BP
- ≥2 BP

**VUS**: 其他情况

## 测试结果

在 HG001_cds 数据集上运行（94 个决赛圈变异）：

```
ACMG分类统计
----------------------------------------
  Pathogenic: 0 (0.0%)
  Likely_pathogenic: 8 (8.5%)
  VUS: 80 (85.1%)
  Likely_benign: 6 (6.4%)
  Benign: 0 (0.0%)
  Conflicting: 0 (0.0%)
  总计: 94
```

### Likely_pathogenic 变异

| 位置 | 基因 | 变异类型 | ACMG 证据 |
|------|------|----------|-----------|
| chr1:65218865 | AK4 | 错义 | PS3,PM1,PM2,PP3 |
| chr2:1830817 | MYT1L | 移码 | PVS1,PM2 |
| chr3:197226031 | DLG1 | 移码 | PVS1,PM2 |
| chr6:31645145 | BAG6 | 移码 | PVS1,PM2 |
| chr12:47980970 | COL2A1 | 剪接 | PVS1,PM2 |
| chr13:77017068 | CLN5 | 移码 | PVS1,PM2 |
| chr15:23646089 | MAGEL2 | 移码 | PVS1,PM2 |
| chrX:45076787 | KDM6A | 移码 | PVS1,PM2 |

### 测试变异验证

4 个添加的测试变异分类结果：

| 变异 | 基因 | 分类 | 证据 |
|------|------|------|------|
| chrX:45076787 | KDM6A (移码) | ✅ Likely_pathogenic | PVS1,PM2 |
| chr15:23646089 | MAGEL2 (移码) | ✅ Likely_pathogenic | PVS1,PM2 |
| chr14:102042083 | DYNC1H1 (错义) | VUS | PM1,PM2,PP3 |
| chr12:47980970 | COL2A1 (剪接) | ✅ Likely_pathogenic | PVS1,PM2 |

## VCF 输出格式

添加到 INFO 字段的新注释：

```
##INFO=<ID=ACMG_Classification,Number=1,Type=String,Description="ACMG classification: Pathogenic, Likely_pathogenic, VUS, Likely_benign, Benign">
##INFO=<ID=ACMG_Evidence,Number=.,Type=String,Description="ACMG evidence codes">
```

示例：
```
...;ACMG_Classification=Likely_pathogenic;ACMG_Evidence=PVS1,PM2
```

## 使用方法

```bash
# 完整流水线运行
./herita-pipeline-go/bin/herita_pipeline \
  -input input_data/HG001_cds.vcf \
  -output results_pipeline \
  -gnomad gnomAD_dataset \
  -clinvar /path/to/clinvar.vcf.gz \
  -dbnsfp /path/to/dbNSFP4.9a_grch38.gz \
  -af 0.01 \
  -threads 8

# 仅运行 ACMG 分类步骤
./herita-pipeline-go/bin/herita_pipeline \
  ... \
  -skip gnomad,clinvar_annotate,clinvar_filter,snpeff,vep,dbnsfp,scoring,merge

# 单独运行 ACMG 分类器
python3 scripts/acmg_classifier.py \
  -i results_pipeline/HG001_cds_finalists_merged.vcf \
  -o results_pipeline/HG001_cds_acmg_classified.vcf \
  -r results_pipeline/HG001_cds_acmg_report.txt \
  --intervar-db InterVar-master/intervardb
```

## 文件变更

1. **新增**: `scripts/acmg_classifier.py` - ACMG 分类器 Python 脚本
2. **修改**: `herita-pipeline-go/cmd/herita_pipeline/main.go` - 添加步骤 9
3. **新增**: `docs/ACMG_CLASSIFIER_INTEGRATION.md` - 本文档

## 日期

2025-01-XX
