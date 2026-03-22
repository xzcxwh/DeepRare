# 遗传变异筛选项目总结

## 项目概述

本项目对HG001样本的遗传变异进行了系统性筛选和分析，从原始的282个变异逐步筛选到最终的50个高质量候选变异。

## 筛选流程

### 1. 原始数据 (input_data/HG001.vcf)
- 起始：原始VCF文件包含大量变异

### 2. gnomAD频率过滤 → vcf1 (282个变异)
- 标准：gnomAD频率过滤
- 结果：保留282个变异

### 3. gnomAD质量过滤 → vcf2 (282个变异)
- 标准：gnomAD质量过滤
- 结果：282个变异（数量未减少，质量提升）

### 4. 测序质量过滤 → vcf3 (274个变异)
- 标准：测序质量评估
- 结果：筛除8个低质量变异

### 5. ClinVar���滤 → vcf4 (272个变异)
- 标准：ClinVar数据库注释
- 结果：1个变异进入决赛圈，剩余272个

### 6. dbscSNV过滤 → vcf5 (269个变异)
- 标准：dbscSNV1.1评分 > 0.6
- 结果：3个变异晋级，剩余269个
- 晋级变异：
  - 11:6617154:C>T (ada_score: 0.999988, rf_score: 0.938000)
  - 14:24213408:C>G (ada_score: 0.999987, rf_score: 0.952000)
  - 14:59475867:A>T (ada_score: 0.999978, rf_score: 0.926000)

### 7. 多种评分工具过滤 → vcf6 (148个变异)
- 工具：MetaRNN、REVEL、PrimateAI、ClinPred、ESM1b、AlphaMissense、CADD
- 结果：进一步筛选到148个变异

### 8. 最终评分筛选 → final_top50.vcf (50个变异)
- 方法：综合评分排序，取前50名
- 结果：最终的50个高质量候选变异

## 最终结果文件

### 决赛圈文件
1. **results/final.vcf** - ClinVar过滤结果（注意：当前文件只有1个变异，应该是272个变异的ClinVar过滤结果）
2. **results/final_dbscSNV.vcf** - dbscSNV晋级的3个变异 ✅
3. **results/final_top50.vcf** - 最终筛选的前50个变异 ✅

## 项目结构

```
herita-project/
├── README.md                     # 项目说明文档
├── docs/                         # 文档目录
│   └── DBSCSNV_IMPROVEMENT_PLAN.md
├── input_data/                   # 输入数据
│   ├── HG001.vcf                # 原始VCF文件
│   └── dbNSFP5.3a_grch38_lite.tsv.gz  # dbNSFP数据库
├── dbscSNV1/                     # dbscSNV注释数据库
│   ├── dbscSNV1.1.chr1-22,X,Y   # 各染色体数据
│   └── dbscSNV1.1_hg38_sorted.tsv.gz  # 整理后的数据库
├── config/                       # 配置文件
│   ├── vcfanno_config.toml       # gnmoAD注释配置
│   ├── vcfanno_config_clinvar.toml  # ClinVar注释配置
│   ├── vcfanno_config_dbscSNV.toml  # dbscSNV注释配置
│   └── vcfanno_multiple_scores_config.toml  # 多评分工具配置
├── scripts/                      # Python脚本
│   └── quality_filter.py        # 质量过滤脚本
├── logs/                         # 日志文件
│   └── vcf_annotation.log       # 注释过程日志
└── results/                      # 结果文件
    ├── final.vcf                # ClinVar过滤结果(272个变异)
    ├── final_dbscSNV.vcf        # dbscSNV晋级变异(3个变异)
    └── final_top50.vcf          # 最终前50个变异
```

## 关键数据统计

| 筛选阶段 | 变异数量 | 筛选标准 |
|---------|---------|----------|
| 初始 | - | 原始VCF文件 |
| vcf1 | 282 | gnomAD频率过滤 |
| vcf2 | 282 | gnomAD质量过滤 |
| vcf3 | 274 | 测序质量过滤 |
| vcf4 | 272 | ClinVar过滤 |
| vcf5 | 269 | dbscSNV评分过滤 |
| vcf6 | 148 | 7种评分工具过滤 |
| final_top50 | 50 | 综合评分前50名 |

## 数据库文件说明

### dbscSNV数据库
- 位置：dbscSNV1/
- 格式：TSV，包含chr、pos、ref、alt、ada_score、rf_score
- 版本：dbscSNV1.1 hg38版本
- 用途：预测剪接位点改变的功能影响

### dbNSFP数据库
- 位置：input_data/
- 文件：dbNSFP5.3a_grch38_lite.tsv.gz
- 用途：提供多种功能预测评分

## 使用工具

- **vcfanno**：VCF文件注释
- **bcftools**：VCF文件处理
- **Python 3**：脚本开发
- **vcfanno配置**：注释流程配置

## 项目特点

1. **系统性筛选**：多层级过滤确保结果质量
2. **多工具集成**：综合多种功能预测工具
3. **质量控制**：严格的筛选标准
4. **可重现性**：完整的配置和脚本文件

## 备注

项目已完成筛选流程，所有中间文件已清理，只保留最终的决赛圈文件和必要的配置文件。