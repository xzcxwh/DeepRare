# HERITA 变异注释与筛选流水线 - 运行总结

## 执行日期: 2024-12-10

## 输入文件
- **文件**: `input_data/HG001_cds.vcf`
- **初始变异数**: 21,716

## 流水线执行步骤

| 步骤 | 名称 | 输入 | 输出 | 耗时 | 结果 |
|:---:|:---|:---:|:---:|:---:|:---|
| 1 | gnomAD 频率过滤 | 21,716 | 741 | 14m19s | AF ≤ 0.01 保留 3.41% |
| 2 | ClinVar 注释 | 741 | 741 | 1s | 276个获得注释 (37.25%) |
| 3 | ClinVar 致病性分流 | 741 | 4 + 640 | 0s | 4个P/LP→决赛圈, 97个B/LB移除 |
| 4 | snpEff HIGH 过滤 | 640 | 22 + 618 | 33s | 22个HIGH→决赛圈 |
| 5 | dbscSNV 剪接过滤 | 618 | 0 + 618 | 26s | 无高分剪接变异 |
| 6 | VEP API 注释 | 618 | 2 + 616 | 2m29s | 2个SpliceAI>0.5→决赛圈 |
| 7 | dbNSFP 深度注释 | 616 | 616 | 0s | 261个获得多预测器分数 |
| 8 | 综合打分筛选 | 616 | 50 | 0s | 选出TOP 50 |
| 9 | 合并决赛圈 | 28 + 50 | 78 | 0s | 去重合并 |

**总耗时**: ~17分48秒

## 决赛圈变异来源

| 来源 | 数量 | 说明 |
|:---|:---:|:---|
| Pathogenic/Likely Pathogenic | 4 | ClinVar 致病/可能致病 |
| HIGH Impact (snpEff) | 22 | 移码/终止密码子/剪接位点 |
| SpliceAI > 0.5 | 2 | 高剪接影响 |
| TOP 50 综合评分 | 50 | 多预测器加权打分 |
| **合并后总计** | **78** | 去重后决赛圈 |

## 决赛圈 HIGH Impact 变异类型

| 类型 | 数量 | 示例 |
|:---|:---:|:---|
| frameshift_variant | 12 | CDK11B, CELA1, MYT1L, DLG1 |
| stop_gained | 6 | ECHDC2, OR4C15, NRIP2, TLR10 |
| splice_donor/acceptor | 4 | AL096870.2, AC003002.2, PIWIL3 |

## TOP 5 综合评分变异

| 排名 | 位置 | 基因 | 变异 | 总分 | 致病性 | 保守性 |
|:---:|:---|:---|:---|:---:|:---:|:---:|
| 1 | chr1:65218865 | AK4 | G>C | 56.00 | 37 | 13 |
| 2 | chr8:11845732 | CTSB | C>G | 56.00 | 37 | 11 |
| 3 | chr3:14821133 | FGD5 | G>A | 54.00 | 35 | 13 |
| 4 | chr22:32408189 | RTCB | G>T | 53.00 | 34 | 13 |
| 5 | chr9:128713304 | PKN3 | G>C | 53.00 | 31 | 11 |

## 输出文件

```
results_pipeline/
├── HG001_cds_gnomad_filtered.vcf       # Step 1: gnomAD 过滤后 (741)
├── HG001_cds_clinvar_annotated.vcf     # Step 2: ClinVar 注释后
├── HG001_cds_next_level.vcf            # Step 3: 待进一步分析 (640)
├── HG001_cds_after_snpeff.vcf          # Step 4: snpEff 过滤后 (618)
├── HG001_cds_after_dbscsnv.vcf         # Step 5: dbscSNV 过滤后 (618)
├── HG001_cds_after_vep.vcf             # Step 6: VEP 注释后 (616)
├── HG001_cds_dbnsfp_annotated.vcf      # Step 7: dbNSFP 注释后 (616)
├── HG001_cds_final_round.vcf           # 决赛圈累积 (28)
├── HG001_cds_top50_finalists.vcf       # Step 8: TOP 50 (50)
├── HG001_cds_finalists_merged.vcf      # Step 9: 最终合并 (78)
└── HG001_cds_scoring_report.txt        # 打分报告
```

## 流水线配置

```bash
AF 阈值:       0.01 (1%)
SpliceAI 阈值: 0.5
TOP N:         50
线程数:        4
```

## 使用的工具

| 工具 | 版本 | 用途 |
|:---|:---|:---|
| gnomad_grpmax_filter | Go | 群体频率过滤 |
| clinvar_annotate | Go | ClinVar 注释 |
| clinvar_filter | Go | 致病性分流 |
| snpeff_filter | Go + Java 25 | 功能影响注释 |
| dbscsnv_filter | Go + vcfanno | 剪接预测 |
| vep_annotate | Go + REST API | Ensembl VEP 注释 |
| dbnsfp_annotate | Go + vcfanno | 多预测器注释 |
| variant_scorer.py | Python 3 | 综合打分 |
| finalists_merger | Go | 决赛圈合并 |

## 注释字段

### ClinVar
- CLNSIG: 临床意义
- CLNREVSTAT: 审核状态

### snpEff
- ANN: 功能注释 (Impact: HIGH/MODERATE/LOW/MODIFIER)

### VEP (通过 REST API)
- AlphaMissense
- SpliceAI (DS_AG, DS_AL, DS_DG, DS_DL)
- CADD
- REVEL
- phyloP

### dbNSFP
- AlphaMissense_score
- REVEL_score
- MetaRNN_score
- PrimateAI_score
- CADD_phred
- ClinPred_score
- ESM1b_score

## 打分公式

**总分 = 致病性分数 + 保守性分数 + 功能分数 + 剪接分数 + 临床分数 + 频率分数**

- 致病性: AlphaMissense, REVEL, MetaRNN, ClinPred, ESM1b (加权)
- 保守性: phyloP, CADD (加权)
- 功能: PrimateAI, 功能影响类型
- 剪接: SpliceAI max score
- 临床: ClinVar 意义 + 审核星级
- 频率: gnomAD AF 惩罚

## 运行命令

```bash
./run_pipeline.sh input_data/HG001_cds.vcf -o results_pipeline -v
```

## 结论

从 21,716 个 CDS 变异中，通过 9 步流水线筛选出 78 个高优先级候选变异，包括：
- 4 个 ClinVar 已知致病/可能致病
- 22 个功能严重影响 (HIGH Impact)
- 2 个高剪接影响
- 50 个综合评分最高的 missense 变异

这些变异值得进一步临床审查和功能验证。
