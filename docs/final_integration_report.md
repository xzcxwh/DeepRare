# 最终整合报告

本报告记录了三个决赛圈文件整合并补全注释的过程。

## 输入文件统计

- final_clinvar.vcf: 4 个变异
- final_dbscSNV.vcf: 2 个变异
- final_top.vcf: 66 个变异

## 整合结果

- 最终整合变异总数: 72 个
- 整合后文件: final_integrated.vcf

## 注释处理

- ✅ 对所有变异重新进行dbNSFP注释
- ✅ 对所有变异重新进行dbscSNV注释
- ✅ 清理重复注释值
- ✅ 移除冗余字段

## 注释值处理规则

### 取最大值的字段（致病性越高分数越大）:
- AlphaMissense_score, REVEL_score, MetaRNN_score
- PrimateAI_score, CADD_phred, ClinPred_score
- DBSCSNV_ADA, DBSCSNV_RF

### 取最小值的字段（影响越大分数越小）:
- ESM1b_score

### 取唯一值的字段:
- 基因名称、ClinVar信息、gnomAD频率等

## 保留的注释字段

- gnomAD4.1_joint_POPMAX_AF
- gnomAD4.1_joint
- DBSCSNV_ADA
- DBSCSNV_RF
- MetaRNN_score
- REVEL_score
- PrimateAI_score
- ClinPred_score
- ESM1b_score
- AlphaMissense_score
- CADD_phred
- aaref
- aaalt
- rs_dbSNP
- genename
- clinvar_id
- clinvar_clnsig
- clinvar_trait
- clinvar_review
- clinvar_hgvs
- clinvar_var_source
- clinvar_MedGen_id
- clinvar_OMIM_id
- clinvar_Orphanet_id

## 分析完成时间

2025年12月 3日 星期三 09时47分44秒 CST
