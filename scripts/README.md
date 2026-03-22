# HERITA 变异筛选流水线

## 版本说明

本项目包含两个版本的流水线：

### 1. 原始版本 (scripts/original/)
- **耗时**: ~48秒
- **结构**: 9个独立步骤
- **运行**: `python scripts/original/run_pipeline.py`

步骤流程:
1. Step 1: gnomAD频率注释与过滤 (vcfanno)
2. Step 2: gnomAD Flag注释
3. Step 3: 测序质量过滤
4. Step 4: ClinVar注释
5. Step 5: dbscSNV注释
6. Step 6: dbNSFP多评分注释与过滤
7. Step 7: 最终Top50筛选
8. Step 8: 最终整合
9. Step 9: InterVar ACMG分类

### 2. 优化版本 (scripts/optimized/)
- **耗时**: ~28秒 (提升42%)
- **结构**: 5个阶段，单脚本执行
- **运行**: `python scripts/optimized/run_pipeline_extreme.py`

优化策略:
1. 使用位置集合(set)快速过滤，避免重复查找
2. 流式扫描dbNSFP，只保留匹配的记录
3. 只查询候选变异的dbscSNV，而非全量加载
4. 一次遍历完成所有过滤、注释和排序
5. 合并多个处理步骤，减少I/O操作

## 输出结果

两个版本的最终输出相同:
- `results/final_intervar_classified.vcf` - ACMG分类后的最终结果

## 性能对比

| 版本 | 耗时 | 提升 |
|------|------|------|
| 原始版本 | ~48秒 | - |
| 优化版本 | ~28秒 | 42% |

## 依赖工具

- Python 3.8+
- vcfanno (原始版本需要)
- bcftools
