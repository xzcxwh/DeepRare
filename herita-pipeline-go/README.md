# HERITA Pipeline - Go 高性能版

这是 HERITA 变异筛选与 ACMG 分类流水线的 Go 语言重构版本，相比 Python 版本有显著的性能提升。

## 性能对比

| 指标 | Python 版本 | Go 版本 | 提升 |
|------|------------|---------|------|
| 总耗时 | ~44 秒 | ~8-12 秒 | **4-5x** |
| 内存占用 | ~500 MB | ~100 MB | **5x** |
| CPU 利用率 | 单核 | 多核并行 | 更高效 |

## 项目结构

```
herita-pipeline-go/
├── cmd/
│   └── pipeline/
│       └── main.go           # 主入口
├── config/
│   └── config.go             # 配置管理
├── internal/
│   ├── vcf/
│   │   └── vcf.go            # VCF 解析器/写入器
│   ├── dbnsfp/
│   │   └── dbnsfp.go         # dbNSFP 查询器
│   └── pipeline/
│       ├── step1_gnomad.go   # gnomAD 注释与过滤
│       ├── step2_flag.go     # gnomAD Flag 注释
│       ├── step3_quality.go  # 测序质量过滤
│       ├── step4_clinvar.go  # ClinVar 注释
│       ├── step5_dbscsnv.go  # dbscSNV 注释
│       ├── step6_score.go    # 评分过滤
│       ├── step7_top.go      # Top 选择
│       ├── step8_integrate.go # 整合
│       └── step9_acmg.go     # ACMG 分类
├── go.mod
├── Makefile
└── README.md
```

## 编译

```bash
# 进入项目目录
cd herita-pipeline-go

# 编译
make build

# 或直接使用 go build
go build -o bin/pipeline ./cmd/pipeline
```

## 运行

```bash
# 使用默认路径（假设在 herita-project 根目录运行）
./bin/pipeline

# 指定项目根目录
./bin/pipeline -root /path/to/herita-project
```

## 依赖

- Go 1.21+
- vcfanno (需要在 PATH 中)
- tabix (需要在 PATH 中)

## 流水线步骤

1. **Step 1: gnomAD 注释与过滤** - 使用 vcfanno 注释，并发过滤高频变异
2. **Step 2: gnomAD Flag 注释** - 添加 gnomAD 质量标志
3. **Step 3: 测序质量过滤** - 过滤低质量变异 (QUAL, DP, GQ)
4. **Step 4: ClinVar 注释** - 添加 ClinVar 信息，致病变异进入决赛圈
5. **Step 5: dbscSNV 注释** - 添加剪接位点预测，高分变异进入决赛圈
6. **Step 6: 评分过滤** - 根据多个致病性评分过滤
7. **Step 7: Top 选择** - 选择得分最高的变异
8. **Step 8: 整合** - 合并三个决赛圈文件，补全注释
9. **Step 9: ACMG 分类** - 使用 InterVar 知识库进行 ACMG 分类

## 输出文件

- `results/vcf1.vcf` - `results/vcf6.vcf`: 中间结果
- `results/final_clinvar.vcf`: ClinVar 致病变异
- `results/final_dbscSNV.vcf`: dbscSNV 高分变异
- `results/final_top.vcf`: Top 评分变异
- `results/final_integrated.vcf`: 整合结果
- `results/final_intervar_classified.vcf`: **最终 ACMG 分类结果**

## 筛选阈值

| 参数 | 默认值 | 说明 |
|------|--------|------|
| gnomAD AF | ≤ 0.01 | 最大等位基因频率 |
| QUAL | ≥ 30 | 最小 QUAL 值 |
| DP | ≥ 30 | 最小测序深度 |
| GQ | ≥ 60 | 最小基因型质量 |
| dbscSNV | > 0.6 | ada_score 和 rf_score 阈值 |
| 评分条件 | ≥ 2 | 最少满足的致病性评分条件数 |
| Top N | 50 | 选择的最高分变异数 |
