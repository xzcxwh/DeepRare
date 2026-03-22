package pipeline

import (
	"fmt"
	"io"
	"os/exec"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step4Result Step4 的运行结果
type Step4Result struct {
	InputCount    int
	OutputCount   int // 保留在 vcf4 中的变异
	ClinVarCount  int // 进入决赛圈的 ClinVar 变异
	DeletedCount  int // 被删除的良性变异
	Duration      time.Duration
}

// Step4ClinVarAnnotation 执行 Step4: ClinVar 注释
func Step4ClinVarAnnotation(cfg *config.Config) (*Step4Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 4: ClinVar 注释")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 打开输入文件
	parser, err := vcf.NewParser(cfg.VCF3)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

	// 创建输出文件 (vcf4 - 继续处理)
	writer4, err := vcf.NewWriter(cfg.VCF4)
	if err != nil {
		return nil, err
	}
	defer writer4.Close()

	// 创建决赛圈文件 (final_clinvar - 致病变异)
	writerClinVar, err := vcf.NewWriter(cfg.FinalClinVar)
	if err != nil {
		return nil, err
	}
	defer writerClinVar.Close()

	// 写入 header
	writer4.WriteHeader(parser.Header)
	writerClinVar.WriteHeader(parser.Header)

	result := &Step4Result{}

	fmt.Println("\n  查询 ClinVar 注释...")

	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		result.InputCount++

		// 使用 tabix 查询 dbNSFP 中的 ClinVar 信息
		clinvarInfo := queryClinVar(cfg.DBNSFPFile, v)

		if clinvarInfo != nil {
			// 添加 ClinVar 注释
			v.SetInfo("clinvar_id", clinvarInfo["clinvar_id"])
			v.SetInfo("clinvar_clnsig", clinvarInfo["clinvar_clnsig"])
			v.SetInfo("clinvar_trait", clinvarInfo["clinvar_trait"])
			v.SetInfo("clinvar_review", clinvarInfo["clinvar_review"])
			v.SetInfo("clinvar_hgvs", clinvarInfo["clinvar_hgvs"])
			v.SetInfo("clinvar_var_source", clinvarInfo["clinvar_var_source"])
			v.SetInfo("clinvar_MedGen_id", clinvarInfo["clinvar_MedGen_id"])
			v.SetInfo("clinvar_OMIM_id", clinvarInfo["clinvar_OMIM_id"])
			v.SetInfo("clinvar_Orphanet_id", clinvarInfo["clinvar_Orphanet_id"])

			clnsig := clinvarInfo["clinvar_clnsig"]  // 保持原始大小写
			review := clinvarInfo["clinvar_review"]

			// 判断是否进入决赛圈（与 Python 版本一致）
			// 条件：致病变异 + (专家评审 或 多个提交者无冲突)
			if shouldBeInFinalClinvar(clnsig, review) {
				result.ClinVarCount++
				writerClinVar.WriteVariant(v)
				fmt.Printf("  变异 %s:%d %s>%s 进入决赛圈: %s, %s\n",
					v.Chrom, v.Pos, v.Ref, v.Alt, clnsig, review)
				continue
			} else if shouldBeDeleted(clnsig, review) {
				// 良性变异 + 专家评审被删除
				result.DeletedCount++
				fmt.Printf("  变异 %s:%d %s>%s 被删除: %s, %s\n",
					v.Chrom, v.Pos, v.Ref, v.Alt, clnsig, review)
				continue
			}
		}

		// 其他变异继续处理
		result.OutputCount++
		writer4.WriteVariant(v)
	}

	result.Duration = time.Since(startTime)

	fmt.Printf("\n✅ Step 4 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   进入决赛圈 (final_clinvar.vcf): %d\n", result.ClinVarCount)
	fmt.Printf("   保留在 vcf4 中: %d\n", result.OutputCount)
	fmt.Printf("   删除的良性变异: %d\n", result.DeletedCount)

	return result, nil
}

// queryClinVar 使用 tabix 查询 ClinVar 信息
func queryClinVar(dbPath string, v *vcf.Variant) map[string]string {
	chromNum := v.Chrom
	if strings.HasPrefix(chromNum, "chr") {
		chromNum = chromNum[3:]
	}

	region := fmt.Sprintf("%s:%d-%d", chromNum, v.Pos, v.Pos)
	cmd := exec.Command("tabix", dbPath, region)
	output, err := cmd.Output()
	if err != nil || len(output) == 0 {
		return nil
	}

	// dbNSFP 列索引 (根据实际 header):
	// 3=ref, 4=alt, 10=clinvar_id, 11=clinvar_clnsig, 12=clinvar_trait
	// 13=clinvar_review, 14=clinvar_hgvs, 15=clinvar_var_source
	// 16=clinvar_MedGen_id, 17=clinvar_OMIM_id, 18=clinvar_Orphanet_id
	
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for _, line := range lines {
		fields := strings.Split(line, "\t")
		if len(fields) < 19 {
			continue
		}

		// 检查 ref/alt 匹配 (列 3 和 4，0-indexed 为 2 和 3)
		ref := fields[2]
		alt := fields[3]

		if ref != v.Ref || alt != v.Alt {
			continue
		}

		// 提取 ClinVar 字段
		result := map[string]string{
			"clinvar_id":         getField(fields, 9),   // 列 10 (0-indexed: 9)
			"clinvar_clnsig":     getField(fields, 10),  // 列 11
			"clinvar_trait":      getField(fields, 11),  // 列 12
			"clinvar_review":     getField(fields, 12),  // 列 13
			"clinvar_hgvs":       getField(fields, 13),  // 列 14
			"clinvar_var_source": getField(fields, 14),  // 列 15
			"clinvar_MedGen_id":  getField(fields, 15),  // 列 16
			"clinvar_OMIM_id":    getField(fields, 16),  // 列 17
			"clinvar_Orphanet_id": getField(fields, 17), // 列 18
		}

		// 检查是否有 ClinVar 信息
		if result["clinvar_id"] != "" && result["clinvar_id"] != "." {
			return result
		}
	}

	return nil
}

func getField(fields []string, idx int) string {
	if idx < len(fields) {
		val := fields[idx]
		if val != "." && val != "" {
			return val
		}
	}
	return ""
}

// shouldBeInFinalClinvar 判断是否应该进入 final_clinvar.vcf（与 Python 版本一致）
// 条件（满足任一即可）：
// 1. pathogenic 变异且 review 为 practice_guideline 或 reviewed_by_expert_panel
// 2. pathogenic 变异且 review 包含 multiple_submitters 且 no_conflicts
func shouldBeInFinalClinvar(clnsig, review string) bool {
	if clnsig == "" || clnsig == "." || review == "" || review == "." {
		return false
	}

	// 只接受确切的致病分类
	pathogenicTerms := []string{
		"Pathogenic",
		"Pathogenic/Likely_pathogenic",
		"Pathogenic|drug_response",
	}

	isPathogenic := false
	for _, term := range pathogenicTerms {
		if clnsig == term {
			isPathogenic = true
			break
		}
	}

	if !isPathogenic {
		return false
	}

	// 检查评审状态
	reviewLower := strings.ToLower(review)

	// 条件1：专家审核
	hasExpertReview := strings.Contains(reviewLower, "practice_guideline") ||
		strings.Contains(reviewLower, "reviewed_by_expert_panel")

	// 条件2：多个提交者且无冲突
	multipleNoConflict := strings.Contains(reviewLower, "multiple_submitters") &&
		strings.Contains(reviewLower, "no_conflicts")

	return hasExpertReview || multipleNoConflict
}

// shouldBeDeleted 判断是否应该从 vcf4 中删除（与 Python 版本一致）
// benign 变异且 review 为 practice_guideline 或 reviewed_by_expert_panel
func shouldBeDeleted(clnsig, review string) bool {
	if clnsig == "" || clnsig == "." || review == "" || review == "." {
		return false
	}

	benignTerms := []string{
		"Likely_benign",
		"Benign",
		"Benign/Likely_benign",
	}

	isBenign := false
	for _, term := range benignTerms {
		if clnsig == term {
			isBenign = true
			break
		}
	}

	if !isBenign {
		return false
	}

	// 检查评审状态
	reviewLower := strings.ToLower(review)
	return strings.Contains(reviewLower, "practice_guideline") ||
		strings.Contains(reviewLower, "reviewed_by_expert_panel")
}

// isPathogenic 判断是否为致病变异（保留用于其他地方）
func isPathogenic(clnsig string) bool {
	pathogenicTerms := []string{
		"Pathogenic",
		"Pathogenic/Likely_pathogenic",
		"Pathogenic|drug_response",
	}
	for _, term := range pathogenicTerms {
		if clnsig == term {
			return true
		}
	}
	return false
}

// isBenign 判断是否为良性变异（保留用于其他地方）
func isBenign(clnsig string) bool {
	benignTerms := []string{
		"Likely_benign",
		"Benign",
		"Benign/Likely_benign",
	}
	for _, term := range benignTerms {
		if clnsig == term {
			return true
		}
	}
	return false
}

// isReviewed 判断是否经过专家评审
func isReviewed(review string) bool {
	review = strings.ToLower(review)
	return strings.Contains(review, "expert") ||
		strings.Contains(review, "practice_guideline") ||
		strings.Contains(review, "multiple_submitters")
}

