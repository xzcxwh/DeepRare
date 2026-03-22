package pipeline

import (
	"fmt"
	"io"
	"sort"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step7Result Step7 的运行结果
type Step7Result struct {
	InputCount  int
	OutputCount int
	MinScore    int
	MaxScore    int
	Duration    time.Duration
}

// ScoredVariant 带分数的变异
type ScoredVariant struct {
	Variant *vcf.Variant
	Score   int
}

// Step7TopSelection 执行 Step7: 最终 Top 变异选择
func Step7TopSelection(cfg *config.Config) (*Step7Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 7: 最终 Top 变异选择")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 读取并打分
	fmt.Println("\n  [1/2] 对变异打分...")
	parser, err := vcf.NewParser(cfg.VCF6)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

	var scoredVariants []ScoredVariant

	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		score := calculateScore(v)
		scoredVariants = append(scoredVariants, ScoredVariant{
			Variant: v,
			Score:   score,
		})
	}

	fmt.Printf("  共 %d 个变异参与排序\n", len(scoredVariants))

	// 按分数排序（降序）
	sort.Slice(scoredVariants, func(i, j int) bool {
		if scoredVariants[i].Score != scoredVariants[j].Score {
			return scoredVariants[i].Score > scoredVariants[j].Score
		}
		// 分数相同时按位置排序
		if scoredVariants[i].Variant.ChromNum() != scoredVariants[j].Variant.ChromNum() {
			return scoredVariants[i].Variant.ChromNum() < scoredVariants[j].Variant.ChromNum()
		}
		return scoredVariants[i].Variant.Pos < scoredVariants[j].Variant.Pos
	})

	// 选择 Top N（确保至少包含 TopN 个，包括同分的）
	fmt.Println("\n  [2/2] 选择高分变异...")

	cutoffScore := 0
	if len(scoredVariants) >= cfg.TopN {
		cutoffScore = scoredVariants[cfg.TopN-1].Score
	}

	var selected []ScoredVariant
	for _, sv := range scoredVariants {
		if sv.Score >= cutoffScore {
			selected = append(selected, sv)
		}
	}

	// 写入结果
	writer, err := vcf.NewWriter(cfg.FinalTop)
	if err != nil {
		return nil, err
	}
	defer writer.Close()

	writer.WriteHeader(parser.Header)

	for _, sv := range selected {
		writer.WriteVariant(sv.Variant)
	}

	result := &Step7Result{
		InputCount:  len(scoredVariants),
		OutputCount: len(selected),
		Duration:    time.Since(startTime),
	}

	if len(selected) > 0 {
		result.MaxScore = selected[0].Score
		result.MinScore = selected[len(selected)-1].Score
	}

	fmt.Printf("\n✅ Step 7 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   输出变异: %d\n", result.OutputCount)
	fmt.Printf("   分数范围: %d - %d\n", result.MaxScore, result.MinScore)

	return result, nil
}

// calculateScore 计算变异得分（与 Python 版本一致）
// 评分规则：
// - AlphaMissense_score ≥ 0.564: +2分
// - ClinPred_score ≥ 0.5: +2分
// - MetaRNN_score ≥ 0.5: +2分
// - REVEL_score ≥ 0.5: +2分
// - CADD_phred ≥ 20: +1分
// - PrimateAI_score ≥ 0.7: +1分
// - ESM1b_score < -3: +1分
func calculateScore(v *vcf.Variant) int {
	score := 0

	// AlphaMissense ≥ 0.564: +2分
	if v.GetInfoFloatMax("AlphaMissense_score") >= 0.564 {
		score += 2
	}

	// ClinPred ≥ 0.5: +2分
	if v.GetInfoFloatMax("ClinPred_score") >= 0.5 {
		score += 2
	}

	// MetaRNN ≥ 0.5: +2分
	if v.GetInfoFloatMax("MetaRNN_score") >= 0.5 {
		score += 2
	}

	// REVEL ≥ 0.5: +2分
	if v.GetInfoFloatMax("REVEL_score") >= 0.5 {
		score += 2
	}

	// CADD ≥ 20: +1分
	if v.GetInfoFloatMax("CADD_phred") >= 20.0 {
		score += 1
	}

	// PrimateAI ≥ 0.7: +1分
	if v.GetInfoFloatMax("PrimateAI_score") >= 0.7 {
		score += 1
	}

	// ESM1b < -3: +1分（取最小值）
	if v.GetInfoFloatMin("ESM1b_score") < -3.0 {
		score += 1
	}

	return score
}
