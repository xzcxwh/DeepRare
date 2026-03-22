package pipeline

import (
	"fmt"
	"io"
	"os"
	"os/exec"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step6Result Step6 的运行结果
type Step6Result struct {
	InputCount  int
	OutputCount int
	Duration    time.Duration
}

// Step6ScoreFilter 执行 Step6: dbNSFP 评分注释与过滤
func Step6ScoreFilter(cfg *config.Config) (*Step6Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 6: dbNSFP 评分注释与过滤")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 1. 使用 vcfanno 添加评分注释
	fmt.Println("\n  [1/2] 使用 vcfanno 添加评分注释...")
	annotatedVCF := cfg.VCF6 + ".annotated.tmp"
	if err := runVcfannoScores(cfg, annotatedVCF); err != nil {
		return nil, fmt.Errorf("vcfanno 注释失败: %w", err)
	}

	// 2. 根据评分过滤
	fmt.Println("\n  [2/2] 根据评分过滤变异...")
	result, err := filterByScores(annotatedVCF, cfg.VCF6, cfg)
	if err != nil {
		return nil, err
	}

	// 清理临时文件
	os.Remove(annotatedVCF)

	result.Duration = time.Since(startTime)

	fmt.Printf("\n✅ Step 6 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   输出变异: %d (%.2f%%)\n", result.OutputCount,
		float64(result.OutputCount)/float64(result.InputCount)*100)

	return result, nil
}

// runVcfannoScores 运行 vcfanno 添加评分
func runVcfannoScores(cfg *config.Config, outputFile string) error {
	// 评分字段列号（根据 dbNSFP 实际结构）
	// 列号: 5=aaref, 6=aaalt, 7=rs_dbSNP, 9=genename
	// 19=MetaRNN, 20=REVEL, 21=PrimateAI, 22=ClinPred, 23=ESM1b, 24=AlphaMissense, 25=CADD
	tomlContent := fmt.Sprintf(`[[annotation]]
file="%s"
columns=[19, 20, 21, 22, 23, 24, 25, 5, 6, 7, 9]
names=["MetaRNN_score", "REVEL_score", "PrimateAI_score", "ClinPred_score", "ESM1b_score", "AlphaMissense_score", "CADD_phred", "aaref", "aaalt", "rs_dbSNP", "genename"]
ops=["self", "self", "self", "self", "self", "self", "self", "self", "self", "self", "self"]
`, cfg.DBNSFPFile)

	tomlPath := cfg.VCF6 + ".toml"
	if err := os.WriteFile(tomlPath, []byte(tomlContent), 0644); err != nil {
		return err
	}
	defer os.Remove(tomlPath)

	cmd := exec.Command("vcfanno", tomlPath, cfg.VCF5)
	outFile, err := os.Create(outputFile)
	if err != nil {
		return err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr

	return cmd.Run()
}

// filterByScores 根据评分过滤变异
func filterByScores(inputFile, outputFile string, cfg *config.Config) (*Step6Result, error) {
	parser, err := vcf.NewParser(inputFile)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

	writer, err := vcf.NewWriter(outputFile)
	if err != nil {
		return nil, err
	}
	defer writer.Close()

	writer.WriteHeader(parser.Header)

	result := &Step6Result{}

	// 评分阈值（与 Python 版本一致）
	thresholds := map[string]float64{
		"AlphaMissense_score": 0.564, // >= 0.564
		"REVEL_score":         0.5,   // >= 0.5
		"MetaRNN_score":       0.5,   // >= 0.5
		"PrimateAI_score":     0.7,   // >= 0.7 (Python 版本)
		"CADD_phred":          20.0,  // >= 20
		"ClinPred_score":      0.5,   // >= 0.5
	}

	// ESM1b_score 是越小越致病，阈值 < -3 (Python 版本)
	esm1bThreshold := -3.0

	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		result.InputCount++

		// 统计满足条件的评分数量
		passCount := 0

		for field, threshold := range thresholds {
			// 使用 max 值判断（处理多值字段）
			maxScore := v.GetInfoFloatMax(field)
			if maxScore >= threshold {
				passCount++
			}
		}

		// ESM1b_score 特殊处理：取 min 值，< -3 为致病
		esm1bMin := v.GetInfoFloatMin("ESM1b_score")
		if esm1bMin < esm1bThreshold {
			passCount++
		}

		// 至少满足 N 个条件
		if passCount >= cfg.MinScoreCount {
			result.OutputCount++
			writer.WriteVariant(v)
		}
	}

	return result, nil
}
