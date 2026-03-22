package pipeline

import (
	"fmt"
	"io"
	"strconv"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step3Result Step3 的运行结果
type Step3Result struct {
	InputCount   int
	OutputCount  int
	FailedQUAL   int
	FailedFilter int
	FailedDP     int
	FailedGQ     int
	Duration     time.Duration
}

// Step3QualityFilter 执行 Step3: 测序质量过滤
func Step3QualityFilter(cfg *config.Config) (*Step3Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 3: 测序质量过滤")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	fmt.Printf("\n  筛选条件:\n")
	fmt.Printf("    - QUAL >= %.0f\n", cfg.MinQUAL)
	fmt.Printf("    - FILTER = PASS\n")
	fmt.Printf("    - DP >= %d\n", cfg.MinDP)
	fmt.Printf("    - GQ >= %d\n", cfg.MinGQ)

	// 打开输入文件
	parser, err := vcf.NewParser(cfg.VCF2)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

	// 创建输出文件
	writer, err := vcf.NewWriter(cfg.VCF3)
	if err != nil {
		return nil, err
	}
	defer writer.Close()

	// 写入 header
	writer.WriteHeader(parser.Header)

	result := &Step3Result{}

	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		result.InputCount++

		// 检查 QUAL
		if v.Qual < cfg.MinQUAL {
			result.FailedQUAL++
			continue
		}

		// 检查 FILTER
		if v.Filter != "PASS" && v.Filter != "." {
			result.FailedFilter++
			continue
		}

		// 提取 DP 和 GQ（从 FORMAT/样本字段）
		dp, gq := extractDPGQ(v)

		if dp < cfg.MinDP {
			result.FailedDP++
			continue
		}

		if gq < cfg.MinGQ {
			result.FailedGQ++
			continue
		}

		// 通过所有过滤
		result.OutputCount++
		writer.WriteVariant(v)
	}

	result.Duration = time.Since(startTime)

	fmt.Printf("\n✅ Step 3 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   输出变异: %d (%.2f%%)\n", result.OutputCount,
		float64(result.OutputCount)/float64(result.InputCount)*100)
	fmt.Printf("   失败原因:\n")
	fmt.Printf("     - QUAL < %.0f: %d\n", cfg.MinQUAL, result.FailedQUAL)
	fmt.Printf("     - FILTER ≠ PASS: %d\n", result.FailedFilter)
	fmt.Printf("     - DP < %d: %d\n", cfg.MinDP, result.FailedDP)
	fmt.Printf("     - GQ < %d: %d\n", cfg.MinGQ, result.FailedGQ)

	return result, nil
}

// extractDPGQ 从 FORMAT/样本字段提取 DP 和 GQ
func extractDPGQ(v *vcf.Variant) (dp int, gq int) {
	if v.Format == "" || len(v.Samples) == 0 {
		return 0, 0
	}

	formatFields := strings.Split(v.Format, ":")
	sampleFields := strings.Split(v.Samples[0], ":")

	dpIdx := -1
	gqIdx := -1
	for i, f := range formatFields {
		switch f {
		case "DP":
			dpIdx = i
		case "GQ":
			gqIdx = i
		}
	}

	if dpIdx >= 0 && dpIdx < len(sampleFields) {
		dp, _ = strconv.Atoi(sampleFields[dpIdx])
	}
	if gqIdx >= 0 && gqIdx < len(sampleFields) {
		gq, _ = strconv.Atoi(sampleFields[gqIdx])
	}

	return dp, gq
}
