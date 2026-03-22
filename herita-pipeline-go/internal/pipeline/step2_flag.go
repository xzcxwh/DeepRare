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

// Step2Result Step2 的运行结果
type Step2Result struct {
	InputCount  int
	OutputCount int
	Duration    time.Duration
}

// Step2GnomADFlag 执行 Step2: gnomAD Flag 注释
func Step2GnomADFlag(cfg *config.Config) (*Step2Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 2: gnomAD Flag 注释")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 创建 vcfanno 配置
	tomlContent := fmt.Sprintf(`[[annotation]]
file="%s"
columns=[26]
names=["gnomAD4.1_joint_flag"]
ops=["self"]
`, cfg.DBNSFPFile)

	tomlPath := cfg.VCF2 + ".toml"
	if err := os.WriteFile(tomlPath, []byte(tomlContent), 0644); err != nil {
		return nil, err
	}
	defer os.Remove(tomlPath)

	// 运行 vcfanno
	fmt.Println("\n  运行 vcfanno 添加 gnomAD flag...")
	cmd := exec.Command("vcfanno", tomlPath, cfg.VCF1)
	outFile, err := os.Create(cfg.VCF2)
	if err != nil {
		return nil, err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr

	if err := cmd.Run(); err != nil {
		return nil, fmt.Errorf("vcfanno 执行失败: %w", err)
	}

	// 统计变异数量
	inputCount := countVariants(cfg.VCF1)
	outputCount := countVariants(cfg.VCF2)

	result := &Step2Result{
		InputCount:  inputCount,
		OutputCount: outputCount,
		Duration:    time.Since(startTime),
	}

	fmt.Printf("\n✅ Step 2 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   输入变异: %d, 输出变异: %d\n", inputCount, outputCount)

	return result, nil
}

// countVariants 统计 VCF 文件中的变异数量
func countVariants(filename string) int {
	parser, err := vcf.NewParser(filename)
	if err != nil {
		return 0
	}
	defer parser.Close()

	count := 0
	for {
		_, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}
		count++
	}
	return count
}
