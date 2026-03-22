package pipeline

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step1Result Step1 的运行结果
type Step1Result struct {
	InputCount   int64
	OutputCount  int64
	Duration     time.Duration
	VariantsRate float64
}

// Step1GnomADFilter 执行 Step1: gnomAD 注释和频率过滤
// 这是整个流水线最耗时的步骤，使用并发处理优化
func Step1GnomADFilter(cfg *config.Config) (*Step1Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 1: gnomAD 注释与频率过滤 (Go 并发版)")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 1. 首先使用 vcfanno 进行注释
	fmt.Println("\n[1/2] 使用 vcfanno 注释 gnomAD 频率...")
	annotatedVCF := cfg.VCF1 + ".annotated.tmp"
	if err := runVcfannoGnomAD(cfg, annotatedVCF); err != nil {
		return nil, fmt.Errorf("vcfanno 注释失败: %w", err)
	}
	vcfannoTime := time.Since(startTime)
	fmt.Printf("  vcfanno 完成，耗时: %.2f秒\n", vcfannoTime.Seconds())

	// 2. 并发过滤
	fmt.Println("\n[2/2] 并发过滤低频变异...")
	filterStart := time.Now()
	result, err := parallelFilter(annotatedVCF, cfg)
	if err != nil {
		return nil, fmt.Errorf("过滤失败: %w", err)
	}
	filterTime := time.Since(filterStart)

	// 清理临时文件
	os.Remove(annotatedVCF)

	result.Duration = time.Since(startTime)
	result.VariantsRate = float64(result.InputCount) / result.Duration.Seconds()

	// 打印结果
	fmt.Println("\n" + strings.Repeat("-", 60))
	fmt.Printf("✅ Step 1 完成！\n")
	fmt.Printf("   总耗时: %.2f秒 (vcfanno: %.2f秒, 过滤: %.2f秒)\n",
		result.Duration.Seconds(), vcfannoTime.Seconds(), filterTime.Seconds())
	fmt.Printf("   处理速度: %.0f 变异/秒\n", result.VariantsRate)
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   输出变异: %d (%.2f%%)\n", result.OutputCount,
		float64(result.OutputCount)/float64(result.InputCount)*100)
	fmt.Printf("   输出文件: %s\n", cfg.VCF1)

	return result, nil
}

// runVcfannoGnomAD 运行 vcfanno 进行 gnomAD 注释
func runVcfannoGnomAD(cfg *config.Config, outputFile string) error {
	// 创建配置文件
	// 列号: 27=gnomAD4.1_joint_POPMAX_AF, 26=gnomAD4.1_joint_flag（实际不是AF，是flag）
	// 根据实际 dbNSFP 结构调整
	tomlContent := fmt.Sprintf(`[[annotation]]
file="%s"
columns=[27]
names=["gnomAD4.1_joint_POPMAX_AF"]
ops=["self"]
`, cfg.DBNSFPFile)

	tomlPath := cfg.VCF1 + ".toml"
	if err := os.WriteFile(tomlPath, []byte(tomlContent), 0644); err != nil {
		return err
	}
	defer os.Remove(tomlPath)

	// 运行 vcfanno
	cmd := exec.Command("vcfanno", "-p", fmt.Sprintf("%d", cfg.NumWorkers), tomlPath, cfg.InputVCF)
	outFile, err := os.Create(outputFile)
	if err != nil {
		return err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr

	return cmd.Run()
}

// parallelFilter 并发过滤 VCF 文件（保持顺序版）
func parallelFilter(inputFile string, cfg *config.Config) (*Step1Result, error) {
	// 打开输入文件
	file, err := os.Open(inputFile)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := bufio.NewReaderSize(file, 4*1024*1024) // 4MB 缓冲

	// 读取 header
	var sampleLine string
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			break
		}
		if strings.HasPrefix(line, "##") {
			continue
		}
		if strings.HasPrefix(line, "#CHROM") {
			sampleLine = strings.TrimRight(line, "\r\n")
			break
		}
	}

	// 顺序读取所有变异
	var allLines []string
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		line = strings.TrimRight(line, "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		allLines = append(allLines, line)
	}

	inputCount := int64(len(allLines))
	fmt.Printf("  读取 %d 个变异\n", inputCount)

	// 并发处理，但保持结果顺序
	numWorkers := runtime.NumCPU()
	chunkSize := (len(allLines) + numWorkers - 1) / numWorkers

	type indexedResult struct {
		index   int
		variant *vcf.Variant
		passed  bool
	}

	results := make([]indexedResult, len(allLines))
	var wg sync.WaitGroup
	var processedCount int64

	// 分块并发处理
	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if end > len(allLines) {
			end = len(allLines)
		}
		if start >= len(allLines) {
			break
		}

		wg.Add(1)
		go func(startIdx, endIdx int) {
			defer wg.Done()
			for i := startIdx; i < endIdx; i++ {
				line := allLines[i]
				v, err := vcf.ParseLine(line)
				if err != nil {
					results[i] = indexedResult{index: i, passed: false}
					atomic.AddInt64(&processedCount, 1)
					continue
				}

				// 检查 gnomAD 频率
				// 逻辑: 必须有AF且 < 0.01 才保留
				af, hasAF := v.GetInfoFloat("gnomAD4.1_joint_POPMAX_AF")
				if !hasAF || af >= cfg.GnomADMaxAF {
					// 没有AF或频率太高，过滤掉
					results[i] = indexedResult{index: i, passed: false}
				} else {
					// 有AF且频率低于阈值，保留
					results[i] = indexedResult{index: i, variant: v, passed: true}
				}
				
				count := atomic.AddInt64(&processedCount, 1)
				if count%100000 == 0 {
					fmt.Printf("  已处理 %d 个变异\r", count)
				}
			}
		}(start, end)
	}

	wg.Wait()
	fmt.Printf("  已处理 %d 个变异\n", processedCount)

	// 创建输出文件
	outFile, err := os.Create(cfg.VCF1)
	if err != nil {
		return nil, err
	}
	defer outFile.Close()

	writer := bufio.NewWriterSize(outFile, 4*1024*1024)
	defer writer.Flush()

	// 写入 header
	cleanHeader := vcf.CreateCleanHeader(sampleLine)
	cleanHeader.Write(writer)

	// 按原始顺序写入通过的变异
	var outputCount int64
	for _, r := range results {
		if r.passed && r.variant != nil {
			writer.WriteString(r.variant.ToLine() + "\n")
			outputCount++
		}
	}

	return &Step1Result{
		InputCount:  inputCount,
		OutputCount: outputCount,
	}, nil
}

// Step1Direct 直接处理（不使用 vcfanno，使用 tabix 按需查询）
// 这个版本更快，但需要 dbNSFP 有 tabix 索引
func Step1Direct(cfg *config.Config) (*Step1Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 1: gnomAD 注释与频率过滤 (Go 直接查询版)")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 检查 tabix 索引
	if _, err := os.Stat(cfg.DBNSFPFile + ".tbi"); os.IsNotExist(err) {
		fmt.Println("⚠️ dbNSFP 索引不存在，回退到 vcfanno 模式")
		return Step1GnomADFilter(cfg)
	}

	// 打开 dbNSFP 获取列索引
	dbFile, err := os.Open(cfg.DBNSFPFile)
	if err != nil {
		return nil, err
	}
	gzReader, err := gzip.NewReader(dbFile)
	if err != nil {
		dbFile.Close()
		return nil, err
	}

	headerReader := bufio.NewReader(gzReader)
	headerLine, err := headerReader.ReadString('\n')
	if err != nil {
		gzReader.Close()
		dbFile.Close()
		return nil, err
	}
	gzReader.Close()
	dbFile.Close()

	// 解析 header 获取列索引
	headers := strings.Split(strings.TrimRight(headerLine, "\r\n"), "\t")
	colIndex := make(map[string]int)
	for i, h := range headers {
		colIndex[strings.TrimPrefix(h, "#")] = i
	}

	gnomadCol := colIndex["gnomAD_genomes_POPMAX_AF"]
	gnomadJointCol := colIndex["gnomAD_genomes_AF"]
	if gnomadCol == 0 {
		gnomadCol = 24 // 默认列
		gnomadJointCol = 23
	}

	// 打开输入 VCF
	vcfParser, err := vcf.NewParser(cfg.InputVCF)
	if err != nil {
		return nil, err
	}
	defer vcfParser.Close()

	// 创建输出
	outWriter, err := vcf.NewWriter(cfg.VCF1)
	if err != nil {
		return nil, err
	}
	defer outWriter.Close()

	// 写入 header
	cleanHeader := vcf.CreateCleanHeader(vcfParser.Header.SampleLine)
	outWriter.WriteHeader(cleanHeader)

	var inputCount, outputCount int64
	lastReport := time.Now()

	for {
		v, err := vcfParser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		inputCount++

		// 使用 tabix 查询
		chromNum := v.Chrom
		if strings.HasPrefix(chromNum, "chr") {
			chromNum = chromNum[3:]
		}

		region := fmt.Sprintf("%s:%d-%d", chromNum, v.Pos, v.Pos)
		cmd := exec.Command("tabix", cfg.DBNSFPFile, region)
		output, _ := cmd.Output()

		if len(output) == 0 {
			continue // 没有 gnomAD 数据
		}

		// 解析 dbNSFP 输出
		lines := strings.Split(strings.TrimSpace(string(output)), "\n")
		var gnomadAF float64
		var found bool

		for _, line := range lines {
			fields := strings.Split(line, "\t")
			if len(fields) <= gnomadCol {
				continue
			}

			// 检查 ref/alt 匹配
			refCol := colIndex["ref"]
			altCol := colIndex["alt"]
			if refCol < len(fields) && altCol < len(fields) {
				if fields[refCol] != v.Ref || fields[altCol] != v.Alt {
					continue
				}
			}

			// 获取 gnomAD 频率
			afStr := fields[gnomadCol]
			if afStr != "" && afStr != "." {
				af, err := strconv.ParseFloat(afStr, 64)
				if err == nil {
					gnomadAF = af
					found = true

					// 也获取 joint AF
					if gnomadJointCol < len(fields) {
						jointAF := fields[gnomadJointCol]
						if jointAF != "" && jointAF != "." {
							v.SetInfo("gnomAD4.1_joint", jointAF)
						}
					}
					break
				}
			}
		}

		if !found {
			continue
		}

		// 过滤高频变异 (与 Python 一致: af < 0.01)
		if gnomadAF >= cfg.GnomADMaxAF {
			continue
		}

		v.SetInfo("gnomAD4.1_joint_POPMAX_AF", fmt.Sprintf("%.6e", gnomadAF))
		outputCount++
		outWriter.WriteVariant(v)

		// 进度报告
		if time.Since(lastReport) > time.Second {
			fmt.Printf("  已处理 %d 个变异，通过 %d 个\r", inputCount, outputCount)
			lastReport = time.Now()
		}
	}

	fmt.Println()

	duration := time.Since(startTime)
	result := &Step1Result{
		InputCount:   inputCount,
		OutputCount:  outputCount,
		Duration:     duration,
		VariantsRate: float64(inputCount) / duration.Seconds(),
	}

	fmt.Println("\n" + strings.Repeat("-", 60))
	fmt.Printf("✅ Step 1 完成！\n")
	fmt.Printf("   总耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   处理速度: %.0f 变异/秒\n", result.VariantsRate)
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   输出变异: %d (%.4f%%)\n", result.OutputCount,
		float64(result.OutputCount)/float64(result.InputCount)*100)

	return result, nil
}
