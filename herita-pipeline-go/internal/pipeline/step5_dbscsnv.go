package pipeline

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step5Result Step5 的运行结果
type Step5Result struct {
	InputCount   int
	OutputCount  int
	DbscSNVCount int // 进入决赛圈的 dbscSNV 变异
	MatchedCount int // 匹配到 dbscSNV 的变异数
	Duration     time.Duration
}

// DbscSNVRecord dbscSNV 记录
type DbscSNVRecord struct {
	Chrom    string
	Pos      int
	Ref      string
	Alt      string
	AdaScore float64
	RFScore  float64
}

// Step5DbscSNVAnnotation 执行 Step5: dbscSNV 注释
func Step5DbscSNVAnnotation(cfg *config.Config) (*Step5Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 5: dbscSNV 注释")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 1. 首先读取所有变异位置
	fmt.Println("\n  [1/3] 读取变异位置...")
	positions, variants, header, err := loadVariantPositions(cfg.VCF4)
	if err != nil {
		return nil, err
	}
	fmt.Printf("  需要查询 %d 个位点\n", len(positions))

	// 2. 从 dbscSNV 文件查找匹配
	fmt.Println("\n  [2/3] 查询 dbscSNV 数据库...")
	dbscsnvData, err := queryDbscSNV(cfg.DBscSNVFile, positions)
	if err != nil {
		return nil, err
	}
	fmt.Printf("  匹配到 %d 个位点\n", len(dbscsnvData))

	// 3. 写入结果
	fmt.Println("\n  [3/3] 写入结果...")

	// 创建输出文件
	writer5, err := vcf.NewWriter(cfg.VCF5)
	if err != nil {
		return nil, err
	}
	defer writer5.Close()

	writerDbscSNV, err := vcf.NewWriter(cfg.FinalDbscSNV)
	if err != nil {
		return nil, err
	}
	defer writerDbscSNV.Close()

	// 写入 header
	writer5.WriteHeader(header)
	writerDbscSNV.WriteHeader(header)

	result := &Step5Result{
		InputCount:   len(variants),
		MatchedCount: len(dbscsnvData),
	}

	for _, v := range variants {
		key := fmt.Sprintf("%s_%d_%s_%s", stripChr(v.Chrom), v.Pos, v.Ref, v.Alt)

		if data, ok := dbscsnvData[key]; ok {
			// 添加注释
			v.SetInfo("DBSCSNV_ADA", fmt.Sprintf("%.6f", data.AdaScore))
			v.SetInfo("DBSCSNV_RF", fmt.Sprintf("%.6f", data.RFScore))

			// 检查是否进入决赛圈
			if data.AdaScore > cfg.DbscSNVThreshold && data.RFScore > cfg.DbscSNVThreshold {
				result.DbscSNVCount++
				writerDbscSNV.WriteVariant(v)
				fmt.Printf("  变异 %s:%d %s>%s 进入决赛圈: ada=%.3f, rf=%.3f\n",
					v.Chrom, v.Pos, v.Ref, v.Alt, data.AdaScore, data.RFScore)
				continue
			}
		}

		result.OutputCount++
		writer5.WriteVariant(v)
	}

	result.Duration = time.Since(startTime)

	fmt.Printf("\n✅ Step 5 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   输入变异: %d\n", result.InputCount)
	fmt.Printf("   dbscSNV 匹配: %d (%.1f%%)\n", result.MatchedCount,
		float64(result.MatchedCount)/float64(result.InputCount)*100)
	fmt.Printf("   进入决赛圈 (final_dbscSNV.vcf): %d\n", result.DbscSNVCount)
	fmt.Printf("   保留在 vcf5 中: %d\n", result.OutputCount)

	return result, nil
}

// loadVariantPositions 加载变异位置
func loadVariantPositions(filename string) (map[string]bool, []*vcf.Variant, *vcf.Header, error) {
	parser, err := vcf.NewParser(filename)
	if err != nil {
		return nil, nil, nil, err
	}
	defer parser.Close()

	positions := make(map[string]bool)
	var variants []*vcf.Variant

	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}

		key := fmt.Sprintf("%s_%d_%s_%s", stripChr(v.Chrom), v.Pos, v.Ref, v.Alt)
		positions[key] = true
		variants = append(variants, v)
	}

	return positions, variants, parser.Header, nil
}

// stripChr 移除 chr 前缀
func stripChr(chrom string) string {
	if strings.HasPrefix(chrom, "chr") {
		return chrom[3:]
	}
	return chrom
}

// queryDbscSNV 查询 dbscSNV 文件
func queryDbscSNV(dbPath string, positions map[string]bool) (map[string]*DbscSNVRecord, error) {
	file, err := os.Open(dbPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := bufio.NewReaderSize(file, 4*1024*1024)

	// 读取 header
	headerLine, _ := reader.ReadString('\n')
	headers := strings.Split(strings.TrimRight(headerLine, "\r\n"), "\t")

	// 找到列索引
	colIndex := make(map[string]int)
	for i, h := range headers {
		colIndex[h] = i
	}

	chrCol := colIndex["chr"]
	posCol := colIndex["pos"]
	refCol := colIndex["ref"]
	altCol := colIndex["alt"]
	adaCol := colIndex["ada_score"]
	rfCol := colIndex["rf_score"]

	result := make(map[string]*DbscSNVRecord)
	remaining := len(positions)

	for remaining > 0 {
		line, err := reader.ReadString('\n')
		if err != nil {
			break
		}

		fields := strings.Split(strings.TrimRight(line, "\r\n"), "\t")
		if len(fields) <= adaCol || len(fields) <= rfCol {
			continue
		}

		chrom := fields[chrCol]
		pos := fields[posCol]
		ref := fields[refCol]
		alt := fields[altCol]

		key := fmt.Sprintf("%s_%s_%s_%s", chrom, pos, ref, alt)

		if positions[key] {
			ada, _ := strconv.ParseFloat(fields[adaCol], 64)
			rf, _ := strconv.ParseFloat(fields[rfCol], 64)

			posInt, _ := strconv.Atoi(pos)
			result[key] = &DbscSNVRecord{
				Chrom:    chrom,
				Pos:      posInt,
				Ref:      ref,
				Alt:      alt,
				AdaScore: ada,
				RFScore:  rf,
			}

			remaining--
		}
	}

	return result, nil
}
