package pipeline

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"sort"
	"strconv"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step8Result Step8 的运行结果
type Step8Result struct {
	ClinVarCount int
	DbscSNVCount int
	TopCount     int
	TotalCount   int // 去重后总数
	Duration     time.Duration
}

// Step8Integrate 执行 Step8: 整合三个决赛圈文件
func Step8Integrate(cfg *config.Config) (*Step8Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 8: 整合三个决赛圈文件")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 1. 合并三个决赛圈文件（去重）
	fmt.Println("\n  [1/3] 合并决赛圈文件...")
	allVariants := make(map[string]*vcf.Variant)
	var sampleLine string

	// 读取 ClinVar 决赛圈
	clinvarVariants, sl, _ := readVCFToMap(cfg.FinalClinVar)
	if sl != "" {
		sampleLine = sl
	}
	for k, v := range clinvarVariants {
		allVariants[k] = v
	}
	clinvarCount := len(clinvarVariants)
	fmt.Printf("    final_clinvar.vcf: %d 个变异\n", clinvarCount)

	// 读取 dbscSNV 决赛圈
	dbscsnvVariants, sl, _ := readVCFToMap(cfg.FinalDbscSNV)
	if sl != "" && sampleLine == "" {
		sampleLine = sl
	}
	for k, v := range dbscsnvVariants {
		if _, exists := allVariants[k]; !exists {
			allVariants[k] = v
		} else {
			// 合并 INFO
			mergeInfo(allVariants[k], v)
		}
	}
	dbscsnvCount := len(dbscsnvVariants)
	fmt.Printf("    final_dbscSNV.vcf: %d 个变异\n", dbscsnvCount)

	// 读取 Top 决赛圈
	topVariants, sl, _ := readVCFToMap(cfg.FinalTop)
	if sl != "" && sampleLine == "" {
		sampleLine = sl
	}
	for k, v := range topVariants {
		if _, exists := allVariants[k]; !exists {
			allVariants[k] = v
		} else {
			mergeInfo(allVariants[k], v)
		}
	}
	topCount := len(topVariants)
	fmt.Printf("    final_top.vcf: %d 个变异\n", topCount)
	fmt.Printf("    合并后去重: %d 个变异\n", len(allVariants))

	// 2. 补全注释
	fmt.Println("\n  [2/3] 补全缺失注释...")
	if err := annotateVariants(allVariants, cfg); err != nil {
		fmt.Printf("  警告: 补全注释时出错: %v\n", err)
	}
	
	// 确保所有必需字段都存在
	for _, v := range allVariants {
		ensureAllFields(v)
	}

	// 3. 排序并写入
	fmt.Println("\n  [3/3] 写入最终结果...")
	sortedVariants := sortVariants(allVariants)

	writer, err := vcf.NewWriter(cfg.FinalIntegrated)
	if err != nil {
		return nil, err
	}
	defer writer.Close()

	// 写入干净的 header
	cleanHeader := vcf.CreateCleanHeader(sampleLine)
	writer.WriteHeader(cleanHeader)

	for _, v := range sortedVariants {
		writer.WriteVariant(v)
	}

	result := &Step8Result{
		ClinVarCount: clinvarCount,
		DbscSNVCount: dbscsnvCount,
		TopCount:     topCount,
		TotalCount:   len(sortedVariants),
		Duration:     time.Since(startTime),
	}

	fmt.Printf("\n✅ Step 8 完成！耗时: %.2f秒\n", result.Duration.Seconds())
	fmt.Printf("   ClinVar 致病变异: %d\n", result.ClinVarCount)
	fmt.Printf("   dbscSNV 高分变异: %d\n", result.DbscSNVCount)
	fmt.Printf("   Top 评分变异: %d\n", result.TopCount)
	fmt.Printf("   最终整合结果: %d\n", result.TotalCount)

	return result, nil
}

// readVCFToMap 读取 VCF 文件到 map
func readVCFToMap(filename string) (map[string]*vcf.Variant, string, error) {
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		return make(map[string]*vcf.Variant), "", nil
	}

	parser, err := vcf.NewParser(filename)
	if err != nil {
		return nil, "", err
	}
	defer parser.Close()

	variants := make(map[string]*vcf.Variant)
	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}
		variants[v.Key()] = v
	}

	return variants, parser.Header.SampleLine, nil
}

// mergeInfo 合并两个变异的 INFO 字段
func mergeInfo(target, source *vcf.Variant) {
	for k, v := range source.Info {
		if _, exists := target.Info[k]; !exists || target.Info[k] == "." {
			target.Info[k] = v
		}
	}
}

// annotateVariants 补全变异注释
func annotateVariants(variants map[string]*vcf.Variant, cfg *config.Config) error {
	// 收集需要查询的位置
	positions := make([]struct {
		Chrom string
		Pos   int
		Ref   string
		Alt   string
	}, 0, len(variants))

	for _, v := range variants {
		positions = append(positions, struct {
			Chrom string
			Pos   int
			Ref   string
			Alt   string
		}{v.Chrom, v.Pos, v.Ref, v.Alt})
	}

	// 查询 dbscSNV
	dbscsnvData := queryDbscSNVBatch(cfg.DBscSNVFile, positions)
	matched := 0
	for key, data := range dbscsnvData {
		if v, ok := variants[key]; ok {
			if v.GetInfoString("DBSCSNV_ADA") == "" {
				v.SetInfo("DBSCSNV_ADA", fmt.Sprintf("%.6f", data.AdaScore))
				v.SetInfo("DBSCSNV_RF", fmt.Sprintf("%.6f", data.RFScore))
				matched++
			}
		}
	}
	fmt.Printf("    dbscSNV 补充匹配: %d 个\n", matched)

	// 使用 tabix 查询 dbNSFP 补全评分
	dbnsfpMatched := 0
	for _, v := range variants {
		if annotateFromDbNSFP(v, cfg.DBNSFPFile) {
			dbnsfpMatched++
		}
	}
	fmt.Printf("    dbNSFP 补充匹配: %d 个\n", dbnsfpMatched)

	return nil
}

// queryDbscSNVBatch 批量查询 dbscSNV
func queryDbscSNVBatch(dbPath string, positions []struct {
	Chrom string
	Pos   int
	Ref   string
	Alt   string
}) map[string]*DbscSNVRecord {
	// 构建查询集合
	posSet := make(map[string]bool)
	for _, p := range positions {
		key := fmt.Sprintf("%s_%d_%s_%s", stripChr(p.Chrom), p.Pos, p.Ref, p.Alt)
		posSet[key] = true
	}

	file, err := os.Open(dbPath)
	if err != nil {
		return nil
	}
	defer file.Close()

	reader := bufio.NewReaderSize(file, 4*1024*1024)
	headerLine, _ := reader.ReadString('\n')
	headers := strings.Split(strings.TrimRight(headerLine, "\r\n"), "\t")

	colIndex := make(map[string]int)
	for i, h := range headers {
		colIndex[h] = i
	}

	result := make(map[string]*DbscSNVRecord)

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			break
		}

		fields := strings.Split(strings.TrimRight(line, "\r\n"), "\t")
		if len(fields) < 6 {
			continue
		}

		chrom := fields[colIndex["chr"]]
		pos := fields[colIndex["pos"]]
		ref := fields[colIndex["ref"]]
		alt := fields[colIndex["alt"]]

		key := fmt.Sprintf("%s_%s_%s_%s", chrom, pos, ref, alt)
		// 也尝试带 chr 前缀的 key
		keyWithChr := fmt.Sprintf("chr%s_%s_%s_%s", chrom, pos, ref, alt)

		if posSet[key] || posSet[keyWithChr] {
			ada, _ := strconv.ParseFloat(fields[colIndex["ada_score"]], 64)
			rf, _ := strconv.ParseFloat(fields[colIndex["rf_score"]], 64)
			posInt, _ := strconv.Atoi(pos)

			record := &DbscSNVRecord{
				Chrom:    chrom,
				Pos:      posInt,
				Ref:      ref,
				Alt:      alt,
				AdaScore: ada,
				RFScore:  rf,
			}

			result[key] = record
			result[keyWithChr] = record
		}
	}

	return result
}

// annotateFromDbNSFP 从 dbNSFP 补全注释
func annotateFromDbNSFP(v *vcf.Variant, dbPath string) bool {
	chromNum := v.Chrom
	if strings.HasPrefix(chromNum, "chr") {
		chromNum = chromNum[3:]
	}

	region := fmt.Sprintf("%s:%d-%d", chromNum, v.Pos, v.Pos)
	cmd := exec.Command("tabix", dbPath, region)
	output, err := cmd.Output()
	if err != nil || len(output) == 0 {
		return false
	}

	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for _, line := range lines {
		fields := strings.Split(line, "\t")
		if len(fields) < 50 {
			continue
		}

		// 检查 ref/alt
		if fields[2] != v.Ref || fields[3] != v.Alt {
			continue
		}

		// 补全字段
		updated := false
		fieldMap := map[int]string{
			6:  "aaref",
			7:  "aaalt",
			8:  "rs_dbSNP",
			9:  "genename",
			35: "MetaRNN_score",
			36: "REVEL_score",
			37: "PrimateAI_score",
			38: "ClinPred_score",
			39: "ESM1b_score",
			40: "AlphaMissense_score",
			41: "CADD_phred",
		}

		for col, name := range fieldMap {
			if col < len(fields) && v.GetInfoString(name) == "" {
				val := fields[col]
				if val != "" && val != "." {
					v.SetInfo(name, val)
					updated = true
				}
			}
		}

		return updated
	}

	return false
}

// sortVariants 按染色体和位置排序变异
func sortVariants(variants map[string]*vcf.Variant) []*vcf.Variant {
	result := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		result = append(result, v)
	}

	sort.Slice(result, func(i, j int) bool {
		if result[i].ChromNum() != result[j].ChromNum() {
			return result[i].ChromNum() < result[j].ChromNum()
		}
		return result[i].Pos < result[j].Pos
	})

	return result
}

// ensureAllFields 确保变异包含所有必需字段，缺失的设为 "."
// 同时清理多值字段（评分取最大/最小值，其他字段去重）
func ensureAllFields(v *vcf.Variant) {
	requiredFields := []string{
		"gnomAD4.1_joint_POPMAX_AF", "gnomAD4.1_joint",
		"DBSCSNV_ADA", "DBSCSNV_RF",
		"MetaRNN_score", "REVEL_score", "PrimateAI_score", "ClinPred_score",
		"ESM1b_score", "AlphaMissense_score", "CADD_phred",
		"aaref", "aaalt", "rs_dbSNP", "genename",
		"clinvar_id", "clinvar_clnsig", "clinvar_trait", "clinvar_review",
		"clinvar_hgvs", "clinvar_var_source", "clinvar_MedGen_id",
		"clinvar_OMIM_id", "clinvar_Orphanet_id",
	}

	for _, field := range requiredFields {
		val := v.GetInfoString(field)
		if val == "" {
			v.SetInfo(field, ".")
		} else {
			// 清理多值字段
			cleanedVal := vcf.CleanInfoValue(field, val)
			v.SetInfo(field, cleanedVal)
		}
	}
}
