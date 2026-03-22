package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"regexp"
	"sort"
	"strings"
)

// Variant 表示一个变异
type Variant struct {
	Chrom   string
	Pos     string
	ID      string
	Ref     string
	Alt     string
	Qual    string
	Filter  string
	Info    string
	Format  string
	Samples []string
	Key     string // chr:pos:ref:alt 用于去重
}

// 有用的 INFO 字段（保留的）
var usefulInfoFields = map[string]bool{
	// 基因和功能注释
	"ANN":                true,
	"LOF":                true,
	"NMD":                true,
	"VEP_Gene":           true,
	"VEP_Consequence":    true,
	"VEP_AlphaMissense":  true,
	"VEP_AM_Pathogenicity": true,
	"VEP_CADD":           true,
	"VEP_REVEL":          true,
	"VEP_PhyloP":         true,
	"VEP_SpliceAI_Max":   true,
	"VEP_SpliceAI_Details": true,
	// ClinVar
	"ClinVar_CLNSIG":     true,
	"ClinVar_CLNDN":      true,
	"ClinVar_CLNREVSTAT": true,
	"ClinVar_ALLELEID":   true,
	"ClinVar_CLNDISDB":   true,
	// dbNSFP
	"dbNSFP_AlphaMissense": true,
	"dbNSFP_REVEL":         true,
	"dbNSFP_MetaRNN":       true,
	"dbNSFP_PrimateAI":     true,
	"dbNSFP_CADD":          true,
	"dbNSFP_ClinPred":      true,
	"dbNSFP_ESM1b":         true,
	// dbscSNV
	"dbscSNV_ada_score": true,
	"dbscSNV_rf_score":  true,
	// gnomAD
	"gnomAD_AF_grpmax": true,
	// 打分
	"HERITA_SCORE": true,
	// 来源标记
	"FINALIST_SOURCE": true,
}

func main() {
	// 命令行参数
	finalRound := flag.String("final-round", "", "当前决赛圈 VCF 文件")
	top50 := flag.String("top50", "", "TOP 50 VCF 文件")
	output := flag.String("output", "", "输出合并后的 VCF 文件")
	flag.Parse()

	if *finalRound == "" || *top50 == "" || *output == "" {
		fmt.Println("Usage: finalists_merger -final-round <vcf> -top50 <vcf> -output <vcf>")
		os.Exit(1)
	}

	fmt.Println("============================================================")
	fmt.Println("HERITA Finalists Merger - 决赛圈合并工具")
	fmt.Println("============================================================")
	fmt.Printf("决赛圈文件: %s\n", *finalRound)
	fmt.Printf("TOP 50 文件: %s\n", *top50)
	fmt.Printf("输出文件: %s\n", *output)
	fmt.Println("============================================================")

	// 读取两个文件
	variants1, header1 := readVCF(*finalRound, "PATHOGENIC_HIGH")
	variants2, _ := readVCF(*top50, "TOP50_SCORED")

	fmt.Printf("\n原决赛圈变异: %d\n", len(variants1))
	fmt.Printf("TOP 50 变异: %d\n", len(variants2))

	// 合并并去重
	merged := mergeVariants(variants1, variants2)
	fmt.Printf("合并后变异 (去重): %d\n", len(merged))

	// 清理注释
	for i := range merged {
		merged[i].Info = cleanInfo(merged[i].Info)
	}

	// 按位置排序
	sortVariants(merged)

	// 写入输出
	writeVCF(*output, merged, header1)

	fmt.Println("\n============================================================")
	fmt.Println("合并完成!")
	fmt.Println("============================================================")
}

// readVCF 读取 VCF 文件
func readVCF(filename string, source string) ([]Variant, []string) {
	file, err := os.Open(filename)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening %s: %v\n", filename, err)
		os.Exit(1)
	}
	defer file.Close()

	var variants []Variant
	var headers []string

	scanner := bufio.NewScanner(file)
	// 增加缓冲区大小
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()

		if strings.HasPrefix(line, "#") {
			headers = append(headers, line)
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 8 {
			continue
		}

		v := Variant{
			Chrom:  fields[0],
			Pos:    fields[1],
			ID:     fields[2],
			Ref:    fields[3],
			Alt:    fields[4],
			Qual:   fields[5],
			Filter: fields[6],
			Info:   fields[7],
			Key:    fmt.Sprintf("%s:%s:%s:%s", fields[0], fields[1], fields[3], fields[4]),
		}

		if len(fields) > 8 {
			v.Format = fields[8]
		}
		if len(fields) > 9 {
			v.Samples = fields[9:]
		}

		// 添加来源标记
		if !strings.Contains(v.Info, "FINALIST_SOURCE=") {
			v.Info = fmt.Sprintf("FINALIST_SOURCE=%s;%s", source, v.Info)
		}

		variants = append(variants, v)
	}

	return variants, headers
}

// mergeVariants 合并变异并去重
func mergeVariants(v1, v2 []Variant) []Variant {
	seen := make(map[string]int) // key -> index in result
	var result []Variant

	// 先添加原决赛圈的变异
	for _, v := range v1 {
		if _, exists := seen[v.Key]; !exists {
			seen[v.Key] = len(result)
			result = append(result, v)
		}
	}

	// 添加 TOP 50 中不重复的变异
	duplicates := 0
	for _, v := range v2 {
		if idx, exists := seen[v.Key]; exists {
			// 合并来源信息
			if !strings.Contains(result[idx].Info, "TOP50_SCORED") {
				result[idx].Info = strings.Replace(result[idx].Info, 
					"FINALIST_SOURCE=PATHOGENIC_HIGH", 
					"FINALIST_SOURCE=PATHOGENIC_HIGH,TOP50_SCORED", 1)
			}
			duplicates++
		} else {
			seen[v.Key] = len(result)
			result = append(result, v)
		}
	}

	if duplicates > 0 {
		fmt.Printf("重复变异 (已合并): %d\n", duplicates)
	}

	return result
}

// cleanInfo 清理 INFO 字段
func cleanInfo(info string) string {
	// 解析 INFO
	fields := strings.Split(info, ";")
	cleaned := make(map[string]string)
	fieldOrder := []string{}

	for _, field := range fields {
		if field == "" {
			continue
		}

		parts := strings.SplitN(field, "=", 2)
		key := parts[0]
		value := ""
		if len(parts) > 1 {
			value = parts[1]
		}

		// 只保留有用的字段
		if !usefulInfoFields[key] {
			continue
		}

		// 清理多值字段中的点号
		if strings.HasPrefix(key, "dbNSFP_") {
			value = cleanMultiValue(key, value, cleaned)
		}

		if _, exists := cleaned[key]; !exists {
			fieldOrder = append(fieldOrder, key)
		}
		cleaned[key] = value
	}

	// 同步 VEP 和 dbNSFP 的值
	syncAnnotations(cleaned)

	// 重建 INFO 字符串
	var result []string
	for _, key := range fieldOrder {
		if value, ok := cleaned[key]; ok {
			if value == "" {
				result = append(result, key)
			} else {
				result = append(result, fmt.Sprintf("%s=%s", key, value))
			}
		}
	}

	if len(result) == 0 {
		return "."
	}
	return strings.Join(result, ";")
}

// cleanMultiValue 清理多值字段，只取第一个有效值
func cleanMultiValue(key, value string, existing map[string]string) string {
	if value == "" || value == "." {
		return ""
	}

	// 分割多值，只取第一个有效值
	values := strings.Split(value, ",")

	for _, v := range values {
		v = strings.TrimSpace(v)
		if v != "" && v != "." {
			return v // 只返回第一个有效值
		}
	}

	return ""
}

// cleanMultiValueOld 清理多值字段，提取所有有效值（已废弃）
func cleanMultiValueOld(key, value string, existing map[string]string) string {
	if value == "" || value == "." {
		return ""
	}

	// 分割多值
	values := strings.Split(value, ",")
	var validValues []string

	for _, v := range values {
		v = strings.TrimSpace(v)
		if v != "" && v != "." {
			validValues = append(validValues, v)
		}
	}

	if len(validValues) == 0 {
		return ""
	}

	// 如果只有一个有效值，直接返回
	if len(validValues) == 1 {
		return validValues[0]
	}

	// 如果有多个有效值，去重并返回
	seen := make(map[string]bool)
	var unique []string
	for _, v := range validValues {
		if !seen[v] {
			seen[v] = true
			unique = append(unique, v)
		}
	}

	// 如果去重后只有一个值，直接返回
	if len(unique) == 1 {
		return unique[0]
	}

	// 返回逗号分隔的唯一值
	return strings.Join(unique, ",")
}

// syncAnnotations 同步 VEP 和 dbNSFP 注释
func syncAnnotations(info map[string]string) {
	// 定义对应关系
	pairs := []struct {
		vep    string
		dbnsfp string
	}{
		{"VEP_AlphaMissense", "dbNSFP_AlphaMissense"},
		{"VEP_REVEL", "dbNSFP_REVEL"},
		{"VEP_CADD", "dbNSFP_CADD"},
	}

	for _, pair := range pairs {
		vepVal := info[pair.vep]
		dbVal := info[pair.dbnsfp]

		// 如果 VEP 有值，dbNSFP 没有，复制过去
		if vepVal != "" && dbVal == "" {
			info[pair.dbnsfp] = vepVal
		}
		// 如果 dbNSFP 有值，VEP 没有，复制过去
		if dbVal != "" && vepVal == "" {
			info[pair.vep] = dbVal
		}
		// 如果两者都有值，优先使用 VEP 的值更新 dbNSFP
		if vepVal != "" && dbVal != "" {
			// 检查 dbNSFP 的值是否包含 VEP 的值
			dbValues := strings.Split(dbVal, ",")
			foundMatch := false
			for _, dv := range dbValues {
				if dv == vepVal {
					foundMatch = true
					break
				}
			}
			// 如果匹配，使用 VEP 的值
			if foundMatch || len(dbValues) > 1 {
				info[pair.dbnsfp] = vepVal
			}
		}
	}
}

// sortVariants 按染色体和位置排序
func sortVariants(variants []Variant) {
	chromOrder := map[string]int{
		"chr1": 1, "chr2": 2, "chr3": 3, "chr4": 4, "chr5": 5,
		"chr6": 6, "chr7": 7, "chr8": 8, "chr9": 9, "chr10": 10,
		"chr11": 11, "chr12": 12, "chr13": 13, "chr14": 14, "chr15": 15,
		"chr16": 16, "chr17": 17, "chr18": 18, "chr19": 19, "chr20": 20,
		"chr21": 21, "chr22": 22, "chrX": 23, "chrY": 24, "chrM": 25,
	}

	sort.Slice(variants, func(i, j int) bool {
		ci := chromOrder[variants[i].Chrom]
		cj := chromOrder[variants[j].Chrom]
		if ci != cj {
			return ci < cj
		}
		// 比较位置
		var pi, pj int
		fmt.Sscanf(variants[i].Pos, "%d", &pi)
		fmt.Sscanf(variants[j].Pos, "%d", &pj)
		return pi < pj
	})
}

// writeVCF 写入清理后的 VCF 文件
func writeVCF(filename string, variants []Variant, originalHeaders []string) {
	file, err := os.Create(filename)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating %s: %v\n", filename, err)
		os.Exit(1)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)

	// 写入精简的头部
	writeCleanHeader(writer, originalHeaders)

	// 写入变异
	for _, v := range variants {
		line := fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			v.Chrom, v.Pos, v.ID, v.Ref, v.Alt, v.Qual, v.Filter, v.Info)

		if v.Format != "" {
			line += "\t" + v.Format
		}
		if len(v.Samples) > 0 {
			line += "\t" + strings.Join(v.Samples, "\t")
		}

		writer.WriteString(line + "\n")
	}

	writer.Flush()
}

// writeCleanHeader 写入精简的头部
func writeCleanHeader(writer *bufio.Writer, headers []string) {
	// 必需的头部
	writer.WriteString("##fileformat=VCFv4.2\n")
	writer.WriteString("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")

	// 只保留主要染色体的 contig
	mainChroms := regexp.MustCompile(`^##contig=<ID=chr([0-9]+|X|Y|M),`)
	for _, h := range headers {
		if mainChroms.MatchString(h) {
			writer.WriteString(h + "\n")
		}
	}

	// 添加精简的 INFO 定义
	infoHeaders := []string{
		`##INFO=<ID=FINALIST_SOURCE,Number=.,Type=String,Description="Source of finalist variant: PATHOGENIC_HIGH (ClinVar pathogenic or snpEff HIGH impact) or TOP50_SCORED (top 50 by HERITA scoring)">`,
		`##INFO=<ID=HERITA_SCORE,Number=1,Type=Float,Description="HERITA comprehensive pathogenicity score">`,
		`##INFO=<ID=ANN,Number=.,Type=String,Description="snpEff functional annotation">`,
		`##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects">`,
		`##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects">`,
		`##INFO=<ID=VEP_Gene,Number=1,Type=String,Description="Gene symbol from VEP">`,
		`##INFO=<ID=VEP_Consequence,Number=1,Type=String,Description="Most severe consequence from VEP">`,
		`##INFO=<ID=VEP_AlphaMissense,Number=1,Type=Float,Description="AlphaMissense score (0-1, higher=more pathogenic)">`,
		`##INFO=<ID=VEP_AM_Pathogenicity,Number=1,Type=String,Description="AlphaMissense pathogenicity classification">`,
		`##INFO=<ID=VEP_CADD,Number=1,Type=Float,Description="CADD phred score (>20 deleterious, >30 highly deleterious)">`,
		`##INFO=<ID=VEP_REVEL,Number=1,Type=Float,Description="REVEL score (0-1, >0.5 pathogenic)">`,
		`##INFO=<ID=VEP_PhyloP,Number=1,Type=Float,Description="PhyloP conservation score">`,
		`##INFO=<ID=VEP_SpliceAI_Max,Number=1,Type=Float,Description="Maximum SpliceAI delta score">`,
		`##INFO=<ID=VEP_SpliceAI_Details,Number=1,Type=String,Description="SpliceAI delta scores (AG,AL,DG,DL)">`,
		`##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="ClinVar clinical significance">`,
		`##INFO=<ID=ClinVar_CLNDN,Number=.,Type=String,Description="ClinVar disease name">`,
		`##INFO=<ID=ClinVar_CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status">`,
		`##INFO=<ID=ClinVar_ALLELEID,Number=1,Type=Integer,Description="ClinVar allele ID">`,
		`##INFO=<ID=ClinVar_CLNDISDB,Number=.,Type=String,Description="ClinVar disease database IDs">`,
		`##INFO=<ID=dbNSFP_AlphaMissense,Number=1,Type=Float,Description="AlphaMissense score from dbNSFP">`,
		`##INFO=<ID=dbNSFP_REVEL,Number=1,Type=Float,Description="REVEL score from dbNSFP">`,
		`##INFO=<ID=dbNSFP_MetaRNN,Number=1,Type=Float,Description="MetaRNN score from dbNSFP">`,
		`##INFO=<ID=dbNSFP_PrimateAI,Number=1,Type=Float,Description="PrimateAI score from dbNSFP">`,
		`##INFO=<ID=dbNSFP_CADD,Number=1,Type=Float,Description="CADD phred score from dbNSFP">`,
		`##INFO=<ID=dbNSFP_ClinPred,Number=1,Type=Float,Description="ClinPred score from dbNSFP">`,
		`##INFO=<ID=dbNSFP_ESM1b,Number=1,Type=Float,Description="ESM1b score from dbNSFP (negative=deleterious)">`,
		`##INFO=<ID=dbscSNV_ada_score,Number=1,Type=Float,Description="dbscSNV ada score for splicing">`,
		`##INFO=<ID=dbscSNV_rf_score,Number=1,Type=Float,Description="dbscSNV rf score for splicing">`,
		`##INFO=<ID=gnomAD_AF_grpmax,Number=1,Type=Float,Description="Maximum allele frequency across gnomAD populations">`,
	}

	for _, h := range infoHeaders {
		writer.WriteString(h + "\n")
	}

	// 找到原来的 #CHROM 行
	for _, h := range headers {
		if strings.HasPrefix(h, "#CHROM") {
			writer.WriteString(h + "\n")
			break
		}
	}
}
