package pipeline

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"herita-pipeline-go/config"
	"herita-pipeline-go/internal/vcf"
)

// Step9Result Step9 的运行结果
type Step9Result struct {
	TotalCount      int
	PathogenicCount int
	LikelyPathCount int
	VUSCount        int
	LikelyBenCount  int
	BenignCount     int
	Duration        time.Duration
}

// ACMGKnowledgeBase ACMG 知识库
type ACMGKnowledgeBase struct {
	LOFGenes    map[string]bool            // PVS1: LoF 不耐受基因
	PS1Changes  map[string]bool            // PS1: 已知致病氨基酸改变
	PM1Domains  map[string][]string        // PM1: 功能域
	PP2Genes    map[string]bool            // PP2: 错义致病基因
	BP1Genes    map[string]bool            // BP1: 截断良性基因
	BS2Variants map[string]bool            // BS2: 健康人群变异
	PS4Variants map[string]bool            // PS4: 病例对照显著变异
}

// Step9ACMGClassification 执行 Step9: ACMG 分类
func Step9ACMGClassification(cfg *config.Config) (*Step9Result, error) {
	fmt.Println("\n" + strings.Repeat("=", 60))
	fmt.Println("Step 9: InterVar-style ACMG 分类")
	fmt.Println(strings.Repeat("=", 60))

	startTime := time.Now()

	// 1. 加载知识库
	fmt.Println("\n  [1/3] 加载 InterVar 知识库...")
	kb, err := loadKnowledgeBase(cfg.InterVarDB)
	if err != nil {
		return nil, fmt.Errorf("加载知识库失败: %w", err)
	}

	// 2. 读取变异并评估
	fmt.Println("\n  [2/3] 评估 ACMG 证据...")
	parser, err := vcf.NewParser(cfg.FinalIntegrated)
	if err != nil {
		return nil, err
	}
	defer parser.Close()

	var variants []*vcf.Variant
	for {
		v, err := parser.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			continue
		}
		variants = append(variants, v)
	}
	fmt.Printf("    读取 %d 个变异\n", len(variants))

	// 评估每个变异
	type ClassifiedVariant struct {
		Variant        *vcf.Variant
		PathEvidence   []string
		BenEvidence    []string
		Classification string
		Reason         string
	}

	var results []ClassifiedVariant
	result := &Step9Result{TotalCount: len(variants)}

	for _, v := range variants {
		pathEv, benEv := evaluateACMGEvidence(v, kb)
		classification, reason := classifyACMG(pathEv, benEv)

		results = append(results, ClassifiedVariant{
			Variant:        v,
			PathEvidence:   pathEv,
			BenEvidence:    benEv,
			Classification: classification,
			Reason:         reason,
		})

		switch classification {
		case "Pathogenic":
			result.PathogenicCount++
		case "Likely_pathogenic":
			result.LikelyPathCount++
		case "Uncertain_significance":
			result.VUSCount++
		case "Likely_benign":
			result.LikelyBenCount++
		case "Benign":
			result.BenignCount++
		}
	}

	// 3. 写入结果
	fmt.Println("\n  [3/3] 写入结果...")
	writer, err := vcf.NewWriter(cfg.FinalACMG)
	if err != nil {
		return nil, err
	}
	defer writer.Close()

	// 写入 header
	cleanHeader := vcf.CreateCleanHeader(parser.Header.SampleLine)
	writer.WriteHeader(cleanHeader)

	for _, cv := range results {
		v := cv.Variant
		v.SetInfo("ACMG_InterVar", cv.Classification)
		v.SetInfo("ACMG_PathEvidence", strings.Join(cv.PathEvidence, "|"))
		v.SetInfo("ACMG_BenEvidence", strings.Join(cv.BenEvidence, "|"))
		v.SetInfo("ACMG_Reason", strings.ReplaceAll(cv.Reason, " ", "_"))
		writer.WriteVariant(v)
	}

	result.Duration = time.Since(startTime)

	// 打印统计
	fmt.Println("\n" + strings.Repeat("-", 50))
	fmt.Println("ACMG 分类结果统计:")
	fmt.Println(strings.Repeat("-", 50))
	printClassificationBar("Pathogenic", result.PathogenicCount, result.TotalCount)
	printClassificationBar("Likely_pathogenic", result.LikelyPathCount, result.TotalCount)
	printClassificationBar("Uncertain_significance", result.VUSCount, result.TotalCount)
	printClassificationBar("Likely_benign", result.LikelyBenCount, result.TotalCount)
	printClassificationBar("Benign", result.BenignCount, result.TotalCount)
	fmt.Println(strings.Repeat("-", 50))
	fmt.Printf("总计: %d 个变异\n", result.TotalCount)

	fmt.Printf("\n✅ Step 9 完成！耗时: %.2f秒\n", result.Duration.Seconds())

	return result, nil
}

func printClassificationBar(name string, count, total int) {
	pct := float64(count) / float64(total) * 100
	barLen := int(pct / 2) // 50个字符表示100%
	bar := strings.Repeat("█", barLen)
	fmt.Printf("  %-25s: %3d (%5.1f%%) %s\n", name, count, pct, bar)
}

// loadKnowledgeBase 加载 InterVar 知识库
func loadKnowledgeBase(dbDir string) (*ACMGKnowledgeBase, error) {
	kb := &ACMGKnowledgeBase{
		LOFGenes:    make(map[string]bool),
		PS1Changes:  make(map[string]bool),
		PM1Domains:  make(map[string][]string),
		PP2Genes:    make(map[string]bool),
		BP1Genes:    make(map[string]bool),
		BS2Variants: make(map[string]bool),
		PS4Variants: make(map[string]bool),
	}

	// 加载 PVS1 基因列表
	if lines, err := readLines(filepath.Join(dbDir, "PVS1.LOF.genes.hg38")); err == nil {
		for _, line := range lines {
			kb.LOFGenes[strings.TrimSpace(line)] = true
		}
		fmt.Printf("    加载 %d 个 LOF 不耐受基因 (PVS1)\n", len(kb.LOFGenes))
	}

	// 加载 PS1 氨基酸改变
	if lines, err := readLines(filepath.Join(dbDir, "PS1.AA.change.patho.hg38")); err == nil {
		for _, line := range lines {
			kb.PS1Changes[strings.TrimSpace(line)] = true
		}
		fmt.Printf("    加载 %d 个已知致病氨基酸改变 (PS1)\n", len(kb.PS1Changes))
	}

	// 加载 PM1 功能域
	if file, err := os.Open(filepath.Join(dbDir, "PM1.domain.bindingsites")); err == nil {
		defer file.Close()
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			fields := strings.Split(scanner.Text(), "\t")
			if len(fields) >= 2 {
				gene := fields[0]
				kb.PM1Domains[gene] = append(kb.PM1Domains[gene], fields[1])
			}
		}
		fmt.Printf("    加载 %d 个功能域记录 (PM1)\n", len(kb.PM1Domains))
	}

	// 加载 PP2 基因
	if lines, err := readLines(filepath.Join(dbDir, "PP2.genes.hg38")); err == nil {
		for _, line := range lines {
			kb.PP2Genes[strings.TrimSpace(line)] = true
		}
		fmt.Printf("    加载 %d 个错义致病基因 (PP2)\n", len(kb.PP2Genes))
	}

	// 加载 BP1 基因
	if lines, err := readLines(filepath.Join(dbDir, "BP1.genes.hg38")); err == nil {
		for _, line := range lines {
			kb.BP1Genes[strings.TrimSpace(line)] = true
		}
		fmt.Printf("    加载 %d 个截断良性基因 (BP1)\n", len(kb.BP1Genes))
	}

	// 加载 BS2 变异
	if file, err := os.Open(filepath.Join(dbDir, "BS2_hom_het.hg38")); err == nil {
		defer file.Close()
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			kb.BS2Variants[strings.TrimSpace(scanner.Text())] = true
		}
		fmt.Printf("    加载 %d 个健康人群变异 (BS2)\n", len(kb.BS2Variants))
	}

	// 加载 PS4 变异
	if file, err := os.Open(filepath.Join(dbDir, "PS4.variants.hg38")); err == nil {
		defer file.Close()
		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			kb.PS4Variants[strings.TrimSpace(scanner.Text())] = true
		}
		fmt.Printf("    加载 %d 个病例对照显著变异 (PS4)\n", len(kb.PS4Variants))
	}

	return kb, nil
}

func readLines(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" && !strings.HasPrefix(line, "#") {
			lines = append(lines, line)
		}
	}
	return lines, scanner.Err()
}

// evaluateACMGEvidence 评估 ACMG 证据
func evaluateACMGEvidence(v *vcf.Variant, kb *ACMGKnowledgeBase) (pathogenic []string, benign []string) {
	gene := v.GetInfoString("genename")
	aaref := v.GetInfoString("aaref")
	aaalt := v.GetInfoString("aaalt")

	// 致病证据
	// PVS1: LoF 变异在 LoF 不耐受基因中
	if kb.LOFGenes[gene] {
		// 检查是否是 LoF 变异 (简化判断)
		if aaalt == "*" || aaalt == "X" {
			pathogenic = append(pathogenic, "PVS1")
		}
	}

	// PS1: 相同氨基酸改变已知致病
	if aaref != "" && aaalt != "" {
		aaChange := fmt.Sprintf("%s:%s:%s:%s", gene, aaref, strconv.Itoa(v.Pos), aaalt)
		if kb.PS1Changes[aaChange] {
			pathogenic = append(pathogenic, "PS1")
		}
	}

	// PS4: 病例对照显著
	varKey := fmt.Sprintf("%s:%d:%s:%s", v.Chrom, v.Pos, v.Ref, v.Alt)
	if kb.PS4Variants[varKey] {
		pathogenic = append(pathogenic, "PS4")
	}

	// PM1: 功能域
	if _, ok := kb.PM1Domains[gene]; ok {
		pathogenic = append(pathogenic, "PM1")
	}

	// PM2: 低频变异
	if af, ok := v.GetInfoFloat("gnomAD4.1_joint_POPMAX_AF"); ok && af < 0.0001 {
		pathogenic = append(pathogenic, "PM2")
	}

	// PP2: 错义致病基因
	if kb.PP2Genes[gene] && aaref != "" && aaalt != "" && aaalt != "*" {
		pathogenic = append(pathogenic, "PP2")
	}

	// PP3: 多个计算预测为致病
	ppCount := 0
	if score, ok := v.GetInfoFloat("REVEL_score"); ok && score >= 0.5 {
		ppCount++
	}
	if score, ok := v.GetInfoFloat("MetaRNN_score"); ok && score >= 0.5 {
		ppCount++
	}
	if score, ok := v.GetInfoFloat("AlphaMissense_score"); ok && score >= 0.564 {
		ppCount++
	}
	if ppCount >= 2 {
		pathogenic = append(pathogenic, "PP3")
	}

	// PP5: ClinVar 报告致病
	clnsig := strings.ToLower(v.GetInfoString("clinvar_clnsig"))
	if strings.Contains(clnsig, "pathogenic") && !strings.Contains(clnsig, "benign") {
		pathogenic = append(pathogenic, "PP5")
	}

	// 良性证据
	// BA1: 高频变异
	if af, ok := v.GetInfoFloat("gnomAD4.1_joint_POPMAX_AF"); ok && af > 0.05 {
		benign = append(benign, "BA1")
	}

	// BS1: 频率高于预期
	if af, ok := v.GetInfoFloat("gnomAD4.1_joint_POPMAX_AF"); ok && af > 0.01 {
		benign = append(benign, "BS1")
	}

	// BS2: 健康人群中观察到
	if kb.BS2Variants[varKey] {
		benign = append(benign, "BS2")
	}

	// BP1: 截断良性基因
	if kb.BP1Genes[gene] {
		benign = append(benign, "BP1")
	}

	// BP4: 计算预测为良性
	bpCount := 0
	if score, ok := v.GetInfoFloat("REVEL_score"); ok && score < 0.3 {
		bpCount++
	}
	if score, ok := v.GetInfoFloat("MetaRNN_score"); ok && score < 0.3 {
		bpCount++
	}
	if bpCount >= 2 {
		benign = append(benign, "BP4")
	}

	// BP6: ClinVar 报告良性
	if strings.Contains(clnsig, "benign") && !strings.Contains(clnsig, "pathogenic") {
		benign = append(benign, "BP6")
	}

	return pathogenic, benign
}

// classifyACMG 根据证据分类
func classifyACMG(pathEvidence, benEvidence []string) (string, string) {
	// 计算证据强度
	pvs := countPrefix(pathEvidence, "PVS")
	ps := countPrefix(pathEvidence, "PS")
	pm := countPrefix(pathEvidence, "PM")
	pp := countPrefix(pathEvidence, "PP")

	ba := countPrefix(benEvidence, "BA")
	bs := countPrefix(benEvidence, "BS")
	bp := countPrefix(benEvidence, "BP")

	// ACMG 分类规则
	// Pathogenic
	if pvs >= 1 && (ps >= 1 || pm >= 2 || (pm >= 1 && pp >= 1) || pp >= 2) {
		return "Pathogenic", "PVS1 + supporting evidence"
	}
	if ps >= 2 {
		return "Pathogenic", "2+ Strong pathogenic"
	}
	if ps >= 1 && (pm >= 3 || (pm >= 2 && pp >= 2) || (pm >= 1 && pp >= 4)) {
		return "Pathogenic", "PS + multiple moderate/supporting"
	}

	// Likely Pathogenic
	if pvs >= 1 && pm >= 1 {
		return "Likely_pathogenic", "PVS1 + PM"
	}
	if ps >= 1 && (pm >= 1 || pm >= 2) {
		return "Likely_pathogenic", "PS + PM"
	}
	if ps >= 1 && pp >= 2 {
		return "Likely_pathogenic", "PS + 2PP"
	}
	if pm >= 3 {
		return "Likely_pathogenic", "3+ Moderate pathogenic"
	}
	if pm >= 2 && pp >= 2 {
		return "Likely_pathogenic", "2PM + 2PP"
	}
	if pm >= 1 && pp >= 4 {
		return "Likely_pathogenic", "PM + 4PP"
	}

	// Benign
	if ba >= 1 {
		return "Benign", "BA1 (high frequency)"
	}
	if bs >= 2 {
		return "Benign", "2+ Strong benign"
	}

	// Likely Benign
	if bs >= 1 && bp >= 1 {
		return "Likely_benign", "BS + BP"
	}
	if bp >= 2 {
		return "Likely_benign", "2+ Supporting benign"
	}

	// VUS
	return "Uncertain_significance", "Insufficient evidence"
}

func countPrefix(evidence []string, prefix string) int {
	count := 0
	for _, e := range evidence {
		if strings.HasPrefix(e, prefix) {
			count++
		}
	}
	return count
}
