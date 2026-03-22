package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"
)

// PipelineConfig 流水线配置
type PipelineConfig struct {
	// 输入输出
	InputVCF  string
	OutputDir string
	WorkDir   string

	// 数据集路径
	GnomADDir  string
	ClinVarVCF string
	DbNSFPFile string

	// 阈值参数
	AFThreshold       float64
	SpliceAIThreshold float64
	Top50Count        int

	// 运行选项
	Threads   int
	SkipSteps []string
	DryRun    bool
	Verbose   bool
}

// StepResult 步骤结果
type StepResult struct {
	Name        string
	InputCount  int
	OutputCount int
	FinalCount  int
	Duration    time.Duration
	Success     bool
	Error       error
}

var config PipelineConfig

func main() {
	// 解析命令行参数
	flag.StringVar(&config.InputVCF, "input", "", "输入 VCF 文件 (必需)")
	flag.StringVar(&config.OutputDir, "output", "results", "输出目录")
	flag.StringVar(&config.WorkDir, "workdir", "", "工作目录 (默认: 当前目录)")

	flag.StringVar(&config.GnomADDir, "gnomad", "gnomAD_dataset", "gnomAD 数据目录")
	flag.StringVar(&config.ClinVarVCF, "clinvar", "input_data/clinvar_vcf_GRCh38.vcf.gz", "ClinVar VCF 文件")
	flag.StringVar(&config.DbNSFPFile, "dbnsfp", "input_data/dbNSFP5.3a_grch38_lite.tsv.gz", "dbNSFP 文件")

	flag.Float64Var(&config.AFThreshold, "af", 0.01, "gnomAD AF 阈值")
	flag.Float64Var(&config.SpliceAIThreshold, "spliceai", 0.5, "SpliceAI 阈值")
	flag.IntVar(&config.Top50Count, "top", 50, "打分筛选 TOP N")

	flag.IntVar(&config.Threads, "threads", 4, "线程数")
	skipSteps := flag.String("skip", "", "跳过的步骤 (逗号分隔)")
	flag.BoolVar(&config.DryRun, "dry-run", false, "仅显示将执行的命令")
	flag.BoolVar(&config.Verbose, "verbose", false, "显示详细输出")

	flag.Parse()

	if config.InputVCF == "" {
		fmt.Println("错误: 必须指定输入 VCF 文件")
		flag.Usage()
		os.Exit(1)
	}

	if *skipSteps != "" {
		config.SkipSteps = strings.Split(*skipSteps, ",")
	}

	if config.WorkDir == "" {
		config.WorkDir, _ = os.Getwd()
	}

	// 运行流水线
	runPipeline()
}

func runPipeline() {
	startTime := time.Now()

	printHeader()

	// 检查依赖
	if !checkDependencies() {
		os.Exit(1)
	}

	// 创建输出目录
	os.MkdirAll(config.OutputDir, 0755)

	// 获取基础文件名
	baseName := getBaseName(config.InputVCF)

	// 定义流水线步骤
	var results []StepResult
	var currentInput string
	var finalRoundVCF string

	// ============================================
	// 步骤 1: gnomAD 人群频率过滤
	// ============================================
	step1Output := filepath.Join(config.OutputDir, baseName+"_gnomad_filtered.vcf")
	if !shouldSkip("gnomad") {
		result := runStep("1. gnomAD 频率过滤", func() error {
			return runGnomadFilter(config.InputVCF, step1Output)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 1 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 1: gnomAD 频率过滤")
	}
	currentInput = step1Output

	// ============================================
	// 步骤 2: ClinVar 注释
	// ============================================
	step2Output := filepath.Join(config.OutputDir, baseName+"_clinvar_annotated.vcf")
	if !shouldSkip("clinvar_annotate") {
		result := runStep("2. ClinVar 注释", func() error {
			return runClinvarAnnotate(currentInput, step2Output)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 2 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 2: ClinVar 注释")
	}
	currentInput = step2Output

	// ============================================
	// 步骤 3: ClinVar 致病性分流
	// ============================================
	finalRoundVCF = filepath.Join(config.OutputDir, baseName+"_final_round.vcf")
	step3NextLevel := filepath.Join(config.OutputDir, baseName+"_next_level.vcf")
	if !shouldSkip("clinvar_filter") {
		result := runStep("3. ClinVar 致病性分流", func() error {
			return runClinvarFilter(currentInput, finalRoundVCF, step3NextLevel)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 3 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 3: ClinVar 致病性分流")
	}
	currentInput = step3NextLevel

	// ============================================
	// 步骤 4: snpEff HIGH impact 过滤
	// ============================================
	step4Output := filepath.Join(config.OutputDir, baseName+"_after_snpeff.vcf")
	if !shouldSkip("snpeff") {
		result := runStep("4. snpEff HIGH 过滤", func() error {
			return runSnpeffFilter(currentInput, step4Output, finalRoundVCF)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 4 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 4: snpEff HIGH 过滤")
	}
	currentInput = step4Output

	// ============================================
	// 步骤 5: VEP REST API 注释
	// ============================================
	step5Output := filepath.Join(config.OutputDir, baseName+"_after_vep.vcf")
	if !shouldSkip("vep") {
		result := runStep("5. VEP API 注释", func() error {
			return runVepAnnotate(currentInput, step5Output, finalRoundVCF)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 5 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 5: VEP API 注释")
	}
	currentInput = step5Output

	// ============================================
	// 步骤 6: dbNSFP 深度注释
	// ============================================
	step6Output := filepath.Join(config.OutputDir, baseName+"_dbnsfp_annotated.vcf")
	if !shouldSkip("dbnsfp") {
		result := runStep("6. dbNSFP 深度注释", func() error {
			return runDbnsfpAnnotate(currentInput, step6Output)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 6 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 6: dbNSFP 深度注释")
	}
	currentInput = step6Output

	// ============================================
	// 步骤 7: 综合打分筛选 TOP 50
	// ============================================
	step7Output := filepath.Join(config.OutputDir, baseName+"_top50_finalists.vcf")
	scoringReport := filepath.Join(config.OutputDir, baseName+"_scoring_report.txt")
	if !shouldSkip("scoring") {
		result := runStep("7. 综合打分筛选", func() error {
			return runVariantScorer(currentInput, step7Output, scoringReport)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 7 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 7: 综合打分筛选")
	}

	// ============================================
	// 步骤 8: 合并决赛圈
	// ============================================
	step8Output := filepath.Join(config.OutputDir, baseName+"_finalists_merged.vcf")
	if !shouldSkip("merge") {
		result := runStep("8. 合并决赛圈", func() error {
			return runFinalistsMerger(finalRoundVCF, step7Output, step8Output)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 8 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 8: 合并决赛圈")
	}

	// ============================================
	// 步骤 9: 补充缺失的 dbNSFP 注释
	// ============================================
	step9Output := filepath.Join(config.OutputDir, baseName+"_dbnsfp_supplemented.vcf")
	if !shouldSkip("dbnsfp_supplement") {
		result := runStep("9. 补充 dbNSFP 注释", func() error {
			return runDbnsfpSupplement(step8Output, step9Output)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 9 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 9: 补充 dbNSFP 注释")
	}

	// ============================================
	// 步骤 10: ACMG 分类
	// ============================================
	finalOutput := filepath.Join(config.OutputDir, baseName+"_acmg_classified.vcf")
	acmgReport := filepath.Join(config.OutputDir, baseName+"_acmg_report.md")
	acmgDetailedLog := filepath.Join(config.OutputDir, baseName+"_acmg_detailed_log.txt")
	if !shouldSkip("acmg") {
		result := runStep("10. ACMG 分类", func() error {
			return runAcmgClassifier(step9Output, finalOutput, acmgReport, acmgDetailedLog)
		})
		results = append(results, result)
		if !result.Success {
			printError("步骤 10 失败，终止流水线")
			os.Exit(1)
		}
	} else {
		printSkip("步骤 10: ACMG 分类")
	}

	// 打印总结
	printSummary(results, finalOutput, time.Since(startTime))
}

// ============================================
// 各步骤执行函数
// ============================================

func runGnomadFilter(input, output string) error {
	binPath := getBinPath("gnomad_grpmax_filter")
	args := []string{
		"-input", input,
		"-output", output,
		"-gnomad", config.GnomADDir,
		"-max-af", fmt.Sprintf("%f", config.AFThreshold),
	}
	if config.Threads > 0 {
		args = append(args, "-threads", fmt.Sprintf("%d", config.Threads))
	}
	return runCommand(binPath, args)
}

func runClinvarAnnotate(input, output string) error {
	binPath := getBinPath("clinvar_annotate")
	args := []string{
		"-input", input,
		"-output", output,
		"-clinvar", config.ClinVarVCF,
		"-threads", fmt.Sprintf("%d", config.Threads),
	}
	return runCommand(binPath, args)
}

func runClinvarFilter(input, finalRound, nextLevel string) error {
	binPath := getBinPath("clinvar_filter")
	args := []string{
		"-input", input,
		"-final", finalRound,
		"-next", nextLevel,
	}
	return runCommand(binPath, args)
}

func runSnpeffFilter(input, output, finalRound string) error {
	binPath := getBinPath("snpeff_filter")
	args := []string{
		"-input", input,
		"-next", output,
		"-final", finalRound,
	}
	return runCommand(binPath, args)
}

func runVepAnnotate(input, output, finalRound string) error {
	binPath := getBinPath("vep_annotate")
	args := []string{
		"-input", input,
		"-next", output,
		"-final", finalRound,
		"-spliceai-thresh", fmt.Sprintf("%f", config.SpliceAIThreshold),
	}
	return runCommand(binPath, args)
}

func runDbnsfpAnnotate(input, output string) error {
	binPath := getBinPath("dbnsfp_annotate")
	args := []string{
		"-input", input,
		"-output", output,
		"-dbnsfp", config.DbNSFPFile,
		"-threads", fmt.Sprintf("%d", config.Threads),
	}
	return runCommand(binPath, args)
}

func runVariantScorer(input, output, report string) error {
	// 使用 Python 脚本
	scriptPath := filepath.Join(config.WorkDir, "scripts", "variant_scorer.py")
	args := []string{
		scriptPath,
		"-i", input,
		"--vcf-output", output,
		"-o", report,
		"-n", fmt.Sprintf("%d", config.Top50Count),
		"--verbose",
	}
	return runCommand("python3", args)
}

func runFinalistsMerger(finalRound, top50, output string) error {
	binPath := getBinPath("finalists_merger")
	args := []string{
		"-final-round", finalRound,
		"-top50", top50,
		"-output", output,
	}
	return runCommand(binPath, args)
}

func runDbnsfpSupplement(input, output string) error {
	// 使用 Python 脚本补充缺失的 dbNSFP 注释
	scriptPath := filepath.Join(config.WorkDir, "scripts", "dbnsfp_supplement.py")
	args := []string{
		scriptPath,
		"-i", input,
		"-o", output,
		"-d", config.DbNSFPFile,
	}
	if config.Verbose {
		args = append(args, "-v")
	}
	return runCommand("python3", args)
}

func runAcmgClassifier(input, output, report, detailedLog string) error {
	// 使用 Python 脚本进行 ACMG 分类 (v2版本)
	scriptPath := filepath.Join(config.WorkDir, "scripts", "acmg_classifier_v2.py")
	intervarDbPath := filepath.Join(config.WorkDir, "InterVar-master", "intervardb")
	args := []string{
		scriptPath,
		"-i", input,
		"-o", output,
		"-r", report,
		"--intervar-db", intervarDbPath,
	}
	if detailedLog != "" {
		args = append(args, "--detailed-log", detailedLog)
	}
	if config.Verbose {
		args = append(args, "-v")
	}
	return runCommand("python3", args)
}

// ============================================
// 辅助函数
// ============================================

func runStep(name string, fn func() error) StepResult {
	fmt.Printf("\n%s\n", strings.Repeat("─", 60))
	fmt.Printf("▶ %s\n", name)
	fmt.Printf("%s\n", strings.Repeat("─", 60))

	start := time.Now()
	err := fn()
	duration := time.Since(start)

	result := StepResult{
		Name:     name,
		Duration: duration,
		Success:  err == nil,
		Error:    err,
	}

	if err != nil {
		fmt.Printf("✗ 失败: %v (耗时 %s)\n", err, duration.Round(time.Second))
	} else {
		fmt.Printf("✓ 完成 (耗时 %s)\n", duration.Round(time.Second))
	}

	return result
}

func runCommand(name string, args []string) error {
	if config.Verbose || config.DryRun {
		fmt.Printf("  命令: %s %s\n", name, strings.Join(args, " "))
	}

	if config.DryRun {
		return nil
	}

	cmd := exec.Command(name, args...)
	cmd.Dir = config.WorkDir
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	return cmd.Run()
}

func getBinPath(tool string) string {
	return filepath.Join(config.WorkDir, "herita-pipeline-go", "bin", tool)
}

func getBaseName(path string) string {
	base := filepath.Base(path)
	ext := filepath.Ext(base)
	return strings.TrimSuffix(base, ext)
}

func shouldSkip(step string) bool {
	for _, s := range config.SkipSteps {
		if strings.TrimSpace(s) == step {
			return true
		}
	}
	return false
}

func checkDependencies() bool {
	fmt.Println("\n检查依赖...")

	tools := []string{
		"gnomad_grpmax_filter",
		"clinvar_annotate",
		"clinvar_filter",
		"snpeff_filter",
		"vep_annotate",
		"dbnsfp_annotate",
		"finalists_merger",
	}

	allFound := true
	for _, tool := range tools {
		binPath := getBinPath(tool)
		if _, err := os.Stat(binPath); os.IsNotExist(err) {
			fmt.Printf("  ✗ 缺少: %s\n", tool)
			allFound = false
		} else {
			if config.Verbose {
				fmt.Printf("  ✓ 找到: %s\n", tool)
			}
		}
	}

	// 检查 Python 脚本
	scriptPath := filepath.Join(config.WorkDir, "scripts", "variant_scorer.py")
	if _, err := os.Stat(scriptPath); os.IsNotExist(err) {
		fmt.Printf("  ✗ 缺少: scripts/variant_scorer.py\n")
		allFound = false
	}

	dbnsfpSupplementPath := filepath.Join(config.WorkDir, "scripts", "dbnsfp_supplement.py")
	if _, err := os.Stat(dbnsfpSupplementPath); os.IsNotExist(err) {
		fmt.Printf("  ✗ 缺少: scripts/dbnsfp_supplement.py\n")
		allFound = false
	}

	acmgScriptPath := filepath.Join(config.WorkDir, "scripts", "acmg_classifier_v2.py")
	if _, err := os.Stat(acmgScriptPath); os.IsNotExist(err) {
		fmt.Printf("  ✗ 缺少: scripts/acmg_classifier_v2.py\n")
		allFound = false
	}

	// 检查 InterVar 知识库
	intervarDbPath := filepath.Join(config.WorkDir, "InterVar-master", "intervardb")
	if _, err := os.Stat(intervarDbPath); os.IsNotExist(err) {
		fmt.Printf("  ✗ 缺少: InterVar-master/intervardb 知识库\n")
		allFound = false
	}

	// 检查外部工具
	externalTools := []string{"vcfanno", "bcftools", "python3"}
	for _, tool := range externalTools {
		if _, err := exec.LookPath(tool); err != nil {
			fmt.Printf("  ✗ 缺少系统工具: %s\n", tool)
			allFound = false
		}
	}

	if !allFound {
		fmt.Println("\n请先编译所有工具或安装缺失的依赖")
		return false
	}

	fmt.Println("  ✓ 所有依赖检查通过")
	return true
}

func printHeader() {
	fmt.Println()
	fmt.Println("╔════════════════════════════════════════════════════════════╗")
	fmt.Println("║         HERITA 变异注释与筛选流水线 v1.0                    ║")
	fmt.Println("╠════════════════════════════════════════════════════════════╣")
	fmt.Printf("║ 输入文件: %-48s ║\n", truncateString(config.InputVCF, 48))
	fmt.Printf("║ 输出目录: %-48s ║\n", config.OutputDir)
	fmt.Printf("║ AF 阈值:  %-48s ║\n", fmt.Sprintf("%.4f", config.AFThreshold))
	fmt.Printf("║ 线程数:   %-48s ║\n", fmt.Sprintf("%d", config.Threads))
	fmt.Println("╚════════════════════════════════════════════════════════════╝")
}

func printSkip(step string) {
	fmt.Printf("\n⏭ 跳过: %s\n", step)
}

func printError(msg string) {
	fmt.Printf("\n❌ 错误: %s\n", msg)
}

func printSummary(results []StepResult, finalOutput string, totalDuration time.Duration) {
	fmt.Println()
	fmt.Println("╔════════════════════════════════════════════════════════════╗")
	fmt.Println("║                    流水线执行总结                           ║")
	fmt.Println("╠════════════════════════════════════════════════════════════╣")

	successCount := 0
	for _, r := range results {
		status := "✓"
		if !r.Success {
			status = "✗"
		} else {
			successCount++
		}
		fmt.Printf("║ %s %-45s %8s ║\n", status, truncateString(r.Name, 45), r.Duration.Round(time.Second))
	}

	fmt.Println("╠════════════════════════════════════════════════════════════╣")
	fmt.Printf("║ 成功步骤: %d/%d                                              ║\n", successCount, len(results))
	fmt.Printf("║ 总耗时: %-50s ║\n", totalDuration.Round(time.Second))
	fmt.Println("╠════════════════════════════════════════════════════════════╣")

	// 统计最终结果
	if _, err := os.Stat(finalOutput); err == nil {
		count := countVariants(finalOutput)
		fmt.Printf("║ 最终决赛圈变异: %-42d ║\n", count)
		fmt.Printf("║ 输出文件: %-48s ║\n", truncateString(finalOutput, 48))
	}

	fmt.Println("╚════════════════════════════════════════════════════════════╝")
}

func countVariants(vcfPath string) int {
	file, err := os.Open(vcfPath)
	if err != nil {
		return 0
	}
	defer file.Close()

	count := 0
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if !strings.HasPrefix(line, "#") {
			count++
		}
	}
	return count
}

func truncateString(s string, maxLen int) string {
	if len(s) <= maxLen {
		return s
	}
	return "..." + s[len(s)-maxLen+3:]
}
