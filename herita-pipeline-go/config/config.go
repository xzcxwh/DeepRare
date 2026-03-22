package config

import (
	"path/filepath"
	"runtime"
)

// Config 流水线配置
type Config struct {
	// 项目根目录
	ProjectRoot string

	// 输入文件
	InputVCF    string
	DBNSFPFile  string
	DBscSNVFile string

	// 中间结果文件
	VCF1 string // Step1 输出
	VCF2 string // Step2 输出
	VCF3 string // Step3 输出
	VCF4 string // Step4 输出
	VCF5 string // Step5 输出
	VCF6 string // Step6 输出

	// 决赛圈文件
	FinalClinVar string
	FinalDbscSNV string
	FinalTop     string

	// 最终输出
	FinalIntegrated string
	FinalACMG       string

	// InterVar 数据库目录
	InterVarDB string

	// 并发设置
	NumWorkers int

	// 筛选阈值
	GnomADMaxAF      float64 // gnomAD 最大等位基因频率 (Python版本使用 0.01)
	MinQUAL          float64 // 最小 QUAL 值
	MinDP            int     // 最小深度
	MinGQ            int     // 最小基因型质量
	DbscSNVThreshold float64 // dbscSNV 阈值
	MinScoreCount    int     // 最少满足的评分条件数
	TopN             int     // Top N 变异数
}

// NewConfig 创建默认配置
func NewConfig(projectRoot string) *Config {
	resultsDir := filepath.Join(projectRoot, "results")
	inputDir := filepath.Join(projectRoot, "input_data")

	return &Config{
		ProjectRoot: projectRoot,

		// 输入文件
		InputVCF:    filepath.Join(inputDir, "HG001.vcf"),
		DBNSFPFile:  filepath.Join(inputDir, "dbNSFP5.3a_grch38_lite.tsv.gz"),
		DBscSNVFile: filepath.Join(projectRoot, "dbscSNV1", "dbscSNV1.1_hg38_sorted.tsv"),

		// 中间结果
		VCF1: filepath.Join(resultsDir, "vcf1.vcf"),
		VCF2: filepath.Join(resultsDir, "vcf2.vcf"),
		VCF3: filepath.Join(resultsDir, "vcf3.vcf"),
		VCF4: filepath.Join(resultsDir, "vcf4.vcf"),
		VCF5: filepath.Join(resultsDir, "vcf5.vcf"),
		VCF6: filepath.Join(resultsDir, "vcf6.vcf"),

		// 决赛圈
		FinalClinVar: filepath.Join(resultsDir, "final_clinvar.vcf"),
		FinalDbscSNV: filepath.Join(resultsDir, "final_dbscSNV.vcf"),
		FinalTop:     filepath.Join(resultsDir, "final_top.vcf"),

		// 最终输出
		FinalIntegrated: filepath.Join(resultsDir, "final_integrated.vcf"),
		FinalACMG:       filepath.Join(resultsDir, "final_intervar_classified.vcf"),

		// InterVar 数据库
		InterVarDB: filepath.Join(projectRoot, "InterVar-master", "intervardb"),

		// 并发设置
		NumWorkers: runtime.NumCPU(),

		// 筛选阈值
		GnomADMaxAF:      0.01,  // 1%
		MinQUAL:          30.0,
		MinDP:            30,
		MinGQ:            60,
		DbscSNVThreshold: 0.6,
		MinScoreCount:    2,
		TopN:             50,
	}
}
