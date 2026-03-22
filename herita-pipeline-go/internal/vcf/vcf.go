package vcf

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// Variant 表示一个 VCF 变异记录
type Variant struct {
	Chrom   string
	Pos     int
	ID      string
	Ref     string
	Alt     string
	Qual    float64
	Filter  string
	Info    map[string]string
	Format  string
	Samples []string
	RawLine string
}

// Key 返回变异的唯一标识
func (v *Variant) Key() string {
	return fmt.Sprintf("%s_%d_%s_%s", v.Chrom, v.Pos, v.Ref, v.Alt)
}

// ChromNum 返回染色体排序用的数字
func (v *Variant) ChromNum() int {
	chrom := v.Chrom
	if strings.HasPrefix(chrom, "chr") {
		chrom = chrom[3:]
	}
	switch chrom {
	case "X":
		return 23
	case "Y":
		return 24
	case "M", "MT":
		return 25
	default:
		num, err := strconv.Atoi(chrom)
		if err != nil {
			return 26
		}
		return num
	}
}

// GetInfoFloat 获取 INFO 字段的浮点值
func (v *Variant) GetInfoFloat(key string) (float64, bool) {
	val, ok := v.Info[key]
	if !ok || val == "" || val == "." {
		return 0, false
	}
	// 处理多个值的情况，取第一个非点的值
	parts := strings.Split(val, ",")
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" && p != "." {
			f, err := strconv.ParseFloat(p, 64)
			if err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

// GetInfoFloatMax 获取 INFO 字段的最大浮点值（处理多值字段）
func (v *Variant) GetInfoFloatMax(key string) float64 {
	val, ok := v.Info[key]
	if !ok || val == "" || val == "." {
		return -1e30 // 返回极小值表示无效
	}
	
	maxVal := -1e30
	found := false
	parts := strings.Split(val, ",")
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" && p != "." {
			f, err := strconv.ParseFloat(p, 64)
			if err == nil {
				if f > maxVal {
					maxVal = f
				}
				found = true
			}
		}
	}
	if !found {
		return -1e30
	}
	return maxVal
}

// GetInfoFloatMin 获取 INFO 字段的最小浮点值（处理多值字段）
func (v *Variant) GetInfoFloatMin(key string) float64 {
	val, ok := v.Info[key]
	if !ok || val == "" || val == "." {
		return 1e30 // 返回极大值表示无效
	}
	
	minVal := 1e30
	found := false
	parts := strings.Split(val, ",")
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" && p != "." {
			f, err := strconv.ParseFloat(p, 64)
			if err == nil {
				if f < minVal {
					minVal = f
				}
				found = true
			}
		}
	}
	if !found {
		return 1e30
	}
	return minVal
}

// GetInfoString 获取 INFO 字段的字符串值
func (v *Variant) GetInfoString(key string) string {
	val, ok := v.Info[key]
	if !ok || val == "" || val == "." {
		return ""
	}
	return val
}

// SetInfo 设置 INFO 字段
func (v *Variant) SetInfo(key, value string) {
	if v.Info == nil {
		v.Info = make(map[string]string)
	}
	v.Info[key] = value
}

// 定义需要取最大值的评分字段
var maxScoreFields = map[string]bool{
	"AlphaMissense_score": true,
	"REVEL_score":         true,
	"MetaRNN_score":       true,
	"PrimateAI_score":     true,
	"CADD_phred":          true,
	"ClinPred_score":      true,
	"DBSCSNV_ADA":         true,
	"DBSCSNV_RF":          true,
}

// 定义需要取最小值的评分字段（ESM1b_score越小表示影响越大）
var minScoreFields = map[string]bool{
	"ESM1b_score": true,
}

// 定义需要取唯一值（去重）的字段
var uniqueFields = map[string]bool{
	"aaref":                     true,
	"aaalt":                     true,
	"rs_dbSNP":                  true,
	"genename":                  true,
	"clinvar_id":                true,
	"clinvar_clnsig":            true,
	"clinvar_trait":             true,
	"clinvar_review":            true,
	"clinvar_hgvs":              true,
	"clinvar_var_source":        true,
	"clinvar_MedGen_id":         true,
	"clinvar_OMIM_id":           true,
	"clinvar_Orphanet_id":       true,
	"gnomAD4.1_joint":           true,
	"gnomAD4.1_joint_POPMAX_AF": true,
}

// CleanInfoValue 清理 INFO 字段值
// 根据字段类型进行不同处理：
// - 评分字段：取最大值或最小值
// - 唯一值字段：去重，取第一个非空值
// - 其他字段：取第一个非空值
func CleanInfoValue(fieldName, value string) string {
	if value == "" || value == "." {
		return "."
	}

	// 检查是否全是点号
	cleaned := strings.ReplaceAll(value, ",", "")
	cleaned = strings.ReplaceAll(cleaned, ";", "")
	cleaned = strings.ReplaceAll(cleaned, "|", "")
	cleaned = strings.ReplaceAll(cleaned, ".", "")
	cleaned = strings.TrimSpace(cleaned)
	if cleaned == "" {
		return "."
	}

	// 分割值（支持逗号、分号、竖线）
	var values []string
	if strings.Contains(value, ";") {
		values = strings.Split(value, ";")
	} else if strings.Contains(value, ",") {
		values = strings.Split(value, ",")
	} else if strings.Contains(value, "|") {
		values = strings.Split(value, "|")
	} else {
		return value // 单个值直接返回
	}

	// 清理每个值，收集非空值
	var cleanedValues []string
	for _, v := range values {
		v = strings.TrimSpace(v)
		if v != "" && v != "." {
			cleanedValues = append(cleanedValues, v)
		}
	}

	if len(cleanedValues) == 0 {
		return "."
	}

	// 处理需要取最大值的评分字段
	if maxScoreFields[fieldName] {
		maxVal := -1e30
		found := false
		for _, v := range cleanedValues {
			if f, err := strconv.ParseFloat(v, 64); err == nil {
				if f > maxVal {
					maxVal = f
				}
				found = true
			}
		}
		if found {
			return strconv.FormatFloat(maxVal, 'f', -1, 64)
		}
		return "."
	}

	// 处理需要取最小值的评分字段
	if minScoreFields[fieldName] {
		minVal := 1e30
		found := false
		for _, v := range cleanedValues {
			if f, err := strconv.ParseFloat(v, 64); err == nil {
				if f < minVal {
					minVal = f
				}
				found = true
			}
		}
		if found {
			return strconv.FormatFloat(minVal, 'f', -1, 64)
		}
		return "."
	}

	// 处理需要取唯一值的字段
	if uniqueFields[fieldName] {
		seen := make(map[string]bool)
		var unique []string
		for _, v := range cleanedValues {
			if !seen[v] {
				unique = append(unique, v)
				seen[v] = true
			}
		}
		if len(unique) == 1 {
			return unique[0]
		}
		if len(unique) > 1 {
			return strings.Join(unique, "|")
		}
		return "."
	}

	// 默认：取第一个非空值
	if len(cleanedValues) > 0 {
		return cleanedValues[0]
	}
	return "."
}

// BuildInfoString 构建 INFO 字段字符串（按照固定顺序）
func (v *Variant) BuildInfoString() string {
	if len(v.Info) == 0 {
		return "."
	}
	
	// 定义字段输出顺序（与 Python 版本一致）
	fieldOrder := []string{
		"gnomAD4.1_joint_POPMAX_AF", "gnomAD4.1_joint",
		"DBSCSNV_ADA", "DBSCSNV_RF",
		"MetaRNN_score", "REVEL_score", "PrimateAI_score", "ClinPred_score",
		"ESM1b_score", "AlphaMissense_score", "CADD_phred",
		"aaref", "aaalt", "rs_dbSNP", "genename",
		"clinvar_id", "clinvar_clnsig", "clinvar_trait", "clinvar_review",
		"clinvar_hgvs", "clinvar_var_source", "clinvar_MedGen_id",
		"clinvar_OMIM_id", "clinvar_Orphanet_id",
		"ACMG_InterVar", "ACMG_PathEvidence", "ACMG_BenEvidence", "ACMG_Reason",
	}
	
	parts := make([]string, 0, len(v.Info))
	usedKeys := make(map[string]bool)
	
	// 首先按固定顺序输出
	for _, key := range fieldOrder {
		if val, ok := v.Info[key]; ok {
			if val == "" {
				val = "."
			}
			parts = append(parts, key+"="+val)
			usedKeys[key] = true
		}
	}
	
	// 然后输出其他未在顺序列表中的字段
	for k, val := range v.Info {
		if !usedKeys[k] {
			if val == "" {
				parts = append(parts, k)
			} else {
				parts = append(parts, k+"="+val)
			}
		}
	}
	
	if len(parts) == 0 {
		return "."
	}
	return strings.Join(parts, ";")
}

// ToLine 将变异转换为 VCF 行
func (v *Variant) ToLine() string {
	qual := "."
	if v.Qual >= 0 {
		qual = strconv.FormatFloat(v.Qual, 'f', -1, 64)
	}

	fields := []string{
		v.Chrom,
		strconv.Itoa(v.Pos),
		v.ID,
		v.Ref,
		v.Alt,
		qual,
		v.Filter,
		v.BuildInfoString(),
	}

	if v.Format != "" {
		fields = append(fields, v.Format)
		fields = append(fields, v.Samples...)
	}

	return strings.Join(fields, "\t")
}

// Header 表示 VCF 文件头
type Header struct {
	Lines      []string
	InfoFields map[string]string // ID -> Description
	SampleLine string            // #CHROM 行
}

// NewHeader 创建新的 Header
func NewHeader() *Header {
	return &Header{
		Lines:      make([]string, 0),
		InfoFields: make(map[string]string),
	}
}

// AddInfoField 添加 INFO 字段定义
func (h *Header) AddInfoField(id, number, typ, description string) {
	line := fmt.Sprintf("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">", id, number, typ, description)
	h.Lines = append(h.Lines, line)
	h.InfoFields[id] = description
}

// Write 写入 Header 到 writer
func (h *Header) Write(w io.Writer) error {
	for _, line := range h.Lines {
		if _, err := fmt.Fprintln(w, line); err != nil {
			return err
		}
	}
	if h.SampleLine != "" {
		if _, err := fmt.Fprintln(w, h.SampleLine); err != nil {
			return err
		}
	}
	return nil
}

// Parser VCF 解析器
type Parser struct {
	reader   *bufio.Reader
	file     *os.File
	gzReader *gzip.Reader
	Header   *Header
}

// NewParser 创建新的 VCF 解析器
func NewParser(filename string) (*Parser, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("无法打开文件 %s: %w", filename, err)
	}

	var reader io.Reader = file
	var gzReader *gzip.Reader

	// 检测是否为 gzip 文件
	if strings.HasSuffix(filename, ".gz") {
		gzReader, err = gzip.NewReader(file)
		if err != nil {
			file.Close()
			return nil, fmt.Errorf("无法创建 gzip reader: %w", err)
		}
		reader = gzReader
	}

	p := &Parser{
		reader:   bufio.NewReaderSize(reader, 1024*1024), // 1MB 缓冲
		file:     file,
		gzReader: gzReader,
		Header:   NewHeader(),
	}

	// 解析 header
	if err := p.parseHeader(); err != nil {
		p.Close()
		return nil, err
	}

	return p, nil
}

// parseHeader 解析 VCF 头部
func (p *Parser) parseHeader() error {
	for {
		line, err := p.reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return fmt.Errorf("读取 header 失败: %w", err)
		}

		line = strings.TrimRight(line, "\r\n")

		if strings.HasPrefix(line, "##") {
			p.Header.Lines = append(p.Header.Lines, line)
			// 解析 INFO 字段
			if strings.HasPrefix(line, "##INFO=") {
				// 简单提取 ID
				if start := strings.Index(line, "ID="); start != -1 {
					end := strings.Index(line[start:], ",")
					if end != -1 {
						id := line[start+3 : start+end]
						p.Header.InfoFields[id] = line
					}
				}
			}
		} else if strings.HasPrefix(line, "#CHROM") {
			p.Header.SampleLine = line
			break
		}
	}
	return nil
}

// Next 读取下一个变异
func (p *Parser) Next() (*Variant, error) {
	for {
		line, err := p.reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				return nil, io.EOF
			}
			return nil, err
		}

		line = strings.TrimRight(line, "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		return ParseLine(line)
	}
}

// Close 关闭解析器
func (p *Parser) Close() error {
	if p.gzReader != nil {
		p.gzReader.Close()
	}
	return p.file.Close()
}

// ParseLine 解析单行 VCF 记录
func ParseLine(line string) (*Variant, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return nil, fmt.Errorf("VCF 行字段不足: %d", len(fields))
	}

	pos, err := strconv.Atoi(fields[1])
	if err != nil {
		return nil, fmt.Errorf("无法解析位置: %s", fields[1])
	}

	qual := -1.0
	if fields[5] != "." {
		qual, _ = strconv.ParseFloat(fields[5], 64)
	}

	v := &Variant{
		Chrom:   fields[0],
		Pos:     pos,
		ID:      fields[2],
		Ref:     fields[3],
		Alt:     fields[4],
		Qual:    qual,
		Filter:  fields[6],
		Info:    parseInfo(fields[7]),
		RawLine: line,
	}

	if len(fields) > 8 {
		v.Format = fields[8]
		v.Samples = fields[9:]
	}

	return v, nil
}

// parseInfo 解析 INFO 字段
func parseInfo(info string) map[string]string {
	result := make(map[string]string)
	if info == "." || info == "" {
		return result
	}

	parts := strings.Split(info, ";")
	for _, part := range parts {
		if idx := strings.Index(part, "="); idx != -1 {
			key := part[:idx]
			value := part[idx+1:]
			result[key] = value
		} else {
			// Flag 类型
			result[part] = ""
		}
	}
	return result
}

// Writer VCF 写入器
type Writer struct {
	file   *os.File
	writer *bufio.Writer
}

// NewWriter 创建新的 VCF 写入器
func NewWriter(filename string) (*Writer, error) {
	file, err := os.Create(filename)
	if err != nil {
		return nil, fmt.Errorf("无法创建文件 %s: %w", filename, err)
	}

	return &Writer{
		file:   file,
		writer: bufio.NewWriterSize(file, 1024*1024), // 1MB 缓冲
	}, nil
}

// WriteHeader 写入 VCF 头部
func (w *Writer) WriteHeader(header *Header) error {
	return header.Write(w.writer)
}

// WriteVariant 写入变异
func (w *Writer) WriteVariant(v *Variant) error {
	_, err := w.writer.WriteString(v.ToLine() + "\n")
	return err
}

// WriteLine 写入原始行
func (w *Writer) WriteLine(line string) error {
	_, err := w.writer.WriteString(line + "\n")
	return err
}

// Close 关闭写入器
func (w *Writer) Close() error {
	if err := w.writer.Flush(); err != nil {
		w.file.Close()
		return err
	}
	return w.file.Close()
}

// CreateCleanHeader 创建干净的 VCF header
func CreateCleanHeader(sampleLine string) *Header {
	h := NewHeader()

	// 文件格式
	h.Lines = append(h.Lines, "##fileformat=VCFv4.2")

	// INFO 字段定义
	infoFields := []struct {
		ID, Number, Type, Desc string
	}{
		{"gnomAD4.1_joint_POPMAX_AF", "1", "Float", "gnomAD 4.1 joint POPMAX allele frequency"},
		{"gnomAD4.1_joint", "1", "Float", "gnomAD 4.1 joint allele frequency"},
		{"DBSCSNV_ADA", "1", "Float", "dbscSNV ada_score for splice site prediction"},
		{"DBSCSNV_RF", "1", "Float", "dbscSNV rf_score for splice site prediction"},
		{"MetaRNN_score", "1", "Float", "MetaRNN pathogenicity score"},
		{"REVEL_score", "1", "Float", "REVEL pathogenicity score"},
		{"PrimateAI_score", "1", "Float", "PrimateAI pathogenicity score"},
		{"ClinPred_score", "1", "Float", "ClinPred pathogenicity score"},
		{"ESM1b_score", "1", "Float", "ESM1b pathogenicity score (lower = more pathogenic)"},
		{"AlphaMissense_score", "1", "Float", "AlphaMissense pathogenicity score"},
		{"CADD_phred", "1", "Float", "CADD phred-scaled score"},
		{"aaref", "1", "String", "Reference amino acid"},
		{"aaalt", "1", "String", "Alternate amino acid"},
		{"rs_dbSNP", "1", "String", "dbSNP rsID"},
		{"genename", "1", "String", "Gene name"},
		{"clinvar_id", "1", "String", "ClinVar variation ID"},
		{"clinvar_clnsig", "1", "String", "ClinVar clinical significance"},
		{"clinvar_trait", "1", "String", "ClinVar associated trait/disease"},
		{"clinvar_review", "1", "String", "ClinVar review status"},
		{"clinvar_hgvs", "1", "String", "ClinVar HGVS notation"},
		{"clinvar_var_source", "1", "String", "ClinVar variant source"},
		{"clinvar_MedGen_id", "1", "String", "ClinVar MedGen ID"},
		{"clinvar_OMIM_id", "1", "String", "ClinVar OMIM ID"},
		{"clinvar_Orphanet_id", "1", "String", "ClinVar Orphanet ID"},
		{"ACMG_InterVar", "1", "String", "ACMG classification by InterVar-style analysis"},
		{"ACMG_PathEvidence", ".", "String", "Pathogenic evidence codes"},
		{"ACMG_BenEvidence", ".", "String", "Benign evidence codes"},
		{"ACMG_Reason", "1", "String", "Classification reason"},
	}

	for _, f := range infoFields {
		h.AddInfoField(f.ID, f.Number, f.Type, f.Desc)
	}

	h.SampleLine = sampleLine
	return h
}
