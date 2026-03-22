package dbnsfp

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"strings"
)

// DBNSFPRecord 表示一条 dbNSFP 记录
type DBNSFPRecord struct {
	Chrom    string
	Pos      int
	Ref      string
	Alt      string
	Fields   map[string]string
}

// Index 存储 dbNSFP 列索引
type Index struct {
	headers  []string
	colIndex map[string]int
}

// NewIndex 创建列索引
func NewIndex(headerLine string) *Index {
	headers := strings.Split(headerLine, "\t")
	colIndex := make(map[string]int)
	for i, h := range headers {
		colIndex[strings.TrimPrefix(h, "#")] = i
	}
	return &Index{
		headers:  headers,
		colIndex: colIndex,
	}
}

// GetCol 获取列索引
func (idx *Index) GetCol(name string) int {
	if i, ok := idx.colIndex[name]; ok {
		return i
	}
	return -1
}

// Querier 用于查询 dbNSFP 数据库
type Querier struct {
	dbPath string
	index  *Index
}

// NewQuerier 创建查询器
func NewQuerier(dbPath string) (*Querier, error) {
	// 读取 header
	file, err := os.Open(dbPath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	gzReader, err := gzip.NewReader(file)
	if err != nil {
		return nil, err
	}
	defer gzReader.Close()

	reader := bufio.NewReader(gzReader)
	headerLine, err := reader.ReadString('\n')
	if err != nil {
		return nil, err
	}

	return &Querier{
		dbPath: dbPath,
		index:  NewIndex(strings.TrimRight(headerLine, "\r\n")),
	}, nil
}

// QueryByPosition 使用 tabix 查询指定位置
func (q *Querier) QueryByPosition(chrom string, pos int) ([]*DBNSFPRecord, error) {
	// 移除 chr 前缀
	chromNum := chrom
	if strings.HasPrefix(chrom, "chr") {
		chromNum = chrom[3:]
	}

	// 使用 tabix 查询
	region := fmt.Sprintf("%s:%d-%d", chromNum, pos, pos)
	cmd := exec.Command("tabix", q.dbPath, region)
	output, err := cmd.Output()
	if err != nil {
		// tabix 没有结果时返回空
		if exitErr, ok := err.(*exec.ExitError); ok && len(exitErr.Stderr) == 0 {
			return nil, nil
		}
		return nil, err
	}

	if len(output) == 0 {
		return nil, nil
	}

	var records []*DBNSFPRecord
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for _, line := range lines {
		if line == "" {
			continue
		}
		record := q.parseLine(line)
		if record != nil {
			records = append(records, record)
		}
	}

	return records, nil
}

// parseLine 解析 dbNSFP 行
func (q *Querier) parseLine(line string) *DBNSFPRecord {
	fields := strings.Split(line, "\t")
	if len(fields) < 5 {
		return nil
	}

	posCol := q.index.GetCol("pos(1-based)")
	if posCol < 0 {
		posCol = q.index.GetCol("hg38_pos(1-based)")
	}
	if posCol < 0 || posCol >= len(fields) {
		return nil
	}

	pos, err := strconv.Atoi(fields[posCol])
	if err != nil {
		return nil
	}

	record := &DBNSFPRecord{
		Chrom:  fields[q.index.GetCol("chr")],
		Pos:    pos,
		Ref:    fields[q.index.GetCol("ref")],
		Alt:    fields[q.index.GetCol("alt")],
		Fields: make(map[string]string),
	}

	// 提取关键字段
	keyFields := []string{
		"gnomAD_genomes_POPMAX_AF",
		"gnomAD_genomes_AF",
		"clinvar_id",
		"clinvar_clnsig",
		"clinvar_trait",
		"clinvar_review",
		"clinvar_hgvs",
		"clinvar_var_source",
		"clinvar_MedGen_id",
		"clinvar_OMIM_id",
		"clinvar_Orphanet_id",
		"REVEL_score",
		"MetaRNN_score",
		"PrimateAI_score",
		"ClinPred_score",
		"ESM1b_score",
		"AlphaMissense_score",
		"CADD_phred",
		"aaref",
		"aaalt",
		"rs_dbSNP",
		"genename",
	}

	for _, key := range keyFields {
		col := q.index.GetCol(key)
		if col >= 0 && col < len(fields) {
			record.Fields[key] = fields[col]
		}
	}

	return record
}

// BatchQuery 批量查询多个位置
func (q *Querier) BatchQuery(positions []struct{ Chrom string; Pos int; Ref, Alt string }) (map[string]*DBNSFPRecord, error) {
	result := make(map[string]*DBNSFPRecord)

	for _, p := range positions {
		records, err := q.QueryByPosition(p.Chrom, p.Pos)
		if err != nil {
			continue
		}

		for _, r := range records {
			if r.Ref == p.Ref && r.Alt == p.Alt {
				key := fmt.Sprintf("%s_%d_%s_%s", p.Chrom, p.Pos, p.Ref, p.Alt)
				result[key] = r
				break
			}
		}
	}

	return result, nil
}
