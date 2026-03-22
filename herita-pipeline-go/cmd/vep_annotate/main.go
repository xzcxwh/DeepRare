package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"strconv"
	"strings"
	"sync"
	"time"
)

var (
	inputVCF       = flag.String("input", "", "Input VCF file")
	finalVCF       = flag.String("final", "", "Final round VCF file")
	nextVCF        = flag.String("next", "", "Output VCF for next level")
	spliceAIThresh = flag.Float64("spliceai-thresh", 0.5, "SpliceAI delta score threshold")
	batchSize      = flag.Int("batch", 50, "Batch size for VEP API requests")
	workers        = flag.Int("workers", 4, "Number of parallel workers")
	retries        = flag.Int("retries", 3, "Number of retries for failed requests")
)

const (
	vepAPIURL = "https://rest.ensembl.org/vep/human/region"
)

// VEP response structures
type VEPResponse struct {
	Input                   string                  `json:"input"`
	MostSevereConsequence   string                  `json:"most_severe_consequence"`
	TranscriptConsequences  []TranscriptConsequence `json:"transcript_consequences"`
	ColocatedVariants       []ColocatedVariant      `json:"colocated_variants"`
}

type TranscriptConsequence struct {
	GeneSymbol      string   `json:"gene_symbol"`
	Impact          string   `json:"impact"`
	Consequence     []string `json:"consequence_terms"`
	AlphaMissense   *AMScore `json:"alphamissense,omitempty"`
	CADD            *float64 `json:"cadd_phred,omitempty"`
	CADDRaw         *float64 `json:"cadd_raw,omitempty"`
	REVEL           *float64 `json:"revel,omitempty"`
	SpliceAI        *SpliceAIScore `json:"spliceai,omitempty"`
	Conservation    *float64 `json:"conservation,omitempty"`
	Canonical       *int     `json:"canonical,omitempty"`
}

type AMScore struct {
	Class         string  `json:"am_class"`
	Pathogenicity float64 `json:"am_pathogenicity"`
}

type SpliceAIScore struct {
	DeltaScoreAG float64 `json:"DS_AG"`
	DeltaScoreAL float64 `json:"DS_AL"`
	DeltaScoreDG float64 `json:"DS_DG"`
	DeltaScoreDL float64 `json:"DS_DL"`
	Symbol       string  `json:"SYMBOL,omitempty"`
}

type ColocatedVariant struct {
	ID          string                 `json:"id"`
	Frequencies map[string]interface{} `json:"frequencies"`
}

// Variant info
type Variant struct {
	Chrom      string
	Pos        int
	Ref        string
	Alt        string
	OrigLine   string
	VEPResult  *VEPAnnotation
}

type VEPAnnotation struct {
	AlphaMissense    string
	AMPathogenicity  string
	SpliceAIMax      float64
	SpliceAIDetails  string
	CADD             string
	REVEL            string
	PhyloP           string
	Gene             string
	Consequence      string
}

func main() {
	flag.Parse()

	if *inputVCF == "" || *finalVCF == "" || *nextVCF == "" {
		log.Fatal("Usage: vep_annotate -input <input.vcf> -final <final.vcf> -next <next.vcf>")
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA VEP REST API Annotation Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Final Round VCF: %s\n", *finalVCF)
	fmt.Printf("Next Level VCF: %s\n", *nextVCF)
	fmt.Printf("SpliceAI Threshold: %.2f\n", *spliceAIThresh)
	fmt.Printf("Batch Size: %d\n", *batchSize)
	fmt.Printf("Workers: %d\n", *workers)
	fmt.Println("============================================================")
	fmt.Println("Annotating fields: AlphaMissense, SpliceAI, CADD, REVEL, phyloP")
	fmt.Println("============================================================")

	// Read variants from VCF
	fmt.Println("Step 1: Reading input VCF...")
	variants, headerLines, err := readVCF(*inputVCF)
	if err != nil {
		log.Fatalf("Failed to read VCF: %v", err)
	}
	fmt.Printf("Read %d variants\n", len(variants))

	// Annotate variants using VEP API
	fmt.Println("Step 2: Annotating variants using VEP REST API...")
	annotateVariants(variants)

	// Filter and write results
	fmt.Println("Step 3: Filtering and writing results...")
	stats := writeResults(variants, headerLines, *finalVCF, *nextVCF, *spliceAIThresh)

	elapsed := time.Since(startTime)
	fmt.Println("============================================================")
	fmt.Println("VEP Annotation Complete!")
	fmt.Printf("Total Time: %s\n", elapsed.Round(time.Second))
	fmt.Println("============================================================")
	fmt.Printf("Input Variants: %d\n", len(variants))
	fmt.Printf("Successfully Annotated: %d\n", stats.annotated)
	fmt.Printf("With SpliceAI scores: %d\n", stats.withSpliceAI)
	fmt.Println("------------------------------------------------------------")
	fmt.Printf("→ Final Round (SpliceAI > %.1f): %d\n", *spliceAIThresh, stats.toFinal)
	fmt.Printf("→ Next Level: %d\n", stats.toNext)
	fmt.Println("============================================================")

	if len(stats.finalVariants) > 0 {
		fmt.Println("\nHigh SpliceAI Variants Added to Final Round:")
		for i, v := range stats.finalVariants {
			fmt.Printf("  %d. %s\n", i+1, v)
		}
		fmt.Println("============================================================")
	}
}

func readVCF(filename string) ([]*Variant, []string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, err
	}
	defer file.Close()

	var variants []*Variant
	var headerLines []string

	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			headerLines = append(headerLines, line)
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 5 {
			continue
		}

		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}

		variants = append(variants, &Variant{
			Chrom:    fields[0],
			Pos:      pos,
			Ref:      fields[3],
			Alt:      fields[4],
			OrigLine: line,
		})
	}

	return variants, headerLines, scanner.Err()
}

func annotateVariants(variants []*Variant) {
	// Create batches
	var batches [][]*Variant
	for i := 0; i < len(variants); i += *batchSize {
		end := i + *batchSize
		if end > len(variants) {
			end = len(variants)
		}
		batches = append(batches, variants[i:end])
	}

	fmt.Printf("Processing %d batches with %d workers...\n", len(batches), *workers)

	// Process batches with worker pool
	var wg sync.WaitGroup
	semaphore := make(chan struct{}, *workers)
	progress := 0
	var mu sync.Mutex

	for i, batch := range batches {
		wg.Add(1)
		go func(batchNum int, batchVariants []*Variant) {
			defer wg.Done()
			semaphore <- struct{}{}
			defer func() { <-semaphore }()

			processBatch(batchVariants)

			mu.Lock()
			progress++
			if progress%10 == 0 || progress == len(batches) {
				fmt.Printf("  Progress: %d/%d batches (%.1f%%)\n", progress, len(batches), float64(progress)/float64(len(batches))*100)
			}
			mu.Unlock()
		}(i, batch)
	}

	wg.Wait()
}

func processBatch(variants []*Variant) {
	// Build request body using VCF format strings
	var vcfLines []string
	variantMap := make(map[string]*Variant)

	for _, v := range variants {
		// Remove chr prefix for Ensembl API
		chrom := strings.TrimPrefix(v.Chrom, "chr")
		
		// Create VCF format string: "chrom pos . ref alt . . ."
		vcfLine := fmt.Sprintf("%s %d . %s %s . . .", chrom, v.Pos, v.Ref, v.Alt)
		
		vcfLines = append(vcfLines, vcfLine)
		variantMap[vcfLine] = v
	}

	// Make API request with retries
	var results []VEPResponse
	var err error
	
	for attempt := 0; attempt < *retries; attempt++ {
		results, err = callVEPAPI(vcfLines)
		if err == nil {
			break
		}
		time.Sleep(time.Duration(attempt+1) * time.Second)
	}

	if err != nil {
		log.Printf("Warning: Failed to annotate batch after %d retries: %v", *retries, err)
		return
	}

	// Map results back to variants
	for _, result := range results {
		// Find matching variant by input string
		if v, ok := variantMap[result.Input]; ok {
			v.VEPResult = parseVEPResult(&result)
		}
	}
}

func callVEPAPI(vcfLines []string) ([]VEPResponse, error) {
	// Build request body with VCF format variants and plugin parameters
	requestBody := map[string]interface{}{
		"variants":      vcfLines,
		"SpliceAI":      1,           // Enable SpliceAI predictions
		"CADD":          1,           // Enable CADD scores
		"REVEL":         1,           // Enable REVEL scores
		"AlphaMissense": 1,           // Enable AlphaMissense predictions
		"Conservation":  1,           // Enable conservation scores (phyloP)
		"canonical":     1,           // Include canonical transcript flag
		"hgvs":          1,           // Include HGVS notation
	}

	jsonBody, err := json.Marshal(requestBody)
	if err != nil {
		return nil, err
	}

	// Make POST request
	req, err := http.NewRequest("POST", vepAPIURL, bytes.NewBuffer(jsonBody))
	if err != nil {
		return nil, err
	}

	req.Header.Set("Content-Type", "application/json")
	req.Header.Set("Accept", "application/json")

	client := &http.Client{Timeout: 60 * time.Second}
	resp, err := client.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != 200 {
		body, _ := io.ReadAll(resp.Body)
		return nil, fmt.Errorf("API error %d: %s", resp.StatusCode, string(body))
	}

	// Parse response
	var results []VEPResponse
	if err := json.NewDecoder(resp.Body).Decode(&results); err != nil {
		return nil, err
	}

	return results, nil
}

func parseVEPResult(result *VEPResponse) *VEPAnnotation {
	ann := &VEPAnnotation{
		Consequence: result.MostSevereConsequence,
	}

	// Find best transcript annotation
	var bestCADD float64 = -1
	var bestREVEL float64 = -1
	var bestSpliceAI float64 = -1
	var spliceAIDetails string
	var bestConservation float64 = -999

	for _, tc := range result.TranscriptConsequences {
		if ann.Gene == "" && tc.GeneSymbol != "" {
			ann.Gene = tc.GeneSymbol
		}

		// AlphaMissense
		if tc.AlphaMissense != nil {
			ann.AlphaMissense = fmt.Sprintf("%.4f", tc.AlphaMissense.Pathogenicity)
			ann.AMPathogenicity = tc.AlphaMissense.Class
		}

		// CADD
		if tc.CADD != nil && *tc.CADD > bestCADD {
			bestCADD = *tc.CADD
		}

		// REVEL
		if tc.REVEL != nil && *tc.REVEL > bestREVEL {
			bestREVEL = *tc.REVEL
		}

		// SpliceAI - get max delta score
		if tc.SpliceAI != nil {
			maxDelta := max(
				tc.SpliceAI.DeltaScoreAG,
				tc.SpliceAI.DeltaScoreAL,
				tc.SpliceAI.DeltaScoreDG,
				tc.SpliceAI.DeltaScoreDL,
			)
			if maxDelta > bestSpliceAI {
				bestSpliceAI = maxDelta
				spliceAIDetails = fmt.Sprintf("AG=%.2f,AL=%.2f,DG=%.2f,DL=%.2f",
					tc.SpliceAI.DeltaScoreAG, tc.SpliceAI.DeltaScoreAL,
					tc.SpliceAI.DeltaScoreDG, tc.SpliceAI.DeltaScoreDL)
			}
		}

		// Conservation (phyloP)
		if tc.Conservation != nil && *tc.Conservation > bestConservation {
			bestConservation = *tc.Conservation
		}
	}

	if bestCADD >= 0 {
		ann.CADD = fmt.Sprintf("%.2f", bestCADD)
	}
	if bestREVEL >= 0 {
		ann.REVEL = fmt.Sprintf("%.4f", bestREVEL)
	}
	if bestSpliceAI >= 0 {
		ann.SpliceAIMax = bestSpliceAI
		ann.SpliceAIDetails = spliceAIDetails
	}
	if bestConservation > -999 {
		ann.PhyloP = fmt.Sprintf("%.3f", bestConservation)
	}

	return ann
}

func max(values ...float64) float64 {
	m := values[0]
	for _, v := range values[1:] {
		if v > m {
			m = v
		}
	}
	return m
}

type writeStats struct {
	annotated     int
	withSpliceAI  int
	toFinal       int
	toNext        int
	finalVariants []string
}

func writeResults(variants []*Variant, headerLines []string, finalVCF, nextVCF string, spliceAIThresh float64) writeStats {
	stats := writeStats{}

	// Add VEP annotation headers
	newHeaders := []string{
		"##INFO=<ID=VEP_Gene,Number=.,Type=String,Description=\"Gene symbol from VEP\">",
		"##INFO=<ID=VEP_Consequence,Number=.,Type=String,Description=\"Most severe consequence from VEP\">",
		"##INFO=<ID=VEP_AlphaMissense,Number=.,Type=String,Description=\"AlphaMissense score\">",
		"##INFO=<ID=VEP_AM_Pathogenicity,Number=.,Type=String,Description=\"AlphaMissense pathogenicity\">",
		"##INFO=<ID=VEP_SpliceAI_Max,Number=.,Type=Float,Description=\"Maximum SpliceAI delta score\">",
		"##INFO=<ID=VEP_SpliceAI_Details,Number=.,Type=String,Description=\"SpliceAI delta scores (AG,AL,DG,DL)\">",
		"##INFO=<ID=VEP_CADD,Number=.,Type=Float,Description=\"CADD phred score\">",
		"##INFO=<ID=VEP_REVEL,Number=.,Type=Float,Description=\"REVEL score\">",
		"##INFO=<ID=VEP_PhyloP,Number=.,Type=Float,Description=\"phyloP100way conservation score\">",
	}

	// Open final VCF in append mode
	finalFile, err := os.OpenFile(finalVCF, os.O_APPEND|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Failed to open final VCF: %v", err)
	}
	defer finalFile.Close()
	finalWriter := bufio.NewWriter(finalFile)
	defer finalWriter.Flush()

	// Create next level VCF
	nextFile, err := os.Create(nextVCF)
	if err != nil {
		log.Fatalf("Failed to create next VCF: %v", err)
	}
	defer nextFile.Close()
	nextWriter := bufio.NewWriter(nextFile)
	defer nextWriter.Flush()

	// Write headers to next VCF
	for _, h := range headerLines {
		if strings.HasPrefix(h, "#CHROM") {
			// Insert new headers before #CHROM
			for _, nh := range newHeaders {
				nextWriter.WriteString(nh + "\n")
			}
		}
		nextWriter.WriteString(h + "\n")
	}

	// Process variants
	for _, v := range variants {
		// Add VEP annotations to INFO field
		annotatedLine := addVEPAnnotations(v)

		if v.VEPResult != nil {
			stats.annotated++
			if v.VEPResult.SpliceAIMax > 0 {
				stats.withSpliceAI++
			}
		}

		// Check SpliceAI threshold
		if v.VEPResult != nil && v.VEPResult.SpliceAIMax > spliceAIThresh {
			finalWriter.WriteString(annotatedLine + "\n")
			stats.toFinal++
			stats.finalVariants = append(stats.finalVariants,
				fmt.Sprintf("%s:%d %s>%s [SpliceAI=%.2f, Gene=%s]",
					v.Chrom, v.Pos, v.Ref, v.Alt, v.VEPResult.SpliceAIMax, v.VEPResult.Gene))
		} else {
			nextWriter.WriteString(annotatedLine + "\n")
			stats.toNext++
		}
	}

	return stats
}

func addVEPAnnotations(v *Variant) string {
	if v.VEPResult == nil {
		return v.OrigLine
	}

	fields := strings.Split(v.OrigLine, "\t")
	if len(fields) < 8 {
		return v.OrigLine
	}

	// Build annotation string
	var annParts []string

	if v.VEPResult.Gene != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_Gene=%s", v.VEPResult.Gene))
	}
	if v.VEPResult.Consequence != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_Consequence=%s", v.VEPResult.Consequence))
	}
	if v.VEPResult.AlphaMissense != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_AlphaMissense=%s", v.VEPResult.AlphaMissense))
	}
	if v.VEPResult.AMPathogenicity != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_AM_Pathogenicity=%s", v.VEPResult.AMPathogenicity))
	}
	if v.VEPResult.SpliceAIMax > 0 {
		annParts = append(annParts, fmt.Sprintf("VEP_SpliceAI_Max=%.4f", v.VEPResult.SpliceAIMax))
	}
	if v.VEPResult.SpliceAIDetails != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_SpliceAI_Details=%s", v.VEPResult.SpliceAIDetails))
	}
	if v.VEPResult.CADD != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_CADD=%s", v.VEPResult.CADD))
	}
	if v.VEPResult.REVEL != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_REVEL=%s", v.VEPResult.REVEL))
	}
	if v.VEPResult.PhyloP != "" {
		annParts = append(annParts, fmt.Sprintf("VEP_PhyloP=%s", v.VEPResult.PhyloP))
	}

	if len(annParts) > 0 {
		fields[7] = fields[7] + ";" + strings.Join(annParts, ";")
	}

	return strings.Join(fields, "\t")
}
