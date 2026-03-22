package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
	"time"
)

var (
	inputVCF  = flag.String("input", "", "Input VCF file with ClinVar annotations")
	finalVCF  = flag.String("final", "", "Output VCF for variants that pass to final round (pathogenic with high confidence)")
	nextVCF   = flag.String("next", "", "Output VCF for variants that need further evaluation")
)

// 高置信度的 CLNREVSTAT 值
var highConfidenceRevstat = map[string]bool{
	"practice_guideline":                              true,
	"reviewed_by_expert_panel":                        true,
	"criteria_provided,_multiple_submitters,_no_conflicts": true,
}

// 致病性 CLNSIG 值 (进入决赛圈)
var pathogenicSig = map[string]bool{
	"Pathogenic":                  true,
	"Likely_pathogenic":           true,
	"Pathogenic/Likely_pathogenic": true,
}

// 良性 CLNSIG 值 (直接删除)
var benignSig = map[string]bool{
	"Benign":              true,
	"Likely_benign":       true,
	"Benign/Likely_benign": true,
}

func main() {
	flag.Parse()

	if *inputVCF == "" || *finalVCF == "" || *nextVCF == "" {
		log.Fatal("Usage: clinvar_filter -input <input.vcf> -final <final.vcf> -next <next.vcf>")
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA ClinVar Filter Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Final Round Output: %s\n", *finalVCF)
	fmt.Printf("Next Level Output: %s\n", *nextVCF)
	fmt.Println("============================================================")
	fmt.Println("Filter Logic:")
	fmt.Println("  1. Pathogenic/Likely_pathogenic + High Confidence → Final Round")
	fmt.Println("  2. Benign/Likely_benign + High Confidence → Removed")
	fmt.Println("  3. All others → Next Level")
	fmt.Println("============================================================")

	// Open input file
	inFile, err := os.Open(*inputVCF)
	if err != nil {
		log.Fatalf("Failed to open input file: %v", err)
	}
	defer inFile.Close()

	// Create output files
	finalFile, err := os.Create(*finalVCF)
	if err != nil {
		log.Fatalf("Failed to create final output file: %v", err)
	}
	defer finalFile.Close()

	nextFile, err := os.Create(*nextVCF)
	if err != nil {
		log.Fatalf("Failed to create next output file: %v", err)
	}
	defer nextFile.Close()

	finalWriter := bufio.NewWriter(finalFile)
	nextWriter := bufio.NewWriter(nextFile)
	defer finalWriter.Flush()
	defer nextWriter.Flush()

	scanner := bufio.NewScanner(inFile)
	// Increase buffer for long lines
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	// Statistics
	stats := struct {
		total        int
		toFinal      int
		removed      int
		toNext       int
		noAnnotation int
	}{}

	// Track variants going to final for reporting
	var finalVariants []string

	for scanner.Scan() {
		line := scanner.Text()

		// Write headers to both files
		if strings.HasPrefix(line, "#") {
			finalWriter.WriteString(line + "\n")
			nextWriter.WriteString(line + "\n")
			continue
		}

		stats.total++

		// Parse INFO field
		fields := strings.Split(line, "\t")
		if len(fields) < 8 {
			continue
		}

		info := fields[7]
		clnsig := extractField(info, "ClinVar_CLNSIG")
		clnrevstat := extractField(info, "ClinVar_CLNREVSTAT")

		// Determine category
		category := categorize(clnsig, clnrevstat)

		switch category {
		case "final":
			finalWriter.WriteString(line + "\n")
			stats.toFinal++
			// Record for reporting
			chrom := fields[0]
			pos := fields[1]
			ref := fields[3]
			alt := fields[4]
			finalVariants = append(finalVariants, fmt.Sprintf("%s:%s %s>%s [%s, %s]", chrom, pos, ref, alt, clnsig, clnrevstat))
		case "removed":
			stats.removed++
			// Don't write anywhere
		case "next":
			nextWriter.WriteString(line + "\n")
			stats.toNext++
		case "no_annotation":
			nextWriter.WriteString(line + "\n")
			stats.noAnnotation++
			stats.toNext++
		}
	}

	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading input: %v", err)
	}

	elapsed := time.Since(startTime)

	fmt.Println("============================================================")
	fmt.Println("Filtering Complete!")
	fmt.Printf("Total Time: %s\n", elapsed.Round(time.Millisecond))
	fmt.Println("============================================================")
	fmt.Printf("Input Variants: %d\n", stats.total)
	fmt.Println("------------------------------------------------------------")
	fmt.Printf("→ Final Round (Pathogenic + High Confidence): %d\n", stats.toFinal)
	fmt.Printf("✗ Removed (Benign + High Confidence): %d\n", stats.removed)
	fmt.Printf("→ Next Level (needs further evaluation): %d\n", stats.toNext)
	fmt.Printf("  (including %d without ClinVar annotation)\n", stats.noAnnotation)
	fmt.Println("============================================================")

	if len(finalVariants) > 0 {
		fmt.Println("\nVariants in Final Round:")
		for i, v := range finalVariants {
			fmt.Printf("  %d. %s\n", i+1, v)
		}
		fmt.Println("============================================================")
	}
}

func extractField(info, fieldName string) string {
	prefix := fieldName + "="
	for _, part := range strings.Split(info, ";") {
		if strings.HasPrefix(part, prefix) {
			return strings.TrimPrefix(part, prefix)
		}
	}
	return ""
}

func categorize(clnsig, clnrevstat string) string {
	// No ClinVar annotation
	if clnsig == "" {
		return "no_annotation"
	}

	// Check if high confidence review status
	isHighConfidence := highConfidenceRevstat[clnrevstat]

	// Rule 1: Pathogenic + High Confidence → Final
	if pathogenicSig[clnsig] && isHighConfidence {
		return "final"
	}

	// Rule 2: Benign + High Confidence → Removed
	if benignSig[clnsig] && isHighConfidence {
		return "removed"
	}

	// Rule 3: Everything else → Next Level
	return "next"
}
