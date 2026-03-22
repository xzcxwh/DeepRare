package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"
)

var (
	inputVCF   = flag.String("input", "", "Input VCF file")
	finalVCF   = flag.String("final", "", "Final round VCF file (HIGH impact variants will be appended)")
	nextVCF    = flag.String("next", "", "Output VCF for variants needing further evaluation")
	snpEffPath = flag.String("snpeff", "/Volumes/T9/snpEff/snpEff/snpEff.jar", "Path to snpEff.jar")
	javaPath   = flag.String("java", "/opt/homebrew/opt/java/bin/java", "Path to Java executable (requires Java 21+)")
	genome     = flag.String("genome", "GRCh38.99", "Genome version for snpEff")
	memory     = flag.String("memory", "4g", "Java heap memory")
)

func main() {
	flag.Parse()

	if *inputVCF == "" || *finalVCF == "" || *nextVCF == "" {
		log.Fatal("Usage: snpeff_filter -input <input.vcf> -final <final.vcf> -next <next.vcf>")
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA snpEff Annotation & HIGH Impact Filter Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Final Round VCF: %s\n", *finalVCF)
	fmt.Printf("Next Level VCF: %s\n", *nextVCF)
	fmt.Printf("snpEff: %s\n", *snpEffPath)
	fmt.Printf("Genome: %s\n", *genome)
	fmt.Println("============================================================")

	// Verify snpEff exists
	if _, err := os.Stat(*snpEffPath); os.IsNotExist(err) {
		log.Fatalf("snpEff.jar not found at: %s", *snpEffPath)
	}

	// Verify Java exists
	if _, err := os.Stat(*javaPath); os.IsNotExist(err) {
		log.Fatalf("Java not found at: %s", *javaPath)
	}

	// Create temp directory for intermediate files
	tempDir, err := os.MkdirTemp("", "snpeff_")
	if err != nil {
		log.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir)

	// Step 1: Run snpEff annotation
	fmt.Println("Step 1: Running snpEff annotation...")
	annotatedVCF := filepath.Join(tempDir, "annotated.vcf")
	err = runSnpEff(*inputVCF, annotatedVCF)
	if err != nil {
		log.Fatalf("snpEff annotation failed: %v", err)
	}
	fmt.Println("✓ snpEff annotation complete")

	// Step 2: Parse and filter by impact
	fmt.Println("Step 2: Filtering by impact level...")
	stats, highVariants, err := filterByImpact(annotatedVCF, *finalVCF, *nextVCF)
	if err != nil {
		log.Fatalf("Filtering failed: %v", err)
	}

	elapsed := time.Since(startTime)

	fmt.Println("============================================================")
	fmt.Println("snpEff Annotation & Filtering Complete!")
	fmt.Printf("Total Time: %s\n", elapsed.Round(time.Millisecond))
	fmt.Println("============================================================")
	fmt.Printf("Input Variants: %d\n", stats.total)
	fmt.Println("------------------------------------------------------------")
	fmt.Printf("→ Final Round (HIGH impact): %d\n", stats.high)
	fmt.Printf("→ Next Level (MODERATE/LOW/MODIFIER): %d\n", stats.other)
	fmt.Println("============================================================")

	if len(highVariants) > 0 {
		fmt.Println("\nHIGH Impact Variants Added to Final Round:")
		for i, v := range highVariants {
			fmt.Printf("  %d. %s\n", i+1, v)
		}
		fmt.Println("============================================================")
	}
}

func runSnpEff(inputVCF, outputVCF string) error {
	// Get snpEff directory for config
	snpEffDir := filepath.Dir(*snpEffPath)

	// Convert to absolute paths
	absInputVCF, err := filepath.Abs(inputVCF)
	if err != nil {
		return fmt.Errorf("failed to get absolute path for input: %v", err)
	}
	absOutputVCF, err := filepath.Abs(outputVCF)
	if err != nil {
		return fmt.Errorf("failed to get absolute path for output: %v", err)
	}

	args := []string{
		fmt.Sprintf("-Xmx%s", *memory),
		"-jar", *snpEffPath,
		"-v",
		"-noStats",
		"-noLog",
		"-config", filepath.Join(snpEffDir, "snpEff.config"),
		"-dataDir", filepath.Join(snpEffDir, "data"),
		*genome,
		absInputVCF,
	}

	cmd := exec.Command(*javaPath, args...)
	cmd.Dir = snpEffDir

	outFile, err := os.Create(absOutputVCF)
	if err != nil {
		return err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr

	return cmd.Run()
}

type filterStats struct {
	total int
	high  int
	other int
}

func filterByImpact(annotatedVCF, finalVCF, nextVCF string) (filterStats, []string, error) {
	stats := filterStats{}
	var highVariants []string

	// Open annotated VCF
	inFile, err := os.Open(annotatedVCF)
	if err != nil {
		return stats, nil, err
	}
	defer inFile.Close()

	// Open final VCF in append mode
	finalFile, err := os.OpenFile(finalVCF, os.O_APPEND|os.O_WRONLY, 0644)
	if err != nil {
		return stats, nil, err
	}
	defer finalFile.Close()
	finalWriter := bufio.NewWriter(finalFile)
	defer finalWriter.Flush()

	// Create next level VCF
	nextFile, err := os.Create(nextVCF)
	if err != nil {
		return stats, nil, err
	}
	defer nextFile.Close()
	nextWriter := bufio.NewWriter(nextFile)
	defer nextWriter.Flush()

	scanner := bufio.NewScanner(inFile)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()

		// Write headers only to next VCF (final already has headers)
		if strings.HasPrefix(line, "#") {
			nextWriter.WriteString(line + "\n")
			continue
		}

		stats.total++

		// Check if HIGH impact
		if isHighImpact(line) {
			finalWriter.WriteString(line + "\n")
			stats.high++

			// Extract variant info for reporting
			fields := strings.Split(line, "\t")
			if len(fields) >= 5 {
				gene := extractGene(fields[7])
				effect := extractEffect(fields[7])
				highVariants = append(highVariants, 
					fmt.Sprintf("%s:%s %s>%s [%s, %s]", 
						fields[0], fields[1], fields[3], fields[4], gene, effect))
			}
		} else {
			nextWriter.WriteString(line + "\n")
			stats.other++
		}
	}

	return stats, highVariants, scanner.Err()
}

func isHighImpact(line string) bool {
	// Look for HIGH impact in the ANN field
	// Format: ANN=A|effect|HIGH|gene|...
	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return false
	}

	info := fields[7]
	
	// Find ANN field
	for _, part := range strings.Split(info, ";") {
		if strings.HasPrefix(part, "ANN=") {
			// Check each annotation
			annotations := strings.TrimPrefix(part, "ANN=")
			for _, ann := range strings.Split(annotations, ",") {
				annFields := strings.Split(ann, "|")
				if len(annFields) >= 3 {
					impact := annFields[2]
					if impact == "HIGH" {
						return true
					}
				}
			}
		}
	}
	return false
}

func extractGene(info string) string {
	for _, part := range strings.Split(info, ";") {
		if strings.HasPrefix(part, "ANN=") {
			annotations := strings.TrimPrefix(part, "ANN=")
			annFields := strings.Split(annotations, "|")
			if len(annFields) >= 4 {
				return annFields[3] // Gene name is at position 3
			}
		}
	}
	return "unknown"
}

func extractEffect(info string) string {
	for _, part := range strings.Split(info, ";") {
		if strings.HasPrefix(part, "ANN=") {
			annotations := strings.TrimPrefix(part, "ANN=")
			annFields := strings.Split(annotations, "|")
			if len(annFields) >= 2 {
				return annFields[1] // Effect is at position 1
			}
		}
	}
	return "unknown"
}
