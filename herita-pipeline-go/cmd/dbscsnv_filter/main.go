package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"
)

var (
	inputVCF    = flag.String("input", "", "Input VCF file")
	finalVCF    = flag.String("final", "", "Final round VCF file (high scoring variants will be appended)")
	nextVCF     = flag.String("next", "", "Output VCF for variants needing further evaluation")
	dbscSNVFile = flag.String("dbscsnv", "/Volumes/T9/herita-project/dbscSNV1/dbscSNV1.1_hg38_sorted.tsv.gz", "Path to dbscSNV TSV file (bgzipped)")
	adaThresh   = flag.Float64("ada-thresh", 0.6, "Threshold for ada_score")
	rfThresh    = flag.Float64("rf-thresh", 0.6, "Threshold for rf_score")
)

func main() {
	flag.Parse()

	if *inputVCF == "" || *finalVCF == "" || *nextVCF == "" {
		log.Fatal("Usage: dbscsnv_filter -input <input.vcf> -final <final.vcf> -next <next.vcf>")
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA dbscSNV Annotation & Filter Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Final Round VCF: %s\n", *finalVCF)
	fmt.Printf("Next Level VCF: %s\n", *nextVCF)
	fmt.Printf("dbscSNV Data: %s\n", *dbscSNVFile)
	fmt.Printf("Thresholds: ada_score > %.2f AND rf_score > %.2f\n", *adaThresh, *rfThresh)
	fmt.Println("============================================================")

	// Check dependencies
	checkDependency("bcftools")
	checkDependency("vcfanno")
	checkDependency("bgzip")
	checkDependency("tabix")

	// Create temp directory
	tempDir, err := os.MkdirTemp("", "dbscsnv_")
	if err != nil {
		log.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir)

	// Step 1: Prepare dbscSNV data for vcfanno (BED format)
	fmt.Println("Step 1: Preparing dbscSNV annotation data...")
	bedFile := filepath.Join(tempDir, "dbscSNV.bed.gz")
	err = prepareDbscSNVBed(*dbscSNVFile, bedFile)
	if err != nil {
		log.Fatalf("Failed to prepare dbscSNV data: %v", err)
	}
	fmt.Println("✓ dbscSNV data prepared")

	// Step 2: Create vcfanno config
	fmt.Println("Step 2: Creating vcfanno config...")
	tomlFile := filepath.Join(tempDir, "dbscsnv_config.toml")
	createVcfannoConfig(tomlFile, bedFile)

	// Step 3: Prepare input VCF
	fmt.Println("Step 3: Preparing input VCF...")
	workingVCF, err := prepareInputVCF(*inputVCF, tempDir)
	if err != nil {
		log.Fatalf("Failed to prepare input VCF: %v", err)
	}

	// Step 4: Run vcfanno
	fmt.Println("Step 4: Running vcfanno annotation...")
	annotatedVCF := filepath.Join(tempDir, "annotated.vcf")
	err = runVcfanno(tomlFile, workingVCF, annotatedVCF)
	if err != nil {
		log.Fatalf("vcfanno annotation failed: %v", err)
	}
	fmt.Println("✓ vcfanno annotation complete")

	// Step 5: Filter by scores
	fmt.Println("Step 5: Filtering by dbscSNV scores...")
	stats, highScoringVariants, err := filterByScores(annotatedVCF, *finalVCF, *nextVCF, *adaThresh, *rfThresh)
	if err != nil {
		log.Fatalf("Filtering failed: %v", err)
	}

	elapsed := time.Since(startTime)

	fmt.Println("============================================================")
	fmt.Println("dbscSNV Annotation & Filtering Complete!")
	fmt.Printf("Total Time: %s\n", elapsed.Round(time.Millisecond))
	fmt.Println("============================================================")
	fmt.Printf("Input Variants: %d\n", stats.total)
	fmt.Printf("With dbscSNV annotation: %d\n", stats.annotated)
	fmt.Println("------------------------------------------------------------")
	fmt.Printf("→ Final Round (ada>%.1f AND rf>%.1f): %d\n", *adaThresh, *rfThresh, stats.high)
	fmt.Printf("→ Next Level: %d\n", stats.other)
	fmt.Println("============================================================")

	if len(highScoringVariants) > 0 {
		fmt.Println("\nHigh-scoring Splice Variants Added to Final Round:")
		for i, v := range highScoringVariants {
			fmt.Printf("  %d. %s\n", i+1, v)
		}
		fmt.Println("============================================================")
	}
}

func checkDependency(tool string) {
	_, err := exec.LookPath(tool)
	if err != nil {
		log.Fatalf("Error: %s is not installed or not in PATH", tool)
	}
}

func prepareDbscSNVBed(tsvFile, bedFile string) error {
	// Read TSV and convert to BED format for vcfanno
	// BED format: chr, start (0-based), end, ref, alt, ada_score, rf_score
	
	// Create temp uncompressed file first
	tempBed := strings.TrimSuffix(bedFile, ".gz")
	outFile, err := os.Create(tempBed)
	if err != nil {
		return err
	}
	writer := bufio.NewWriter(outFile)

	var cmd *exec.Cmd
	if strings.HasSuffix(tsvFile, ".gz") {
		cmd = exec.Command("gzcat", tsvFile)
	} else {
		cmd = exec.Command("cat", tsvFile)
	}
	
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		outFile.Close()
		return err
	}
	
	if err := cmd.Start(); err != nil {
		outFile.Close()
		return err
	}

	scanner := bufio.NewScanner(stdout)
	// Increase buffer size
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()
		
		// Skip header
		if lineNum == 1 {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 6 {
			continue
		}

		chr := fields[0]
		pos := fields[1]
		ref := fields[2]
		alt := fields[3]
		adaScore := fields[4]
		rfScore := fields[5]

		// Add chr prefix if missing
		if !strings.HasPrefix(chr, "chr") {
			chr = "chr" + chr
		}

		// Convert to 0-based start for BED
		posInt, err := strconv.Atoi(pos)
		if err != nil {
			continue
		}
		start := posInt - 1
		end := posInt

		// Write BED line: chr, start, end, ref, alt, ada_score, rf_score
		fmt.Fprintf(writer, "%s\t%d\t%d\t%s\t%s\t%s\t%s\n", chr, start, end, ref, alt, adaScore, rfScore)
	}

	writer.Flush()
	outFile.Close()
	cmd.Wait()

	// Sort with chromosome order
	sortedBed := tempBed + ".sorted"
	sortCmd := exec.Command("sort", "-k1,1V", "-k2,2n", tempBed)
	sortedFile, err := os.Create(sortedBed)
	if err != nil {
		return err
	}
	sortCmd.Stdout = sortedFile
	if err := sortCmd.Run(); err != nil {
		sortedFile.Close()
		return fmt.Errorf("sort failed: %v", err)
	}
	sortedFile.Close()

	// bgzip
	bgzipCmd := exec.Command("bgzip", "-f", sortedBed)
	if err := bgzipCmd.Run(); err != nil {
		return fmt.Errorf("bgzip failed: %v", err)
	}

	// Rename to target
	os.Rename(sortedBed+".gz", bedFile)

	// Create tabix index
	tabixCmd := exec.Command("tabix", "-f", "-p", "bed", bedFile)
	if err := tabixCmd.Run(); err != nil {
		return fmt.Errorf("tabix failed: %v", err)
	}

	return nil
}

func createVcfannoConfig(tomlPath, bedFile string) {
	config := fmt.Sprintf(`[[annotation]]
file = "%s"
columns = [4, 5, 6, 7]
names = ["dbscSNV_ref", "dbscSNV_alt", "dbscSNV_ada_score", "dbscSNV_rf_score"]
ops = ["self", "self", "self", "self"]
`, bedFile)

	if err := os.WriteFile(tomlPath, []byte(config), 0644); err != nil {
		log.Fatalf("Failed to write vcfanno config: %v", err)
	}
}

func prepareInputVCF(inputPath string, tempDir string) (string, error) {
	if strings.HasSuffix(inputPath, ".gz") {
		// Check for index
		if _, err := os.Stat(inputPath + ".tbi"); os.IsNotExist(err) {
			if _, err := os.Stat(inputPath + ".csi"); os.IsNotExist(err) {
				exec.Command("bcftools", "index", inputPath).Run()
			}
		}
		return inputPath, nil
	}

	// Compress and index
	outputPath := filepath.Join(tempDir, "input.vcf.gz")
	cmd := exec.Command("bgzip", "-c", inputPath)
	outFile, err := os.Create(outputPath)
	if err != nil {
		return "", err
	}
	cmd.Stdout = outFile
	if err := cmd.Run(); err != nil {
		outFile.Close()
		return "", err
	}
	outFile.Close()

	exec.Command("bcftools", "index", outputPath).Run()
	return outputPath, nil
}

func runVcfanno(tomlFile, inputVCF, outputVCF string) error {
	cmd := exec.Command("vcfanno", "-p", "4", tomlFile, inputVCF)
	outFile, err := os.Create(outputVCF)
	if err != nil {
		return err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

type filterStats struct {
	total     int
	annotated int
	high      int
	other     int
}

func filterByScores(annotatedVCF, finalVCF, nextVCF string, adaThresh, rfThresh float64) (filterStats, []string, error) {
	stats := filterStats{}
	var highScoringVariants []string

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

		// Write headers only to next VCF
		if strings.HasPrefix(line, "#") {
			nextWriter.WriteString(line + "\n")
			continue
		}

		stats.total++

		// Extract scores
		adaScore, rfScore := extractScores(line)
		
		if adaScore >= 0 || rfScore >= 0 {
			stats.annotated++
		}

		// Check if both scores exceed thresholds
		if adaScore > adaThresh && rfScore > rfThresh {
			finalWriter.WriteString(line + "\n")
			stats.high++

			// Extract variant info for reporting
			fields := strings.Split(line, "\t")
			if len(fields) >= 5 {
				highScoringVariants = append(highScoringVariants,
					fmt.Sprintf("%s:%s %s>%s [ada=%.3f, rf=%.3f]",
						fields[0], fields[1], fields[3], fields[4], adaScore, rfScore))
			}
		} else {
			nextWriter.WriteString(line + "\n")
			stats.other++
		}
	}

	return stats, highScoringVariants, scanner.Err()
}

func extractScores(line string) (float64, float64) {
	fields := strings.Split(line, "\t")
	if len(fields) < 8 {
		return -1, -1
	}

	info := fields[7]
	adaScore := -1.0
	rfScore := -1.0

	for _, part := range strings.Split(info, ";") {
		if strings.HasPrefix(part, "dbscSNV_ada_score=") {
			val := strings.TrimPrefix(part, "dbscSNV_ada_score=")
			if val != "." && val != "" {
				if f, err := strconv.ParseFloat(val, 64); err == nil {
					adaScore = f
				}
			}
		} else if strings.HasPrefix(part, "dbscSNV_rf_score=") {
			val := strings.TrimPrefix(part, "dbscSNV_rf_score=")
			if val != "." && val != "" {
				if f, err := strconv.ParseFloat(val, 64); err == nil {
					rfScore = f
				}
			}
		}
	}

	return adaScore, rfScore
}
