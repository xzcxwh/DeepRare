package main

import (
	"bufio"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

var (
	inputVCF   = flag.String("input", "", "Input VCF file (e.g., HG001_cds.vcf)")
	outputVCF  = flag.String("output", "", "Output filtered VCF file")
	gnomadDir  = flag.String("gnomad", "/Volumes/T9/herita-project/gnomAD_dataset", "Directory containing gnomAD VCFs")
	maxAF      = flag.Float64("max-af", 0.01, "Maximum AF_grpmax threshold (variants with AF >= this will be filtered out)")
	numWorkers = flag.Int("threads", 0, "Number of threads (0 = auto)")
	keepTemp   = flag.Bool("keep-temp", false, "Keep temporary files for debugging")
)

func main() {
	flag.Parse()

	if *inputVCF == "" || *outputVCF == "" {
		log.Fatal("Usage: gnomad_grpmax_filter -input <input.vcf> -output <output.vcf> [-gnomad <dir>] [-max-af 0.01] [-threads N]")
	}

	if *numWorkers == 0 {
		*numWorkers = runtime.NumCPU()
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA gnomAD AF_grpmax Annotation & Filter Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Output: %s\n", *outputVCF)
	fmt.Printf("GnomAD Dir: %s\n", *gnomadDir)
	fmt.Printf("Max AF_grpmax: %g\n", *maxAF)
	fmt.Printf("Threads: %d\n", *numWorkers)
	fmt.Println("============================================================")

	// 1. Check dependencies
	checkDependency("bcftools")
	checkDependency("vcfanno")
	checkDependency("bgzip")
	checkDependency("tabix")

	// 2. Create temp directory
	tempDir, err := ioutil.TempDir("", "herita_grpmax_")
	if err != nil {
		log.Fatalf("Failed to create temp dir: %v", err)
	}
	if !*keepTemp {
		defer os.RemoveAll(tempDir)
	} else {
		fmt.Printf("Temp directory: %s\n", tempDir)
	}

	// 3. Prepare Input VCF (bgzip and index if needed)
	workingVCF, err := prepareInputVCF(*inputVCF, tempDir)
	if err != nil {
		log.Fatalf("Failed to prepare input VCF: %v", err)
	}

	// 4. Get chromosomes from input VCF
	chromosomes := getChromosomes(workingVCF)
	fmt.Printf("Found %d chromosomes in input VCF\n", len(chromosomes))

	// 5. Process each chromosome in parallel
	type result struct {
		chrom    string
		order    int
		vcfPath  string
		err      error
	}

	chromOrder := make(map[string]int)
	for i, c := range chromosomes {
		chromOrder[c] = i
	}

	resultsChan := make(chan result, len(chromosomes))
	semaphore := make(chan struct{}, *numWorkers)
	var wg sync.WaitGroup

	for _, chrom := range chromosomes {
		wg.Add(1)
		go func(c string) {
			defer wg.Done()
			semaphore <- struct{}{}
			defer func() { <-semaphore }()

			// Find gnomAD file for this chromosome
			gnomadFile := findGnomadFile(*gnomadDir, c)
			if gnomadFile == "" {
				fmt.Printf("⚠ No gnomAD file for %s, skipping annotation\n", c)
				// Extract chromosome without annotation
				chromVCF := filepath.Join(tempDir, fmt.Sprintf("%s_raw.vcf.gz", c))
				err := extractChromosome(workingVCF, c, chromVCF)
				if err != nil {
					resultsChan <- result{chrom: c, order: chromOrder[c], err: err}
					return
				}
				resultsChan <- result{chrom: c, order: chromOrder[c], vcfPath: chromVCF}
				return
			}

			// Create vcfanno config for this chromosome
			tomlFile := filepath.Join(tempDir, fmt.Sprintf("%s_config.toml", c))
			createVcfannoConfig(tomlFile, gnomadFile)

			// Extract chromosome from input
			chromInputVCF := filepath.Join(tempDir, fmt.Sprintf("%s_input.vcf.gz", c))
			err := extractChromosome(workingVCF, c, chromInputVCF)
			if err != nil {
				resultsChan <- result{chrom: c, order: chromOrder[c], err: err}
				return
			}

			// Run vcfanno
			annotatedVCF := filepath.Join(tempDir, fmt.Sprintf("%s_annotated.vcf", c))
			err = runVcfanno(tomlFile, chromInputVCF, annotatedVCF)
			if err != nil {
				resultsChan <- result{chrom: c, order: chromOrder[c], err: err}
				return
			}

			// Compress and index
			finalVCF := filepath.Join(tempDir, fmt.Sprintf("%s_final.vcf.gz", c))
			runCommandSilent("bgzip", "-c", annotatedVCF)
			cmd := exec.Command("bgzip", "-c", annotatedVCF)
			outFile, _ := os.Create(finalVCF)
			cmd.Stdout = outFile
			cmd.Run()
			outFile.Close()
			runCommandSilent("bcftools", "index", finalVCF)

			fmt.Printf("✓ %s processed\n", c)
			resultsChan <- result{chrom: c, order: chromOrder[c], vcfPath: finalVCF}
		}(chrom)
	}

	// Wait for all goroutines to finish
	go func() {
		wg.Wait()
		close(resultsChan)
	}()

	// Collect results
	results := make([]result, 0, len(chromosomes))
	for r := range resultsChan {
		if r.err != nil {
			log.Printf("Error processing %s: %v", r.chrom, r.err)
			continue
		}
		results = append(results, r)
	}

	// Sort by chromosome order
	sort.Slice(results, func(i, j int) bool {
		return results[i].order < results[j].order
	})

	// 6. Merge all chromosome results
	fmt.Println("Merging annotated results...")
	mergedVCF := filepath.Join(tempDir, "merged_annotated.vcf.gz")
	var vcfPaths []string
	for _, r := range results {
		if r.vcfPath != "" {
			vcfPaths = append(vcfPaths, r.vcfPath)
		}
	}

	if len(vcfPaths) == 0 {
		log.Fatal("No results to merge")
	}

	args := append([]string{"concat", "-a", "-O", "z", "-o", mergedVCF}, vcfPaths...)
	runCommand("bcftools", args...)
	runCommand("bcftools", "index", mergedVCF)

	// 7. Filter by AF_grpmax < 0.01 or AF_grpmax is empty
	fmt.Printf("Filtering variants with AF_grpmax < %g or missing...\n", *maxAF)
	filterByAFGrpmax(mergedVCF, *outputVCF, *maxAF)

	// 8. Statistics
	fmt.Println("Calculating statistics...")
	inputCount := countVariants(workingVCF)
	outputCount := countVariants(*outputVCF)
	annotatedCount := countAnnotatedVariants(*outputVCF)

	elapsed := time.Since(startTime)
	fmt.Println("============================================================")
	fmt.Println("Annotation & Filtering Complete!")
	fmt.Printf("Total Time: %s\n", elapsed.Round(time.Millisecond))
	fmt.Printf("Input Variants: %d\n", inputCount)
	fmt.Printf("Output Variants: %d (%.2f%% retained)\n", outputCount, float64(outputCount)/float64(inputCount)*100)
	fmt.Printf("Variants with AF_grpmax annotation: %d (%.2f%%)\n", annotatedCount, float64(annotatedCount)/float64(outputCount)*100)
	fmt.Println("============================================================")
}

func checkDependency(tool string) {
	_, err := exec.LookPath(tool)
	if err != nil {
		log.Fatalf("Error: %s is not installed or not in PATH", tool)
	}
}

func prepareInputVCF(inputPath string, tempDir string) (string, error) {
	if strings.HasSuffix(inputPath, ".gz") {
		// Already compressed, check for index
		if _, err := os.Stat(inputPath + ".tbi"); os.IsNotExist(err) {
			if _, err := os.Stat(inputPath + ".csi"); os.IsNotExist(err) {
				fmt.Println("Indexing input VCF...")
				runCommand("bcftools", "index", inputPath)
			}
		}
		return inputPath, nil
	}

	// Need to compress and index
	fmt.Println("Compressing and indexing input VCF...")
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

	runCommand("bcftools", "index", outputPath)
	return outputPath, nil
}

func getChromosomes(vcfPath string) []string {
	cmd := exec.Command("bcftools", "index", "-s", vcfPath)
	out, err := cmd.Output()
	if err != nil {
		log.Fatalf("Failed to list chromosomes: %v", err)
	}

	var chroms []string
	scanner := bufio.NewScanner(strings.NewReader(string(out)))
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) > 0 {
			chroms = append(chroms, fields[0])
		}
	}
	return chroms
}

func findGnomadFile(gnomadDir, chrom string) string {
	// Try different naming patterns
	patterns := []string{
		fmt.Sprintf("gnomad.exomes.v4.1.sites.%s.vcf.bgz", chrom),
		fmt.Sprintf("gnomad.exomes.v4.1.sites.%s.vcf.gz", chrom),
		fmt.Sprintf("gnomad.genomes.v4.1.sites.%s.vcf.bgz", chrom),
	}

	for _, pattern := range patterns {
		path := filepath.Join(gnomadDir, pattern)
		if _, err := os.Stat(path); err == nil {
			// Check for index
			if _, err := os.Stat(path + ".tbi"); err == nil {
				return path
			}
			if _, err := os.Stat(path + ".csi"); err == nil {
				return path
			}
		}
	}
	return ""
}

func createVcfannoConfig(tomlPath, gnomadFile string) {
	config := fmt.Sprintf(`[[annotation]]
file = "%s"
fields = ["AF_grpmax"]
names = ["gnomAD_AF_grpmax"]
ops = ["self"]
`, gnomadFile)

	if err := ioutil.WriteFile(tomlPath, []byte(config), 0644); err != nil {
		log.Fatalf("Failed to write vcfanno config: %v", err)
	}
}

func extractChromosome(inputVCF, chrom, outputVCF string) error {
	cmd := exec.Command("bcftools", "view", "-r", chrom, "-O", "z", "-o", outputVCF, inputVCF)
	if err := cmd.Run(); err != nil {
		return err
	}
	return exec.Command("bcftools", "index", outputVCF).Run()
}

func runVcfanno(tomlFile, inputVCF, outputVCF string) error {
	cmd := exec.Command("vcfanno", "-p", "2", tomlFile, inputVCF)
	outFile, err := os.Create(outputVCF)
	if err != nil {
		return err
	}
	defer outFile.Close()
	
	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr
	return cmd.Run()
}

func runCommand(name string, args ...string) {
	cmd := exec.Command(name, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		log.Fatalf("Command failed: %s %v\nError: %v", name, args, err)
	}
}

func runCommandSilent(name string, args ...string) error {
	cmd := exec.Command(name, args...)
	return cmd.Run()
}

func filterByAFGrpmax(inputVCF, outputVCF string, maxAF float64) {
	// Read input VCF
	cmd := exec.Command("bcftools", "view", inputVCF)
	stdout, err := cmd.StdoutPipe()
	if err != nil {
		log.Fatalf("Failed to create pipe: %v", err)
	}
	if err := cmd.Start(); err != nil {
		log.Fatalf("Failed to start bcftools: %v", err)
	}

	outFile, err := os.Create(outputVCF)
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	defer outFile.Close()

	writer := bufio.NewWriter(outFile)
	scanner := bufio.NewScanner(stdout)
	
	// Increase buffer size for long lines
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	kept := 0
	filtered := 0

	for scanner.Scan() {
		line := scanner.Text()
		
		if strings.HasPrefix(line, "#") {
			// Add INFO header for gnomAD_AF_grpmax if not present
			if strings.HasPrefix(line, "##INFO") && strings.Contains(line, "gnomAD_AF_grpmax") {
				// Already has it
			}
			if strings.HasPrefix(line, "#CHROM") {
				// Add header before #CHROM line
				writer.WriteString("##INFO=<ID=gnomAD_AF_grpmax,Number=A,Type=Float,Description=\"Maximum allele frequency across genetic ancestry groups from gnomAD v4.1\">\n")
			}
			writer.WriteString(line + "\n")
			continue
		}

		// Parse INFO field to get AF_grpmax
		fields := strings.Split(line, "\t")
		if len(fields) < 8 {
			continue
		}

		info := fields[7]
		afGrpmax := extractAFGrpmax(info)

		// Keep if AF_grpmax is missing (.) or < maxAF
		if afGrpmax < 0 || afGrpmax < maxAF {
			writer.WriteString(line + "\n")
			kept++
		} else {
			filtered++
		}
	}

	writer.Flush()
	cmd.Wait()

	fmt.Printf("  Kept: %d, Filtered: %d\n", kept, filtered)
}

func extractAFGrpmax(info string) float64 {
	for _, field := range strings.Split(info, ";") {
		if strings.HasPrefix(field, "gnomAD_AF_grpmax=") {
			val := strings.TrimPrefix(field, "gnomAD_AF_grpmax=")
			if val == "." || val == "" {
				return -1 // Missing value
			}
			// Handle multiple alleles (take max)
			maxVal := -1.0
			for _, v := range strings.Split(val, ",") {
				if v == "." || v == "" {
					continue
				}
				f, err := strconv.ParseFloat(v, 64)
				if err == nil && f > maxVal {
					maxVal = f
				}
			}
			return maxVal
		}
	}
	return -1 // Not found
}

func countVariants(vcfPath string) int {
	var cmd *exec.Cmd
	if strings.HasSuffix(vcfPath, ".gz") {
		cmd = exec.Command("bcftools", "view", "-H", vcfPath)
	} else {
		cmd = exec.Command("grep", "-v", "^#", vcfPath)
	}
	out, err := cmd.Output()
	if err != nil {
		return 0
	}
	return strings.Count(string(out), "\n")
}

func countAnnotatedVariants(vcfPath string) int {
	cmd := exec.Command("grep", "-c", "gnomAD_AF_grpmax=", vcfPath)
	out, _ := cmd.Output()
	count, _ := strconv.Atoi(strings.TrimSpace(string(out)))
	return count
}
