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
	inputVCF   = flag.String("input", "", "Input VCF file")
	outputVCF  = flag.String("output", "", "Output annotated VCF file")
	clinvarVCF = flag.String("clinvar", "", "ClinVar VCF file (bgzipped and indexed)")
	numWorkers = flag.Int("threads", 0, "Number of threads (0 = auto)")
	keepTemp   = flag.Bool("keep-temp", false, "Keep temporary files for debugging")
)

// ClinVar fields to annotate
var clinvarFields = []string{"CLNSIG", "CLNREVSTAT", "ALLELEID", "CLNDN", "CLNDISDB"}

func main() {
	flag.Parse()

	if *inputVCF == "" || *outputVCF == "" || *clinvarVCF == "" {
		log.Fatal("Usage: clinvar_annotate -input <input.vcf> -output <output.vcf> -clinvar <clinvar.vcf.gz> [-threads N]")
	}

	if *numWorkers == 0 {
		*numWorkers = runtime.NumCPU()
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA ClinVar Annotation Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Output: %s\n", *outputVCF)
	fmt.Printf("ClinVar: %s\n", *clinvarVCF)
	fmt.Printf("Threads: %d\n", *numWorkers)
	fmt.Printf("Fields: %s\n", strings.Join(clinvarFields, ", "))
	fmt.Println("============================================================")

	// 1. Check dependencies
	checkDependency("bcftools")
	checkDependency("vcfanno")
	checkDependency("bgzip")
	checkDependency("tabix")

	// 2. Verify ClinVar file exists and is indexed
	if _, err := os.Stat(*clinvarVCF); os.IsNotExist(err) {
		log.Fatalf("ClinVar file not found: %s", *clinvarVCF)
	}
	ensureIndex(*clinvarVCF)

	// 3. Create temp directory
	tempDir, err := ioutil.TempDir("", "herita_clinvar_")
	if err != nil {
		log.Fatalf("Failed to create temp dir: %v", err)
	}
	if !*keepTemp {
		defer os.RemoveAll(tempDir)
	} else {
		fmt.Printf("Temp directory: %s\n", tempDir)
	}

	// 4. Prepare input VCF (bgzip and index if needed)
	workingVCF, err := prepareInputVCF(*inputVCF, tempDir)
	if err != nil {
		log.Fatalf("Failed to prepare input VCF: %v", err)
	}

	// 5. Get chromosomes from input VCF
	chromosomes := getChromosomes(workingVCF)
	fmt.Printf("Found %d chromosomes in input VCF\n", len(chromosomes))

	// 6. Create vcfanno config
	tomlFile := filepath.Join(tempDir, "clinvar_config.toml")
	createVcfannoConfig(tomlFile, *clinvarVCF)

	// 7. Process each chromosome in parallel
	type result struct {
		chrom   string
		order   int
		vcfPath string
		err     error
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

			// Compress
			finalVCF := filepath.Join(tempDir, fmt.Sprintf("%s_final.vcf.gz", c))
			err = compressVCF(annotatedVCF, finalVCF)
			if err != nil {
				resultsChan <- result{chrom: c, order: chromOrder[c], err: err}
				return
			}

			fmt.Printf("✓ %s processed\n", c)
			resultsChan <- result{chrom: c, order: chromOrder[c], vcfPath: finalVCF}
		}(chrom)
	}

	// Wait for all goroutines
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

	// 8. Merge all chromosome results
	fmt.Println("Merging annotated results...")
	var vcfPaths []string
	for _, r := range results {
		if r.vcfPath != "" {
			vcfPaths = append(vcfPaths, r.vcfPath)
		}
	}

	if len(vcfPaths) == 0 {
		log.Fatal("No results to merge")
	}

	// Merge VCFs
	mergedVCF := filepath.Join(tempDir, "merged.vcf.gz")
	args := append([]string{"concat", "-a", "-O", "z", "-o", mergedVCF}, vcfPaths...)
	runCommand("bcftools", args...)

	// Decompress to final output (user requested plain VCF)
	if strings.HasSuffix(*outputVCF, ".gz") {
		// Output as compressed
		copyFile(mergedVCF, *outputVCF)
		runCommand("bcftools", "index", *outputVCF)
	} else {
		// Output as plain VCF
		runCommand("bcftools", "view", "-o", *outputVCF, mergedVCF)
	}

	// 9. Statistics
	fmt.Println("Calculating statistics...")
	inputCount := countVariants(workingVCF)
	outputCount := countVariants(*outputVCF)
	annotatedCount := countAnnotatedVariants(*outputVCF)

	elapsed := time.Since(startTime)
	fmt.Println("============================================================")
	fmt.Println("ClinVar Annotation Complete!")
	fmt.Printf("Total Time: %s\n", elapsed.Round(time.Millisecond))
	fmt.Printf("Input Variants: %d\n", inputCount)
	fmt.Printf("Output Variants: %d\n", outputCount)
	fmt.Printf("Variants with ClinVar annotation: %d (%.2f%%)\n", annotatedCount, float64(annotatedCount)/float64(outputCount)*100)
	fmt.Println("============================================================")
}

func checkDependency(tool string) {
	_, err := exec.LookPath(tool)
	if err != nil {
		log.Fatalf("Error: %s is not installed or not in PATH", tool)
	}
}

func ensureIndex(vcfPath string) {
	tbiPath := vcfPath + ".tbi"
	csiPath := vcfPath + ".csi"

	if _, err := os.Stat(tbiPath); err == nil {
		return
	}
	if _, err := os.Stat(csiPath); err == nil {
		return
	}

	fmt.Println("Indexing ClinVar VCF...")
	runCommand("bcftools", "index", "-t", vcfPath)
}

func prepareInputVCF(inputPath string, tempDir string) (string, error) {
	if strings.HasSuffix(inputPath, ".gz") {
		// Already compressed, check for index
		ensureIndex(inputPath)
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

func createVcfannoConfig(tomlPath, clinvarFile string) {
	// Build fields and names arrays
	var fields, names []string
	for _, f := range clinvarFields {
		fields = append(fields, fmt.Sprintf(`"%s"`, f))
		names = append(names, fmt.Sprintf(`"ClinVar_%s"`, f))
	}

	// All fields use "self" operation for exact matching
	ops := make([]string, len(clinvarFields))
	for i := range ops {
		ops[i] = `"self"`
	}

	config := fmt.Sprintf(`[[annotation]]
file = "%s"
fields = [%s]
names = [%s]
ops = [%s]
`, clinvarFile, strings.Join(fields, ", "), strings.Join(names, ", "), strings.Join(ops, ", "))

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

func compressVCF(inputVCF, outputVCF string) error {
	cmd := exec.Command("bgzip", "-c", inputVCF)
	outFile, err := os.Create(outputVCF)
	if err != nil {
		return err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	if err := cmd.Run(); err != nil {
		return err
	}

	return exec.Command("bcftools", "index", outputVCF).Run()
}

func copyFile(src, dst string) {
	input, err := ioutil.ReadFile(src)
	if err != nil {
		log.Fatalf("Failed to read file: %v", err)
	}
	if err := ioutil.WriteFile(dst, input, 0644); err != nil {
		log.Fatalf("Failed to write file: %v", err)
	}
}

func runCommand(name string, args ...string) {
	cmd := exec.Command(name, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		log.Fatalf("Command failed: %s %v\nError: %v", name, args, err)
	}
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
	// Count variants that have at least one ClinVar annotation
	cmd := exec.Command("grep", "-c", "ClinVar_CLNSIG=", vcfPath)
	out, _ := cmd.Output()
	count, _ := strconv.Atoi(strings.TrimSpace(string(out)))
	return count
}
