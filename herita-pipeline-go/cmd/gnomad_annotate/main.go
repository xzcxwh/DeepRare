package main

import (
"bufio"
"bytes"
"flag"
"fmt"
"io/ioutil"
"log"
"os"
"os/exec"
"path/filepath"
"runtime"
"strings"
"sync"
"time"
)

var (
inputVCF   = flag.String("input", "", "Input VCF file")
outputVCF  = flag.String("output", "", "Output VCF file")
	gnomadDir  = flag.String("gnomad", "/Volumes/T9/gnomAD/gnomAD_lite", "Directory containing gnomAD VCFs")
numWorkers = flag.Int("threads", runtime.NumCPU(), "Number of threads to use")
)

func main() {
	flag.Parse()

	if *inputVCF == "" || *outputVCF == "" {
		log.Fatal("Please provide -input and -output arguments")
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA gnomAD Annotation Tool (Go version)")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Output: %s\n", *outputVCF)
	fmt.Printf("GnomAD Dir: %s\n", *gnomadDir)
	fmt.Printf("Threads: %d\n", *numWorkers)
	fmt.Println("============================================================")

	// 1. Check dependencies
	checkDependency("bcftools")
	checkDependency("vcfanno")
	checkDependency("bgzip")

	// 2. Create temp directory
	tempDir, err := ioutil.TempDir("", "herita_gnomad_")
	if err != nil {
		log.Fatalf("Failed to create temp dir: %v", err)
	}
	defer os.RemoveAll(tempDir) // Clean up on exit

	// 3. Prepare Input VCF (Compress and Index if needed)
	workingVCF := *inputVCF
	if !strings.HasSuffix(*inputVCF, ".gz") {
		fmt.Println("Input VCF is not compressed. Creating temporary bgzipped version...")
		workingVCF = filepath.Join(tempDir, "input.vcf.gz")
		// Use bgzip -c to compress
		cmd := exec.Command("bgzip", "-c", *inputVCF)
		outfile, err := os.Create(workingVCF)
		if err != nil {
			log.Fatalf("Failed to create temp vcf: %v", err)
		}
		cmd.Stdout = outfile
		if err := cmd.Run(); err != nil {
			log.Fatalf("Failed to compress input VCF: %v", err)
		}
		outfile.Close()
	}
	
	ensureIndex(workingVCF)

	// 4. Get chromosomes present in input VCF
	chroms := getChromosomes(workingVCF)
	fmt.Printf("Found %d chromosomes in input VCF\n", len(chroms))

	// 5. Process chromosomes in parallel
	var wg sync.WaitGroup
	sem := make(chan struct{}, *numWorkers)
	results := make([]string, len(chroms))
	var resultsLock sync.Mutex

	// Map to store order
	chromOrder := make(map[string]int)
	for i, c := range chroms {
		chromOrder[c] = i
	}

	for _, chrom := range chroms {
		wg.Add(1)
		go func(c string) {
			defer wg.Done()
			sem <- struct{}{} // Acquire token
			defer func() { <-sem }() // Release token

			// Define file paths
			chromClean := strings.ReplaceAll(c, "chr", "")
			gnomadFile := filepath.Join(*gnomadDir, fmt.Sprintf("gnomad.exomes.v4.1.sites.chr%s.mini.vcf.gz", chromClean))
			
			// Check if gnomAD file exists
			if _, err := os.Stat(gnomadFile); os.IsNotExist(err) {
				fmt.Printf("⚠️  Warning: gnomAD file not found for %s, skipping annotation\n", c)
				// Just copy the input chunk to output chunk without annotation
				outFile := filepath.Join(tempDir, fmt.Sprintf("%s.vcf.gz", c))
				runCommand("bcftools", "view", "-r", c, "-O", "z", "-o", outFile, workingVCF)
				runCommand("bcftools", "index", outFile)
				
				resultsLock.Lock()
				results[chromOrder[c]] = outFile
				resultsLock.Unlock()
				return
			}

			// Extract chromosome from input
			inputChunk := filepath.Join(tempDir, fmt.Sprintf("input.%s.vcf.gz", c))
			runCommand("bcftools", "view", "-r", c, "-O", "z", "-o", inputChunk, workingVCF)
			runCommand("bcftools", "index", inputChunk)

			// Create vcfanno config
			configFile := filepath.Join(tempDir, fmt.Sprintf("config.%s.toml", c))
			configContent := fmt.Sprintf(`[[annotation]]
file="%s"
fields = ["AF"]
names = ["gnomAD_AF"]
ops = ["self"]
`, gnomadFile)
			ioutil.WriteFile(configFile, []byte(configContent), 0644)

			// Run vcfanno
			annotatedChunk := filepath.Join(tempDir, fmt.Sprintf("annotated.%s.vcf", c))
			// vcfanno writes to stdout
			cmd := exec.Command("vcfanno", "-p", "4", configFile, inputChunk)
			outfile, err := os.Create(annotatedChunk)
			if err != nil {
				log.Printf("Failed to create output file for %s: %v", c, err)
				return
			}
			cmd.Stdout = outfile
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				log.Printf("vcfanno failed for %s: %v", c, err)
				outfile.Close()
				return
			}
			outfile.Close()

			// Compress and index result
			finalChunk := filepath.Join(tempDir, fmt.Sprintf("%s.vcf.gz", c))
			runCommand("bcftools", "view", "-O", "z", "-o", finalChunk, annotatedChunk)
			runCommand("bcftools", "index", finalChunk)

			resultsLock.Lock()
			results[chromOrder[c]] = finalChunk
			resultsLock.Unlock()

			fmt.Printf("✓ %s processed\n", c)
		}(chrom)
	}

	wg.Wait()

	// 6. Concat results
	fmt.Println("Merging results...")
	var validResults []string
	for _, r := range results {
		if r != "" {
			validResults = append(validResults, r)
		}
	}

	if len(validResults) == 0 {
		log.Fatal("No results generated")
	}

	args := []string{"concat", "-a", "-O", "z", "-o", *outputVCF}
	args = append(args, validResults...)
	runCommand("bcftools", args...)
	runCommand("bcftools", "index", *outputVCF)

	// 7. Statistics
	fmt.Println("Calculating statistics...")
	total := countVariants(workingVCF)
	annotated := countAnnotated(*outputVCF)

	elapsed := time.Since(startTime)
	fmt.Println("============================================================")
	fmt.Println("Annotation Complete!")
	fmt.Printf("Total Time: %s\n", elapsed)
	fmt.Printf("Total Variants: %d\n", total)
	fmt.Printf("Annotated Variants (gnomAD_AF found): %d (%.2f%%)\n", annotated, float64(annotated)/float64(total)*100)
	fmt.Println("============================================================")
}

func checkDependency(tool string) {
	_, err := exec.LookPath(tool)
	if err != nil {
		log.Fatalf("Error: %s is not installed or not in PATH", tool)
	}
}

func ensureIndex(vcfPath string) {
	if _, err := os.Stat(vcfPath + ".tbi"); os.IsNotExist(err) {
		if _, err := os.Stat(vcfPath + ".csi"); os.IsNotExist(err) {
			fmt.Println("Indexing input VCF...")
			runCommand("bcftools", "index", vcfPath)
		}
	}
}

func getChromosomes(vcfPath string) []string {
	cmd := exec.Command("bcftools", "index", "-s", vcfPath)
	out, err := cmd.Output()
	if err != nil {
		log.Fatalf("Failed to list chromosomes: %v", err)
	}
	var chroms []string
	scanner := bufio.NewScanner(bytes.NewReader(out))
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) > 0 {
			chroms = append(chroms, fields[0])
		}
	}
	return chroms
}

func runCommand(name string, args ...string) {
	cmd := exec.Command(name, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		log.Fatalf("Command failed: %s %v\nError: %v", name, args, err)
	}
}

func countVariants(vcfPath string) int {
	// bcftools index -n is fast
	cmd := exec.Command("bcftools", "index", "-n", vcfPath)
	out, err := cmd.Output()
	if err == nil {
		var count int
		fmt.Sscanf(string(out), "%d", &count)
		return count
	}
	return 0
}

func countAnnotated(vcfPath string) int {
	// Count variants where gnomAD_AF is not "."
	cmd := exec.Command("bcftools", "query", "-f", "%INFO/gnomAD_AF\n", vcfPath)
	out, err := cmd.Output()
	if err != nil {
		return 0
	}
	count := 0
	scanner := bufio.NewScanner(bytes.NewReader(out))
	for scanner.Scan() {
		line := scanner.Text()
		if line != "." && line != "" {
			count++
		}
	}
	return count
}
