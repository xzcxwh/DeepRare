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
	inputVCF  = flag.String("input", "", "Input VCF file")
	outputVCF = flag.String("output", "", "Output annotated VCF file")
	dbNSFP    = flag.String("dbnsfp", "", "dbNSFP TSV.GZ file path")
	threads   = flag.Int("threads", 4, "Number of threads for vcfanno")
)

func main() {
	flag.Parse()

	if *inputVCF == "" || *outputVCF == "" || *dbNSFP == "" {
		log.Fatal("Usage: dbnsfp_annotate -input <input.vcf> -output <output.vcf> -dbnsfp <dbNSFP.tsv.gz>")
	}

	startTime := time.Now()
	fmt.Println("============================================================")
	fmt.Println("HERITA dbNSFP Annotation Tool")
	fmt.Println("============================================================")
	fmt.Printf("Input: %s\n", *inputVCF)
	fmt.Printf("Output: %s\n", *outputVCF)
	fmt.Printf("dbNSFP: %s\n", *dbNSFP)
	fmt.Printf("Threads: %d\n", *threads)
	fmt.Println("============================================================")
	fmt.Println("Annotating fields:")
	fmt.Println("  - AlphaMissense_score")
	fmt.Println("  - REVEL_score")
	fmt.Println("  - MetaRNN_score")
	fmt.Println("  - PrimateAI_score")
	fmt.Println("  - CADD_phred")
	fmt.Println("  - ClinPred_score")
	fmt.Println("  - ESM1b_score")
	fmt.Println("============================================================")

	// Step 1: Create vcfanno config
	fmt.Println("\nStep 1: Creating vcfanno configuration...")
	configPath, err := createVcfannoConfig(*dbNSFP)
	if err != nil {
		log.Fatalf("Failed to create vcfanno config: %v", err)
	}
	defer os.Remove(configPath)
	fmt.Printf("Config created: %s\n", configPath)

	// Step 2: Prepare input VCF (bgzip and index if needed)
	fmt.Println("\nStep 2: Preparing input VCF...")
	preparedVCF, cleanup, err := prepareVCF(*inputVCF)
	if err != nil {
		log.Fatalf("Failed to prepare VCF: %v", err)
	}
	defer cleanup()
	fmt.Printf("Prepared VCF: %s\n", preparedVCF)

	// Step 3: Run vcfanno
	fmt.Println("\nStep 3: Running vcfanno annotation...")
	tempOutput := *outputVCF + ".temp"
	err = runVcfanno(configPath, preparedVCF, tempOutput, *threads)
	if err != nil {
		log.Fatalf("vcfanno failed: %v", err)
	}

	// Step 4: Post-process to add chr prefix back if needed and clean up
	fmt.Println("\nStep 4: Post-processing output...")
	err = postProcessVCF(tempOutput, *outputVCF, *inputVCF)
	if err != nil {
		log.Fatalf("Post-processing failed: %v", err)
	}
	os.Remove(tempOutput)

	// Step 5: Statistics
	fmt.Println("\nStep 5: Generating statistics...")
	printStats(*outputVCF)

	elapsed := time.Since(startTime)
	fmt.Println("\n============================================================")
	fmt.Printf("dbNSFP Annotation Complete! Time: %s\n", elapsed.Round(time.Second))
	fmt.Println("============================================================")
}

func createVcfannoConfig(dbNSFPPath string) (string, error) {
	// dbNSFP columns (1-based):
	// 1: #chr, 2: pos(1-based), 3: ref, 4: alt
	// 19: MetaRNN_score, 20: REVEL_score, 21: PrimateAI_score
	// 22: ClinPred_score, 23: ESM1b_score, 24: AlphaMissense_score, 25: CADD_phred

	config := fmt.Sprintf(`[[annotation]]
file="%s"
names=["dbNSFP_AlphaMissense", "dbNSFP_REVEL", "dbNSFP_MetaRNN", "dbNSFP_PrimateAI", "dbNSFP_CADD", "dbNSFP_ClinPred", "dbNSFP_ESM1b"]
columns=[24, 20, 19, 21, 25, 22, 23]
ops=["self", "self", "self", "self", "self", "self", "self"]
`, dbNSFPPath)

	configFile, err := os.CreateTemp("", "vcfanno_dbnsfp_*.toml")
	if err != nil {
		return "", err
	}
	defer configFile.Close()

	if _, err := configFile.WriteString(config); err != nil {
		return "", err
	}

	return configFile.Name(), nil
}

func prepareVCF(inputPath string) (string, func(), error) {
	cleanup := func() {}

	// Check if input has chr prefix
	hasChrPrefix, err := checkChrPrefix(inputPath)
	if err != nil {
		return "", cleanup, err
	}

	// dbNSFP uses chromosomes without chr prefix (1, 2, 3, ...)
	// So we need to strip chr prefix if present
	var vcfToProcess string

	if hasChrPrefix {
		// Strip chr prefix
		tempVCF := inputPath + ".nochr.vcf"
		err = stripChrPrefix(inputPath, tempVCF)
		if err != nil {
			return "", cleanup, err
		}
		vcfToProcess = tempVCF
		cleanup = func() { os.Remove(tempVCF) }
	} else {
		vcfToProcess = inputPath
	}

	// bgzip and index
	if !strings.HasSuffix(vcfToProcess, ".gz") {
		gzPath := vcfToProcess + ".gz"
		cmd := exec.Command("bgzip", "-c", vcfToProcess)
		outFile, err := os.Create(gzPath)
		if err != nil {
			return "", cleanup, err
		}
		cmd.Stdout = outFile
		if err := cmd.Run(); err != nil {
			outFile.Close()
			return "", cleanup, fmt.Errorf("bgzip failed: %v", err)
		}
		outFile.Close()

		// Index
		indexCmd := exec.Command("tabix", "-p", "vcf", gzPath)
		if err := indexCmd.Run(); err != nil {
			return "", cleanup, fmt.Errorf("tabix failed: %v", err)
		}

		oldCleanup := cleanup
		cleanup = func() {
			oldCleanup()
			os.Remove(gzPath)
			os.Remove(gzPath + ".tbi")
		}
		return gzPath, cleanup, nil
	}

	// Check if index exists
	if _, err := os.Stat(vcfToProcess + ".tbi"); os.IsNotExist(err) {
		indexCmd := exec.Command("tabix", "-p", "vcf", vcfToProcess)
		if err := indexCmd.Run(); err != nil {
			return "", cleanup, fmt.Errorf("tabix failed: %v", err)
		}
	}

	return vcfToProcess, cleanup, nil
}

func checkChrPrefix(vcfPath string) (bool, error) {
	file, err := os.Open(vcfPath)
	if err != nil {
		return false, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		// First data line
		return strings.HasPrefix(line, "chr"), nil
	}
	return false, nil
}

func stripChrPrefix(inputPath, outputPath string) error {
	inFile, err := os.Open(inputPath)
	if err != nil {
		return err
	}
	defer inFile.Close()

	outFile, err := os.Create(outputPath)
	if err != nil {
		return err
	}
	defer outFile.Close()

	scanner := bufio.NewScanner(inFile)
	writer := bufio.NewWriter(outFile)
	defer writer.Flush()

	// Increase buffer size for long lines
	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 10*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "##contig=<ID=chr") {
			// Update contig lines
			line = strings.Replace(line, "##contig=<ID=chr", "##contig=<ID=", 1)
		} else if !strings.HasPrefix(line, "#") && strings.HasPrefix(line, "chr") {
			// Strip chr from data lines
			line = line[3:]
		}
		fmt.Fprintln(writer, line)
	}

	return scanner.Err()
}

func runVcfanno(configPath, inputVCF, outputVCF string, threads int) error {
	cmd := exec.Command("vcfanno", "-p", fmt.Sprintf("%d", threads), configPath, inputVCF)
	
	outFile, err := os.Create(outputVCF)
	if err != nil {
		return err
	}
	defer outFile.Close()

	cmd.Stdout = outFile
	cmd.Stderr = os.Stderr

	fmt.Printf("Running: vcfanno -p %d %s %s\n", threads, filepath.Base(configPath), filepath.Base(inputVCF))
	
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("vcfanno execution failed: %v", err)
	}

	return nil
}

func postProcessVCF(tempVCF, outputVCF, originalVCF string) error {
	// Check if original had chr prefix
	hasChrPrefix, err := checkChrPrefix(originalVCF)
	if err != nil {
		return err
	}

	inFile, err := os.Open(tempVCF)
	if err != nil {
		return err
	}
	defer inFile.Close()

	outFile, err := os.Create(outputVCF)
	if err != nil {
		return err
	}
	defer outFile.Close()

	scanner := bufio.NewScanner(inFile)
	writer := bufio.NewWriter(outFile)
	defer writer.Flush()

	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 10*1024*1024)

	// Track if we've added INFO headers
	infoHeadersAdded := false

	for scanner.Scan() {
		line := scanner.Text()

		// Add INFO headers before #CHROM line
		if strings.HasPrefix(line, "#CHROM") && !infoHeadersAdded {
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_AlphaMissense,Number=.,Type=String,Description="AlphaMissense pathogenicity score from dbNSFP">`)
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_REVEL,Number=.,Type=String,Description="REVEL pathogenicity score from dbNSFP">`)
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_MetaRNN,Number=.,Type=String,Description="MetaRNN pathogenicity score from dbNSFP">`)
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_PrimateAI,Number=.,Type=String,Description="PrimateAI pathogenicity score from dbNSFP">`)
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_CADD,Number=.,Type=String,Description="CADD phred score from dbNSFP">`)
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_ClinPred,Number=.,Type=String,Description="ClinPred pathogenicity score from dbNSFP">`)
			fmt.Fprintln(writer, `##INFO=<ID=dbNSFP_ESM1b,Number=.,Type=String,Description="ESM1b pathogenicity score from dbNSFP">`)
			infoHeadersAdded = true
		}

		if hasChrPrefix && !strings.HasPrefix(line, "#") {
			// Add chr prefix back to data lines
			fields := strings.SplitN(line, "\t", 2)
			if len(fields) >= 2 && !strings.HasPrefix(fields[0], "chr") {
				line = "chr" + line
			}
		} else if hasChrPrefix && strings.HasPrefix(line, "##contig=<ID=") && !strings.HasPrefix(line, "##contig=<ID=chr") {
			// Add chr prefix to contig lines
			line = strings.Replace(line, "##contig=<ID=", "##contig=<ID=chr", 1)
		}

		fmt.Fprintln(writer, line)
	}

	return scanner.Err()
}

func printStats(vcfPath string) {
	file, err := os.Open(vcfPath)
	if err != nil {
		log.Printf("Warning: could not open output for stats: %v", err)
		return
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	buf := make([]byte, 0, 1024*1024)
	scanner.Buffer(buf, 10*1024*1024)

	var total, withAlphaMissense, withREVEL, withMetaRNN, withPrimateAI, withCADD, withClinPred, withESM1b int

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		total++

		if strings.Contains(line, "dbNSFP_AlphaMissense=") {
			withAlphaMissense++
		}
		if strings.Contains(line, "dbNSFP_REVEL=") {
			withREVEL++
		}
		if strings.Contains(line, "dbNSFP_MetaRNN=") {
			withMetaRNN++
		}
		if strings.Contains(line, "dbNSFP_PrimateAI=") {
			withPrimateAI++
		}
		if strings.Contains(line, "dbNSFP_CADD=") {
			withCADD++
		}
		if strings.Contains(line, "dbNSFP_ClinPred=") {
			withClinPred++
		}
		if strings.Contains(line, "dbNSFP_ESM1b=") {
			withESM1b++
		}
	}

	fmt.Println("\n============================================================")
	fmt.Println("dbNSFP Annotation Statistics")
	fmt.Println("============================================================")
	fmt.Printf("Total Variants: %d\n", total)
	fmt.Println("------------------------------------------------------------")
	fmt.Printf("With AlphaMissense: %d (%.1f%%)\n", withAlphaMissense, float64(withAlphaMissense)/float64(total)*100)
	fmt.Printf("With REVEL: %d (%.1f%%)\n", withREVEL, float64(withREVEL)/float64(total)*100)
	fmt.Printf("With MetaRNN: %d (%.1f%%)\n", withMetaRNN, float64(withMetaRNN)/float64(total)*100)
	fmt.Printf("With PrimateAI: %d (%.1f%%)\n", withPrimateAI, float64(withPrimateAI)/float64(total)*100)
	fmt.Printf("With CADD: %d (%.1f%%)\n", withCADD, float64(withCADD)/float64(total)*100)
	fmt.Printf("With ClinPred: %d (%.1f%%)\n", withClinPred, float64(withClinPred)/float64(total)*100)
	fmt.Printf("With ESM1b: %d (%.1f%%)\n", withESM1b, float64(withESM1b)/float64(total)*100)
}
