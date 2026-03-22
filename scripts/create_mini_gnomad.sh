#!/bin/bash

INPUT_DIR="/Volumes/T9/gnomAD/gnomAD_dataset"
OUTPUT_DIR="/Volumes/T9/gnomAD/gnomAD_lite"
THREADS=4

echo "Starting generation of mini-gnomAD files (keeping only INFO/AF)..."
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"

process_chrom() {
    chrom=$1
    input="$INPUT_DIR/gnomad.exomes.v4.1.sites.chr$chrom.vcf.bgz"
    output="$OUTPUT_DIR/gnomad.exomes.v4.1.sites.chr$chrom.mini.vcf.gz"
    
    if [ ! -f "$input" ]; then
        echo "⚠️  Input file for chr$chrom not found, skipping."
        return
    fi

    if [ -f "$output" ] && [ -f "$output.tbi" ]; then
        echo "✓ chr$chrom already exists, skipping."
        return
    fi

    echo "⏳ Processing chr$chrom..."
    # -x ^INFO/AF: Remove all INFO tags EXCEPT AF
    # --no-version: Don't append bcftools version line to header (cleaner)
    bcftools annotate -x ^INFO/AF --no-version -O z -o "$output" "$input"
    bcftools index -t "$output"
    echo "✓ Finished chr$chrom"
}

# Process chromosomes
# We use a simple semaphore approach for parallelism to avoid 'wait -n' compatibility issues
job_count=0
for chrom in {1..22} X Y; do
    process_chrom "$chrom" &
    ((job_count++))
    
    if [[ $job_count -ge $THREADS ]]; then
        wait
        job_count=0
    fi
done

wait
echo "All done!"
