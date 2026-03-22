#!/bin/bash

INPUT_DIR="/Volumes/T9/gnomAD/gnomAD_dataset"
OUTPUT_DIR="/Volumes/T9/gnomAD/gnomAD_lite"

# Set to 1 to avoid I/O contention on external drive
THREADS=1

echo "Starting generation of mini-gnomAD files (Sequential Mode)..."
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

    # Check if output exists and has a valid index (simple check)
    if [ -f "$output" ] && [ -f "$output.tbi" ]; then
        echo "✓ chr$chrom already exists and indexed, skipping."
        return
    fi

    echo "⏳ Processing chr$chrom..."
    
    # Use temp file to ensure atomicity
    temp_output="${output}.temp"
    
    # -x ^INFO/AF: Keep only AF
    if bcftools annotate -x ^INFO/AF --no-version -O z -o "$temp_output" "$input"; then
        mv "$temp_output" "$output"
        bcftools index -t "$output"
        echo "✓ Finished chr$chrom"
    else
        echo "❌ Failed to process chr$chrom"
        rm -f "$temp_output"
    fi
}

# Process chromosomes sequentially
for chrom in {1..22} X Y; do
    process_chrom "$chrom"
done

echo "All done!"
