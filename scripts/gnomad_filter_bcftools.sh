#!/bin/bash
#
# 使用 bcftools 快速进行 gnomAD AF 注释和筛选
# 筛选条件: gnomAD_AF < 0.01 或 gnomAD_AF 为空
#
# 用法: ./gnomad_filter_bcftools.sh input.vcf output.vcf
#

set -e

INPUT_VCF="${1:-input_data/HG001_exonic.vcf}"
OUTPUT_VCF="${2:-input_data/HG001_exonic_rare.vcf}"
GNOMAD_DIR="gnomAD_dataset"
MAX_AF=0.01

echo "============================================================"
echo "gnomAD AF 注释与筛选 (bcftools)"
echo "============================================================"
echo "输入文件: $INPUT_VCF"
echo "输出文件: $OUTPUT_VCF"
echo "筛选条件: AF < $MAX_AF 或 AF 为空"
echo "============================================================"

# 统计输入
INPUT_COUNT=$(grep -cv "^#" "$INPUT_VCF" || echo 0)
echo "输入变异数: $INPUT_COUNT"

# 创建临时目录
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo ""
echo "[1/4] 排序、压缩并索引输入 VCF..."
# 先排序 VCF
bcftools sort "$INPUT_VCF" -O z -o "$TEMP_DIR/input.vcf.gz"
tabix -p vcf "$TEMP_DIR/input.vcf.gz"

echo "[2/4] 按染色体注释 gnomAD AF..."

# 创建输出头部
bcftools view -h "$TEMP_DIR/input.vcf.gz" | grep -v "^##INFO=<ID=gnomAD_AF" > "$TEMP_DIR/header.txt"
echo '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD v4.1 exomes allele frequency">' >> "$TEMP_DIR/header_info.txt"

# 获取所有染色体
CHROMS=$(bcftools view -H "$TEMP_DIR/input.vcf.gz" | cut -f1 | sort -u)

# 初始化输出
> "$TEMP_DIR/annotated.vcf"

for CHROM in $CHROMS; do
    CHROM_CLEAN="${CHROM#chr}"
    GNOMAD_FILE="$GNOMAD_DIR/gnomad.exomes.v4.1.sites.chr${CHROM_CLEAN}.vcf.bgz"
    
    if [[ ! -f "$GNOMAD_FILE" ]]; then
        echo "   ⚠️ $CHROM: 无 gnomAD 数据，保留所有变异"
        bcftools view -H "$TEMP_DIR/input.vcf.gz" -r "$CHROM" >> "$TEMP_DIR/annotated.vcf"
        continue
    fi
    
    # 提取该染色体的变异
    CHROM_COUNT=$(bcftools view -H "$TEMP_DIR/input.vcf.gz" -r "$CHROM" | wc -l | tr -d ' ')
    
    # 使用 bcftools annotate 添加 gnomAD AF
    bcftools annotate \
        -a "$GNOMAD_FILE" \
        -c INFO/gnomAD_AF:=INFO/AF \
        -r "$CHROM" \
        "$TEMP_DIR/input.vcf.gz" 2>/dev/null | \
        bcftools view -H >> "$TEMP_DIR/annotated.vcf"
    
    echo "   ✓ $CHROM: $CHROM_COUNT 变异已注释"
done

echo ""
echo "[3/4] 合并并添加头部..."
# 获取原始头部并添加 gnomAD_AF INFO
bcftools view -h "$TEMP_DIR/input.vcf.gz" | grep -v "^#CHROM" | grep -v "##INFO=<ID=gnomAD_AF" > "$TEMP_DIR/final_header.txt"
echo '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD v4.1 exomes allele frequency">' >> "$TEMP_DIR/final_header.txt"
bcftools view -h "$TEMP_DIR/input.vcf.gz" | grep "^#CHROM" >> "$TEMP_DIR/final_header.txt"

cat "$TEMP_DIR/final_header.txt" "$TEMP_DIR/annotated.vcf" > "$TEMP_DIR/annotated_full.vcf"

echo "[4/4] 筛选: AF < $MAX_AF 或 AF 为空..."
# 筛选条件: gnomAD_AF 不存在 或 gnomAD_AF < 0.01
bcftools filter \
    -i "INFO/gnomAD_AF<$MAX_AF || INFO/gnomAD_AF=\".\"" \
    -o "$OUTPUT_VCF" \
    "$TEMP_DIR/annotated_full.vcf"

# 统计输出
OUTPUT_COUNT=$(grep -cv "^#" "$OUTPUT_VCF" || echo 0)
FILTERED=$((INPUT_COUNT - OUTPUT_COUNT))

echo ""
echo "============================================================"
echo "结果统计:"
echo "   输入变异: $INPUT_COUNT"
echo "   输出变异: $OUTPUT_COUNT"
echo "   过滤掉 (AF >= $MAX_AF): $FILTERED"
echo "   保留比例: $(echo "scale=2; $OUTPUT_COUNT * 100 / $INPUT_COUNT" | bc)%"
echo "============================================================"
echo "✅ 完成! 输出: $OUTPUT_VCF"
