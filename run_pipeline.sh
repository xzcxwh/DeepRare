#!/bin/bash
#
# HERITA 变异注释与筛选流水线启动脚本
# 使用方法: ./run_pipeline.sh <input.vcf> [选项]
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_BIN="${SCRIPT_DIR}/herita-pipeline-go/bin/herita_pipeline"

# 默认参数
INPUT_VCF=""
OUTPUT_DIR="results"
THREADS=4
AF_THRESHOLD=0.01
SPLICEAI_THRESHOLD=0.5
TOP_N=50
SKIP_STEPS=""
DRY_RUN=""
VERBOSE=""

# 显示帮助
show_help() {
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║         HERITA 变异注释与筛选流水线                         ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    echo "用法: $0 <input.vcf> [选项]"
    echo ""
    echo "必需参数:"
    echo "  <input.vcf>        输入 VCF 文件"
    echo ""
    echo "可选参数:"
    echo "  -o, --output DIR   输出目录 (默认: results)"
    echo "  -t, --threads N    线程数 (默认: 4)"
    echo "  -a, --af FLOAT     gnomAD AF 阈值 (默认: 0.01)"
    echo "  -s, --spliceai F   SpliceAI 阈值 (默认: 0.5)"
    echo "  -n, --top N        TOP N 变异 (默认: 50)"
    echo "  --skip STEPS       跳过步骤 (逗号分隔)"
    echo "                     可选: gnomad,clinvar_annotate,clinvar_filter,"
    echo "                           snpeff,dbscsnv,vep,dbnsfp,scoring,merge"
    echo "  --dry-run          仅显示命令不执行"
    echo "  -v, --verbose      详细输出"
    echo "  -h, --help         显示帮助"
    echo ""
    echo "示例:"
    echo "  $0 input_data/HG001.vcf"
    echo "  $0 input_data/HG001.vcf -o my_results -t 8"
    echo "  $0 input_data/HG001.vcf --skip vep,dbnsfp"
    echo ""
}

# 解析参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -a|--af)
            AF_THRESHOLD="$2"
            shift 2
            ;;
        -s|--spliceai)
            SPLICEAI_THRESHOLD="$2"
            shift 2
            ;;
        -n|--top)
            TOP_N="$2"
            shift 2
            ;;
        --skip)
            SKIP_STEPS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="-dry-run"
            shift
            ;;
        -v|--verbose)
            VERBOSE="-verbose"
            shift
            ;;
        -*)
            echo "未知选项: $1"
            show_help
            exit 1
            ;;
        *)
            if [[ -z "$INPUT_VCF" ]]; then
                INPUT_VCF="$1"
            else
                echo "多余的参数: $1"
                exit 1
            fi
            shift
            ;;
    esac
done

# 检查输入文件
if [[ -z "$INPUT_VCF" ]]; then
    echo "错误: 必须指定输入 VCF 文件"
    show_help
    exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "错误: 输入文件不存在: $INPUT_VCF"
    exit 1
fi

# 检查流水线是否已编译
if [[ ! -x "$PIPELINE_BIN" ]]; then
    echo "流水线未编译，正在编译..."
    cd "${SCRIPT_DIR}/herita-pipeline-go"
    
    # 编译所有工具
    for tool in gnomad_grpmax_filter clinvar_annotate clinvar_filter snpeff_filter dbscsnv_filter vep_annotate dbnsfp_annotate finalists_merger herita_pipeline; do
        if [[ -d "cmd/$tool" ]]; then
            echo "  编译 $tool..."
            go build -o "bin/$tool" "./cmd/$tool/" 2>/dev/null || true
        fi
    done
    
    cd "$SCRIPT_DIR"
fi

# 构建命令
CMD="$PIPELINE_BIN"
CMD="$CMD -input $INPUT_VCF"
CMD="$CMD -output $OUTPUT_DIR"
CMD="$CMD -threads $THREADS"
CMD="$CMD -af $AF_THRESHOLD"
CMD="$CMD -spliceai $SPLICEAI_THRESHOLD"
CMD="$CMD -top $TOP_N"
CMD="$CMD -workdir $SCRIPT_DIR"

if [[ -n "$SKIP_STEPS" ]]; then
    CMD="$CMD -skip $SKIP_STEPS"
fi

if [[ -n "$DRY_RUN" ]]; then
    CMD="$CMD $DRY_RUN"
fi

if [[ -n "$VERBOSE" ]]; then
    CMD="$CMD $VERBOSE"
fi

# 运行流水线
echo ""
echo "开始运行 HERITA 流水线..."
echo "命令: $CMD"
echo ""

exec $CMD
