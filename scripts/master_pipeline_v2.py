#!/usr/bin/env python3
"""
HERITA 双分支流水线 v2
=====================
完整流程:
- Step0: CDS 交集预处理
- Step1: gnomAD 低频注释 (vcfanno)
- Step2: 频率过滤 (MAF < 0.01)
- Step3: dbNSFP 命中判定并拆分 (SNV+dbNSFP命中 -> Branch A; indel+未命中 -> Branch B)
- Branch A: 调用旧 run_pipeline.py 处理
- Branch B: snpEff -> 影响力过滤 -> VEP REST -> 致病性过滤 -> 整合注释
- Step10: 先整合 Branch B (snpEff+VEP), 再与 Branch A 合并 -> 最终输出
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import requests

# ============================================================
# 路径配置
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
CONFIG_DIR = os.path.join(PROJECT_ROOT, "config")
ANNOTATION_DIR = os.path.join(PROJECT_ROOT, "annotation_data")
INPUT_DIR = os.path.join(PROJECT_ROOT, "input_data")
TEMP_DIR = os.path.join(PROJECT_ROOT, "temp")

# snpEff 配置
SNPEFF_JAR = "/Volumes/T9/snpEff/snpEff/snpEff.jar"
JAVA_BIN = "/opt/homebrew/opt/openjdk@21/bin/java"

# dbNSFP 评分字段 (用于判定 dbNSFP 命中)
DBNSFP_SCORE_FIELDS = [
    "REVEL_score", "MetaRNN_score", "PrimateAI_score", "ClinPred_score",
    "AlphaMissense_score", "CADD_phred", "ESM1b_score"
]

# VEP REST API
VEP_API_URL = "https://rest.ensembl.org/vep/human/region"
VEP_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
VEP_PARAMS = {
    "AlphaMissense": 1,
    "CADD": 1,
    "REVEL": 1,
    "SpliceAI": 1,
    "dbNSFP": "phyloP100way_vertebrate",
    "hgvs": 1,
    "canonical": 1,
    "mane": 1,
    "numbers": 1,
    "variant_class": 1,
}


# 全局时间记录
STEP_TIMES = {}

def log(msg: str) -> None:
    timestamp = time.strftime('%H:%M:%S')
    print(f"[{timestamp}] {msg}", flush=True)

def log_step_start(step_name: str) -> None:
    """记录步骤开始时间"""
    STEP_TIMES[step_name] = {'start': time.time()}
    log(f"========== {step_name} 开始 ==========")

def log_step_end(step_name: str) -> None:
    """记录步骤结束时间"""
    if step_name in STEP_TIMES:
        STEP_TIMES[step_name]['end'] = time.time()
        elapsed = STEP_TIMES[step_name]['end'] - STEP_TIMES[step_name]['start']
        log(f"========== {step_name} 完成 (耗时: {elapsed:.2f}秒) ==========")
    else:
        log(f"========== {step_name} 完成 ==========")


def run_cmd(cmd: List[str], desc: str, stdout_path: Optional[str] = None, check: bool = True) -> int:
    """运行命令并记录"""
    log(f"[CMD] {desc}")
    log(f"      {' '.join(cmd[:6])}{'...' if len(cmd) > 6 else ''}")
    sys.stdout.flush()  # 确保日志立即输出
    
    stdout_handle = None
    try:
        if stdout_path:
            os.makedirs(os.path.dirname(stdout_path) or '.', exist_ok=True)
            stdout_handle = open(stdout_path, "w")
        # 让 stderr 直接显示到终端，便于调试
        result = subprocess.run(cmd, cwd=PROJECT_ROOT, stdout=stdout_handle, 
                                stderr=None, text=True, check=check)
        sys.stdout.flush()  # 命令完成后刷新输出
        return result.returncode
    except subprocess.CalledProcessError as e:
        log(f"      命令失败 (返回码 {e.returncode})")
        sys.stdout.flush()
        if check:
            raise RuntimeError(f"命令失败 (返回码 {e.returncode})")
        return e.returncode
    finally:
        if stdout_handle:
            stdout_handle.close()


def ensure_dirs() -> None:
    """确保必要目录存在"""
    for d in [RESULTS_DIR, TEMP_DIR]:
        os.makedirs(d, exist_ok=True)


def count_variants(vcf_path: str) -> int:
    """统计 VCF 中的变异数"""
    if not os.path.exists(vcf_path):
        return 0
    count = 0
    with open(vcf_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count


def parse_info(info: str) -> Dict[str, str]:
    """解析 INFO 字段"""
    out: Dict[str, str] = {}
    for entry in info.split(';'):
        if '=' in entry:
            k, v = entry.split('=', 1)
            out[k] = v
    return out


# ============================================================
# Step0: CDS 交集预处理
# ============================================================
def step0_cds_intersect(input_vcf: str, output_vcf: str) -> None:
    """使用 bedtools 与 CDS 区域取交集"""
    cds_bed = os.path.join(ANNOTATION_DIR, "hg38_cds.bed.gz")
    log(f"Step0: CDS 交集预处理")
    log(f"  输入: {input_vcf}")
    log(f"  CDS BED: {cds_bed}")
    
    run_cmd([
        "bedtools", "intersect", "-header", "-a", input_vcf, "-b", cds_bed
    ], "bedtools intersect", stdout_path=output_vcf)
    
    n = count_variants(output_vcf)
    log(f"  输出: {output_vcf} ({n} 个变异)")


# ============================================================
# Step1: gnomAD 精简版注释
# ============================================================
def step1_vcfanno_gnomad(input_vcf: str, output_vcf: str) -> None:
    """使用 Go 工具进行高效的 gnomAD 频率注释"""
    log(f"Step1: gnomAD 频率注释（精简版）")
    
    gnomad_tool = os.path.join(PROJECT_ROOT, "herita-pipeline-go/cmd/gnomad_annotate/gnomad_annotate")
    gnomad_dir = "/Volumes/T9/gnomAD/gnomAD_lite"
    
    # 使用临时输出文件（Go 工具输出 .vcf.gz）
    temp_output = output_vcf + ".gz"
    temp_log = os.path.join(TEMP_DIR, "gnomad_annotation.log")
    
    log(f"[CMD] gnomAD annotation (Go)")
    log(f"      {gnomad_tool} -input {input_vcf} -output {temp_output}...")
    log(f"      日志将写入: {temp_log}")
    sys.stdout.flush()
    
    cmd = [gnomad_tool, "-input", input_vcf, "-output", temp_output, 
           "-gnomad", gnomad_dir, "-threads", "8"]
    
    # 将输出写入临时文件,避免管道阻塞
    with open(temp_log, 'w') as log_file:
        result = subprocess.run(cmd, cwd=PROJECT_ROOT, stdout=log_file, 
                               stderr=subprocess.STDOUT, text=True)
    
    if result.returncode != 0:
        # 如果失败,打印日志
        with open(temp_log, 'r') as f:
            log("Go工具输出:")
            for line in f:
                print(line, end='')
        raise RuntimeError(f"Go gnomAD annotation failed with code {result.returncode}")
    
    # 打印最后几行日志确认完成
    log("Go工具执行完成,最后输出:")
    with open(temp_log, 'r') as f:
        lines = f.readlines()
        for line in lines[-10:]:
            print("  " + line, end='')
    
    # 解压缩为普通 VCF（保持与后续步骤兼容）
    log("解压缩输出文件...")
    with open(output_vcf, 'w') as fout:
        subprocess.run(['gunzip', '-c', temp_output], stdout=fout, check=True)
    
    # 清理临时文件
    if os.path.exists(temp_output):
        os.remove(temp_output)
    if os.path.exists(temp_output + ".csi"):
        os.remove(temp_output + ".csi")
    
    n = count_variants(output_vcf)
    log(f"  输出: {output_vcf} ({n} 个变异)")


# ============================================================
# Step2: 频率过滤
# ============================================================
def step2_frequency_filter(input_vcf: str, output_vcf: str, max_af: float = 0.01) -> None:
    """过滤高频变异 (gnomAD_AF >= max_af)，保留 AF < 0.01 或 AF 为空的变异"""
    log(f"Step2: 频率过滤 (gnomAD_AF < {max_af} 或为空)")
    
    keep = 0
    drop = 0
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            af_raw = info.get("gnomAD_AF", ".")
            try:
                af_val = float(af_raw) if af_raw not in ("", ".") else None
            except ValueError:
                af_val = None
            
            # 保留：AF < 0.01 或 AF 为空
            if af_val is None or af_val < max_af:
                fout.write(line)
                keep += 1
            else:
                drop += 1
    
    log(f"  保留: {keep} (AF<{max_af}或为空), 过滤: {drop}")
    log(f"  输出: {output_vcf}")


# ============================================================
# Step2.5: dbNSFP 注释（仅评分字段）
# ============================================================
def step2_5_annotate_dbnsfp(input_vcf: str, output_vcf: str) -> int:
    """
    使用 vcfanno 注释 dbNSFP 评分字段，用于后续判定命中情况
    """
    log(f"Step2.5: dbNSFP 评分注释")
    
    dbnsfp_file = os.path.join(INPUT_DIR, "dbNSFP5.3a_grch38_lite.tsv.gz")
    toml_config = os.path.join(CONFIG_DIR, "vcfanno_dbnsfp_scores.toml")
    
    # 创建临时 TOML 配置（仅注释评分字段）
    with open(toml_config, 'w') as f:
        f.write(f"""[[annotation]]
file = "{dbnsfp_file}"
columns = [1, 2, 3, 4, 19, 20, 21, 22, 23, 24, 25]
names = ["dbNSFP_chr", "dbNSFP_pos", "dbNSFP_ref", "dbNSFP_alt", "MetaRNN_score", "REVEL_score", "PrimateAI_score", "ClinPred_score", "ESM1b_score", "AlphaMissense_score", "CADD_phred"]
ops = ["self", "self", "self", "self", "self", "self", "self", "self", "self", "self", "self"]
""")
    
    cmd = [
        'vcfanno',
        '-lua', os.path.join(CONFIG_DIR, 'vcfanno_functions.lua'),
        toml_config,
        input_vcf
    ]
    
    run_cmd(cmd, "vcfanno dbNSFP scores", stdout_path=output_vcf)
    
    count = count_variants(output_vcf)
    log(f"  输出: {output_vcf} ({count} 个变异)")
    return count

# ============================================================
# Step3: 统计变异类型并检查 dbNSFP 命中，然后拆分
# ============================================================
def step3_split_by_dbnsfp(input_vcf: str, branch_a_vcf: str, branch_b_vcf: str) -> Tuple[int, int]:
    """
    统计并拆分变异:
    1. 统计 SNV 和 indel 数量
    2. 检查 dbNSFP 命中情况（检查 INFO 字段中的评分字段）
    3. 拆分:
       - Branch A: SNV 且被 dbNSFP 命中（有任一评分字段）
       - Branch B: indel + SNV 未被 dbNSFP 命中
    """
    log(f"Step3: 统计变异类型并检查 dbNSFP 命中")
    
    total_snv = 0
    total_indel = 0
    snv_dbnsfp_hit = 0
    snv_dbnsfp_miss = 0
    
    branch_a_lines = []
    branch_b_lines = []
    header_lines = []
    
    # dbNSFP 评分字段列表
    score_fields = ['MetaRNN_score', 'REVEL_score', 'PrimateAI_score', 
                    'ClinPred_score', 'ESM1b_score', 'AlphaMissense_score', 'CADD_phred']
    
    with open(input_vcf, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                header_lines.append(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            ref = fields[3]
            alts = fields[4].split(',')
            info = parse_info(fields[7])
            
            # 判断是否为 SNV
            is_snv = len(ref) == 1 and all(len(a) == 1 for a in alts)
            
            if is_snv:
                total_snv += 1
                
                # 检查是否有任一 dbNSFP 评分
                has_score = False
                for field in score_fields:
                    val = info.get(field, '.')
                    if val and val != '.' and val != '':
                        has_score = True
                        break
                
                if has_score:
                    snv_dbnsfp_hit += 1
                    branch_a_lines.append(line)
                else:
                    snv_dbnsfp_miss += 1
                    branch_b_lines.append(line)
            else:
                total_indel += 1
                branch_b_lines.append(line)
    
    log(f"  变异统计: SNV={total_snv}, indel={total_indel}, 总计={total_snv + total_indel}")
    log(f"  SNV dbNSFP 命中: {snv_dbnsfp_hit}/{total_snv}")
    log(f"  SNV dbNSFP 未命中: {snv_dbnsfp_miss}/{total_snv}")
    
    # 写入 Branch A（SNV + dbNSFP 命中）
    with open(branch_a_vcf, 'w') as fa:
        for h in header_lines:
            fa.write(h)
        for line in branch_a_lines:
            fa.write(line)
    
    # 写入 Branch B（indel + SNV 未命中）
    with open(branch_b_vcf, 'w') as fb:
        for h in header_lines:
            fb.write(h)
        for line in branch_b_lines:
            fb.write(line)
    
    branch_a_count = len(branch_a_lines)
    branch_b_count = len(branch_b_lines)
    
    log(f"  Branch A (SNV+dbNSFP命中): {branch_a_count} 个变异")
    log(f"  Branch B (indel+SNV未命中): {branch_b_count} 个变异")
    return branch_a_count, branch_b_count


# ============================================================
# Branch A: 调用旧流程
# ============================================================
def run_branch_a(branch_a_vcf: str, prefix: str, clean: bool = False) -> Optional[str]:
    """运行 Branch A 旧流程"""
    n = count_variants(branch_a_vcf)
    if n == 0:
        log("Branch A: 无变异，跳过")
        return None
    
    log(f"Branch A: 运行旧流程 ({n} 个变异)")
    
    # 备份原始 HG001.vcf
    original_input = os.path.join(INPUT_DIR, "HG001.vcf")
    backup_input = os.path.join(INPUT_DIR, "HG001.vcf.master_backup")
    if os.path.exists(original_input):
        shutil.copy2(original_input, backup_input)
    
    # 复制 Branch A 输入到 HG001.vcf
    shutil.copy2(branch_a_vcf, original_input)
    
    try:
        cmd = [sys.executable, os.path.join(SCRIPT_DIR, "run_pipeline.py")]
        if clean:
            cmd.append("--clean")
        run_cmd(cmd, "run_pipeline.py (Branch A)")
        
        # 复制结果
        legacy_out = os.path.join(RESULTS_DIR, "final_intervar_classified.vcf")
        if os.path.exists(legacy_out):
            dest = os.path.join(RESULTS_DIR, f"{prefix}_branchA_annotated.vcf")
            shutil.copy2(legacy_out, dest)
            log(f"  Branch A 输出: {dest}")
            return dest
        else:
            log("  警告: Branch A 输出文件不存在")
            return None
    finally:
        # 恢复原始输入文件
        if os.path.exists(backup_input):
            shutil.copy2(backup_input, original_input)


# ============================================================
# Branch B: snpEff + VEP 流程
# ============================================================
def run_branch_b_snpeff(input_vcf: str, output_vcf: str, genome: str = "GRCh38.99") -> None:
    """B1: snpEff 注释"""
    n = count_variants(input_vcf)
    log(f"Branch B - Step B1: snpEff 注释 ({n} 个变异)")
    
    if n == 0:
        shutil.copy2(input_vcf, output_vcf)
        return
    
    run_cmd([
        JAVA_BIN, "-Xmx4g", "-jar", SNPEFF_JAR,
        "-canon", "-noStats", "-hgvs", genome, input_vcf
    ], "snpEff annotation", stdout_path=output_vcf)
    
    log(f"  输出: {output_vcf}")


def run_branch_b_impact_filter(input_vcf: str, output_vcf: str, 
                                impacts: Tuple[str, ...] = ("HIGH", "MODERATE")) -> None:
    """B2: 影响力过滤 (保留 HIGH/MODERATE)"""
    log(f"Branch B - Step B2: 影响力过滤 (保留 {impacts})")
    
    keep = 0
    total = 0
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            total += 1
            fields = line.split('\t')
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            ann = info.get('ANN')
            if not ann:
                # 无 ANN 字段，保留 (可能是 indel 等)
                fout.write(line)
                keep += 1
                continue
            
            # 检查是否有 HIGH/MODERATE 影响
            accept = False
            for rec in ann.split(','):
                parts = rec.split('|')
                if len(parts) > 2 and parts[2] in impacts:
                    accept = True
                    break
            
            if accept:
                fout.write(line)
                keep += 1
    
    log(f"  保留: {keep}/{total}")
    log(f"  输出: {output_vcf}")


def run_branch_b_vep(input_vcf: str, output_vcf: str, batch_size: int = 100) -> None:
    """B3: VEP REST API 注释"""
    n = count_variants(input_vcf)
    log(f"Branch B - Step B3: VEP REST 注释 ({n} 个变异, batch={batch_size})")
    
    if n == 0:
        shutil.copy2(input_vcf, output_vcf)
        return
    
    # 解析变异
    variants = []
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0].replace('chr', '')
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(',')[0]
            
            # 构建 VEP 格式
            if len(ref) == 1 and len(alt) == 1:
                vep_str = f"{chrom} {pos} {pos} {ref}/{alt} 1"
            elif len(ref) > len(alt):
                deleted = ref[len(alt):]
                start = pos + len(alt)
                end = pos + len(ref) - 1
                vep_str = f"{chrom} {start} {end} {deleted}/- 1"
            else:
                inserted = alt[len(ref):]
                start = pos + len(ref)
                end = pos + len(ref) - 1
                vep_str = f"{chrom} {start} {end} -/{inserted} 1"
            
            variants.append({
                "key": f"{fields[0]}_{fields[1]}_{ref}_{alt}",
                "vep_string": vep_str,
                "line": line.strip(),
            })
    
    # 批量调用 VEP API
    vep_results: Dict[str, Dict[str, str]] = {}
    
    for i in range(0, len(variants), batch_size):
        batch = variants[i:i + batch_size]
        log(f"    VEP batch {i//batch_size + 1}/{(len(variants)-1)//batch_size + 1}")
        
        payload = {"variants": [v["vep_string"] for v in batch]}
        payload.update(VEP_PARAMS)
        
        for attempt in range(3):
            try:
                resp = requests.post(VEP_API_URL, headers=VEP_HEADERS, 
                                    json=payload, timeout=120)
                if resp.status_code == 200:
                    results = resp.json()
                    for r in results:
                        input_str = r.get('input')
                        if not input_str:
                            continue
                        for v in batch:
                            if v['vep_string'] == input_str:
                                vep_results[v['key']] = extract_vep_annotations(r)
                                break
                    break
                elif resp.status_code == 429:
                    wait = int(resp.headers.get('Retry-After', '5'))
                    log(f"    VEP 限流，等待 {wait}s...")
                    time.sleep(wait)
                else:
                    log(f"    VEP 错误: {resp.status_code}")
                    time.sleep(2)
            except Exception as e:
                log(f"    VEP 异常: {e}")
                time.sleep(2)
        
        time.sleep(1)  # 避免限流
    
    # 写入结果
    vep_info_headers = [
        '##INFO=<ID=VEP_ClinVar,Number=.,Type=String,Description="ClinVar significance from VEP">\n',
        '##INFO=<ID=VEP_AlphaMissense,Number=.,Type=Float,Description="AlphaMissense score">\n',
        '##INFO=<ID=VEP_SpliceAI,Number=.,Type=Float,Description="SpliceAI max delta">\n',
        '##INFO=<ID=VEP_CADD,Number=.,Type=Float,Description="CADD phred">\n',
        '##INFO=<ID=VEP_REVEL,Number=.,Type=Float,Description="REVEL score">\n',
        '##INFO=<ID=VEP_phyloP,Number=.,Type=Float,Description="phyloP100way">\n',
    ]
    
    with open(output_vcf, 'w') as fout:
        with open(input_vcf, 'r') as fin:
            header_written = False
            for line in fin:
                if line.startswith('##') and not header_written:
                    fout.write(line)
                elif line.startswith('#CHROM'):
                    for h in vep_info_headers:
                        fout.write(h)
                    fout.write(line)
                    header_written = True
                elif not line.startswith('#'):
                    fields = line.strip().split('\t')
                    key = f"{fields[0]}_{fields[1]}_{fields[3]}_{fields[4].split(',')[0]}"
                    ann = vep_results.get(key, {})
                    
                    # 添加 VEP 注释到 INFO
                    additions = []
                    for k, v in ann.items():
                        if v and v != '.':
                            additions.append(f"{k}={v}")
                    
                    if additions:
                        fields[7] = fields[7] + ';' + ';'.join(additions)
                    
                    fout.write('\t'.join(fields) + '\n')
    
    log(f"  输出: {output_vcf}")


def extract_vep_annotations(vep_result: Dict) -> Dict[str, str]:
    """从 VEP 结果提取注释"""
    ann = {
        "VEP_ClinVar": ".",
        "VEP_AlphaMissense": ".",
        "VEP_SpliceAI": ".",
        "VEP_CADD": ".",
        "VEP_REVEL": ".",
        "VEP_phyloP": ".",
    }
    
    if not vep_result:
        return ann
    
    # ClinVar from colocated_variants
    if 'colocated_variants' in vep_result:
        for cv in vep_result['colocated_variants']:
            if 'clin_sig' in cv and ann['VEP_ClinVar'] == '.':
                sig = cv['clin_sig']
                ann['VEP_ClinVar'] = ','.join(sig) if isinstance(sig, list) else str(sig)
    
    # 从 transcript_consequences 提取评分
    if 'transcript_consequences' in vep_result:
        for tc in vep_result['transcript_consequences']:
            # AlphaMissense
            if ann['VEP_AlphaMissense'] == '.':
                am = tc.get('alphamissense') or tc.get('am_pathogenicity')
                if am:
                    if isinstance(am, dict):
                        am = am.get('am_pathogenicity')
                    if am:
                        ann['VEP_AlphaMissense'] = str(am)
            
            # CADD
            if ann['VEP_CADD'] == '.' and 'cadd_phred' in tc:
                val = tc['cadd_phred']
                if val is not None:
                    ann['VEP_CADD'] = str(val)
            
            # REVEL
            if ann['VEP_REVEL'] == '.' and 'revel' in tc:
                val = tc['revel']
                if val is not None and str(val) != '.':
                    ann['VEP_REVEL'] = str(val)
            
            # SpliceAI
            if ann['VEP_SpliceAI'] == '.' and 'spliceai' in tc:
                spliceai = tc['spliceai']
                if isinstance(spliceai, dict):
                    scores = []
                    for key in ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']:
                        if key in spliceai and spliceai[key] is not None:
                            try:
                                scores.append(float(spliceai[key]))
                            except (ValueError, TypeError):
                                pass
                    if scores:
                        ann['VEP_SpliceAI'] = f"{max(scores):.4f}"
            
            # phyloP
            if ann['VEP_phyloP'] == '.' and 'phylop100way_vertebrate' in tc:
                val = tc['phylop100way_vertebrate']
                if val is not None:
                    ann['VEP_phyloP'] = str(val)
    
    return ann


def run_branch_b_pathogenic_filter(input_vcf: str, output_vcf: str,
                                    thresholds: Dict[str, float] = None,
                                    min_hits: int = 2) -> None:
    """B4: 致病性过滤 (满足 >= min_hits 个条件)"""
    if thresholds is None:
        thresholds = {
            "VEP_CADD": 20.0,
            "VEP_SpliceAI": 0.5,
            "VEP_REVEL": 0.5,
            "VEP_AlphaMissense": 0.5,
        }
    
    log(f"Branch B - Step B4: 致病性过滤 (>= {min_hits} 个条件)")
    
    keep = 0
    total = 0
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            total += 1
            fields = line.split('\t')
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            
            # 计算命中数
            hits = 0
            for k, threshold in thresholds.items():
                v = info.get(k)
                if v and v != '.':
                    try:
                        if float(v) >= threshold:
                            hits += 1
                    except ValueError:
                        pass
            
            # 如果有 ClinVar 致病/可能致病也算命中
            clinvar = info.get('VEP_ClinVar', '')
            if 'pathogenic' in clinvar.lower():
                hits += 1
            
            if hits >= min_hits:
                fout.write(line)
                keep += 1
    
    log(f"  保留: {keep}/{total}")
    log(f"  输出: {output_vcf}")


def run_branch_b_merge_annotations(snpeff_filtered_vcf: str, vep_pathogenic_vcf: str, output_vcf: str) -> None:
    """B5: 整合 snpEff HIGH 影响力变异 和 VEP 致病性变异"""
    log(f"Branch B - Step B5: 整合 snpEff HIGH + VEP 致病性变异")
    
    # 1. 收集所有 header 行
    headers = []
    with open(snpeff_filtered_vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                if line not in headers:
                    headers.append(line)
            elif line.startswith('#CHROM'):
                header_line = line
                break
    
    # 添加 VEP headers
    with open(vep_pathogenic_vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                if 'VEP' in line and line not in headers:
                    headers.append(line)
            else:
                break
    
    # 2. 提取 snpEff HIGH 影响力变异
    snpeff_high_variants = {}
    with open(snpeff_filtered_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            ann = info.get('ANN', '')
            # 检查是否包含 HIGH 影响力
            if '|HIGH|' in ann:
                key = f"{fields[0]}_{fields[1]}_{fields[3]}_{fields[4]}"
                snpeff_high_variants[key] = line
    
    # 3. 收集 VEP 致病性变异（包含 VEP 注释）
    vep_variants = {}
    with open(vep_pathogenic_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            key = f"{fields[0]}_{fields[1]}_{fields[3]}_{fields[4]}"
            vep_variants[key] = line
    
    # 4. 合并：先写入 VEP 变异（有完整注释），再添加 snpEff HIGH 独有的变异
    log(f"  snpEff HIGH: {len(snpeff_high_variants)} 个变异")
    log(f"  VEP 致病性: {len(vep_variants)} 个变异")
    
    all_variants = {}
    
    # VEP 变异优先（包含更完整的注释）
    for key, line in vep_variants.items():
        all_variants[key] = line
    
    # 添加 snpEff HIGH 独有的变异
    for key, line in snpeff_high_variants.items():
        if key not in all_variants:
            all_variants[key] = line
    
    # 5. 写入输出文件
    with open(output_vcf, 'w') as fout:
        # 写入 headers
        for h in sorted(set(headers)):
            fout.write(h)
        fout.write(header_line)
        
        # 写入变异（按染色体位置排序）
        sorted_keys = sorted(all_variants.keys(), key=lambda x: (x.split('_')[0], int(x.split('_')[1])))
        for key in sorted_keys:
            fout.write(all_variants[key])
    
    n = count_variants(output_vcf)
    log(f"  合并后总计: {n} 个变异")
    log(f"  输出: {output_vcf}")


# ============================================================
# Step10: 最终整合
# ============================================================
def sort_vcf_by_chromosome(input_vcf: str, output_vcf: str) -> None:
    """按染色体和位置排序 VCF 文件"""
    log("  对 VCF 文件按染色体排序...")
    
    # 定义染色体排序顺序
    chrom_order = {f"chr{i}": i for i in range(1, 23)}
    chrom_order['chrX'] = 23
    chrom_order['chrY'] = 24
    chrom_order['chrM'] = 25
    
    def chrom_sort_key(chrom: str) -> int:
        return chrom_order.get(chrom, 999)
    
    # 读取所有内容
    headers = []
    variants = []
    
    with open(input_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line)
            else:
                variants.append(line)
    
    # 按染色体和位置排序
    variants.sort(key=lambda x: (chrom_sort_key(x.split('\t')[0]), int(x.split('\t')[1])))
    
    # 写入排序后的文件
    with open(output_vcf, 'w') as f:
        for h in headers:
            f.write(h)
        for v in variants:
            f.write(v)
    
    log(f"  排序完成: {len(variants)} 个变异")

def step10_final_merge(branch_a_vcf: Optional[str], branch_b_vcf: Optional[str], 
                       output_vcf: str) -> None:
    """合并 Branch A 和 Branch B 结果"""
    log(f"Step10: 最终整合")
    
    # 收集所有变异
    variants: Dict[Tuple[str, str, str, str], str] = {}
    
    # 先读取 Branch B (优先)
    if branch_b_vcf and os.path.exists(branch_b_vcf):
        n_b = 0
        with open(branch_b_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    key = (fields[0], fields[1], fields[3], fields[4])
                    variants[key] = line.strip()
                    n_b += 1
        log(f"  Branch B: {n_b} 个变异")
    
    # 再读取 Branch A
    if branch_a_vcf and os.path.exists(branch_a_vcf):
        n_a = 0
        with open(branch_a_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    key = (fields[0], fields[1], fields[3], fields[4])
                    if key not in variants:  # Branch B 优先
                        variants[key] = line.strip()
                        n_a += 1
        log(f"  Branch A: {n_a} 个新变异")
    
    # 写入合并结果
    # 使用 Branch A 或 Branch B 的 header
    header_source = branch_a_vcf if branch_a_vcf and os.path.exists(branch_a_vcf) else branch_b_vcf
    
    with open(output_vcf, 'w') as fout:
        if header_source and os.path.exists(header_source):
            with open(header_source, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        fout.write(line)
                    else:
                        break
        
        # 先写入未排序的临时文件
        temp_vcf = output_vcf + '.unsorted'
        with open(temp_vcf, 'w') as ftmp:
            if header_source and os.path.exists(header_source):
                with open(header_source, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            ftmp.write(line)
                        else:
                            break
            
            for key in variants.keys():
                ftmp.write(variants[key] + '\n')
        
        # 排序后写入最终文件
        sort_vcf_by_chromosome(temp_vcf, output_vcf)
        os.remove(temp_vcf)
    
    n = count_variants(output_vcf)
    log(f"  最终输出: {output_vcf} ({n} 个变异)")


# ============================================================
# 验证结果
# ============================================================
def validate_results(output_vcf: str) -> None:
    """验证最终结果"""
    log(f"验证结果: {output_vcf}")
    
    if not os.path.exists(output_vcf):
        log("  错误: 输出文件不存在")
        return
    
    total = 0
    ann_count = 0
    vep_cadd = 0
    vep_spliceai = 0
    vep_clinvar = 0
    
    with open(output_vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            total += 1
            if 'ANN=' in line:
                ann_count += 1
            if 'VEP_CADD=' in line:
                vep_cadd += 1
            if 'VEP_SpliceAI=' in line:
                vep_spliceai += 1
            if 'VEP_ClinVar=' in line:
                vep_clinvar += 1
    
    log(f"  总变异数: {total}")
    log(f"  snpEff ANN: {ann_count}/{total}")
    log(f"  VEP_CADD: {vep_cadd}/{total}")
    log(f"  VEP_SpliceAI: {vep_spliceai}/{total}")
    log(f"  VEP_ClinVar: {vep_clinvar}/{total}")


# ============================================================
# 主函数
# ============================================================
def main() -> None:
    start_time = time.time()
    
    parser = argparse.ArgumentParser(description="HERITA 双分支流水线 v2")
    parser.add_argument("input_vcf", nargs="?", 
                        default=os.path.join(INPUT_DIR, "HG001.vcf"),
                        help="输入 VCF 文件")
    parser.add_argument("--prefix", default=None, help="输出文件前缀")
    parser.add_argument("--max-af", type=float, default=0.01, help="MAF 阈值")
    parser.add_argument("--snpeff-genome", default="GRCh38.99", help="snpEff 基因组")
    parser.add_argument("--vep-batch", type=int, default=100, help="VEP 批次大小")
    parser.add_argument("--pathogenic-min-hits", type=int, default=1, help="致病性最小命中数")
    parser.add_argument("--clean-branch-a", action="store_true", help="Branch A 清理旧结果")
    parser.add_argument("--skip-branch-a", action="store_true", help="跳过 Branch A")
    parser.add_argument("--skip-branch-b", action="store_true", help="跳过 Branch B")
    args = parser.parse_args()
    
    input_vcf = os.path.abspath(args.input_vcf)
    prefix = args.prefix or os.path.splitext(os.path.basename(input_vcf))[0]
    
    ensure_dirs()
    
    log("=" * 60)
    log("HERITA 双分支流水线 v2")
    log(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"输入: {input_vcf}")
    log(f"前缀: {prefix}")
    log(f"参数: MAF<{args.max_af}, VEP批次={args.vep_batch}, 致病性最小命中={args.pathogenic_min_hits}")
    log("=" * 60)
    
    # 定义中间文件路径
    step0_out = os.path.join(RESULTS_DIR, f"{prefix}_cds.vcf")
    step1_out = os.path.join(RESULTS_DIR, f"{prefix}_cds_annotated.vcf")
    step2_out = os.path.join(RESULTS_DIR, f"{prefix}_cds_filtered.vcf")
    branch_a_in = os.path.join(RESULTS_DIR, f"{prefix}_branchA_input.vcf")
    branch_b_in = os.path.join(RESULTS_DIR, f"{prefix}_branchB_input.vcf")
    
    # Branch B 中间文件
    b1_snpeff = os.path.join(RESULTS_DIR, f"{prefix}_branchB_snpeff.vcf")
    b2_filtered = os.path.join(RESULTS_DIR, f"{prefix}_branchB_snpeff_filtered.vcf")
    b3_vep = os.path.join(RESULTS_DIR, f"{prefix}_branchB_vep.vcf")
    b4_pathogenic = os.path.join(RESULTS_DIR, f"{prefix}_branchB_pathogenic.vcf")
    b5_annotated = os.path.join(RESULTS_DIR, f"{prefix}_branchB_annotated.vcf")
    
    # Branch A 输出
    branch_a_out = os.path.join(RESULTS_DIR, f"{prefix}_branchA_annotated.vcf")
    
    # 最终输出
    final_out = os.path.join(RESULTS_DIR, f"{prefix}_final_annotated.vcf")
    
    # 执行流程
    log_step_start("Step0: CDS交集预处理")
    step0_cds_intersect(input_vcf, step0_out)
    log_step_end("Step0: CDS交集预处理")
    
    log_step_start("Step1: gnomAD频率注释")
    step1_vcfanno_gnomad(step0_out, step1_out)
    log_step_end("Step1: gnomAD频率注释")
    
    log_step_start("Step2: 频率过滤")
    step2_frequency_filter(step1_out, step2_out, args.max_af)
    log_step_end("Step2: 频率过滤")
    
    log_step_start("Step2.5: dbNSFP评分注释")
    step2_5_out = os.path.join(RESULTS_DIR, f"{prefix}_cds_filtered_dbnsfp.vcf")
    step2_5_annotate_dbnsfp(step2_out, step2_5_out)
    log_step_end("Step2.5: dbNSFP评分注释")
    
    log_step_start("Step3: 变异分类拆分")
    step3_split_by_dbnsfp(step2_5_out, branch_a_in, branch_b_in)
    log_step_end("Step3: 变异分类拆分")
    
    # Branch A
    branch_a_result = None
    if not args.skip_branch_a:
        log_step_start("Branch A: 旧流程处理")
        branch_a_result = run_branch_a(branch_a_in, prefix, args.clean_branch_a)
        log_step_end("Branch A: 旧流程处理")
    
    # Branch B
    branch_b_result = None
    if not args.skip_branch_b:
        n_b = count_variants(branch_b_in)
        if n_b > 0:
            log_step_start("Branch B: snpEff+VEP处理")
            run_branch_b_snpeff(branch_b_in, b1_snpeff, args.snpeff_genome)
            run_branch_b_impact_filter(b1_snpeff, b2_filtered)
            run_branch_b_vep(b2_filtered, b3_vep, args.vep_batch)
            run_branch_b_pathogenic_filter(b3_vep, b4_pathogenic, min_hits=args.pathogenic_min_hits)
            run_branch_b_merge_annotations(b2_filtered, b4_pathogenic, b5_annotated)
            branch_b_result = b5_annotated
            log_step_end("Branch B: snpEff+VEP处理")
        else:
            log("Branch B: 无变异，跳过")
    
    # Step10: 最终整合
    log_step_start("Step10: 最终整合与排序")
    step10_final_merge(branch_a_result, branch_b_result, final_out)
    log_step_end("Step10: 最终整合与排序")
    
    # 验证
    validate_results(final_out)
    
    # 生成运行报告
    total_time = time.time() - start_time
    log("=" * 60)
    log("流水线执行完成!")
    log(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"总耗时: {total_time:.2f}秒 ({total_time/60:.2f}分钟)")
    log("")
    log("各步骤耗时详情:")
    for step_name, times in STEP_TIMES.items():
        if 'end' in times:
            elapsed = times['end'] - times['start']
            log(f"  {step_name}: {elapsed:.2f}秒")
    log("")
    log(f"最终输出: {final_out}")
    log("=" * 60)


if __name__ == "__main__":
    main()
