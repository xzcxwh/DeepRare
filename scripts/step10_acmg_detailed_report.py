#!/usr/bin/env python3
"""
Step 10: Generate Detailed ACMG Analysis Log
生成详细的 ACMG 分析日志，包含每个变异的证据详情和判定原因
"""
import os
import sys

# Ensure we can import step9
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

try:
    import step9_intervar_acmg as acmg
except ImportError:
    # Fallback if running from project root
    sys.path.append(os.path.join(os.getcwd(), 'scripts'))
    import step9_intervar_acmg as acmg

def generate_detailed_log(input_vcf, output_log):
    print(f"\n" + "="*70)
    print("Step 10: Generating Detailed ACMG Analysis Log")
    print("="*70)
    print(f"Input VCF: {input_vcf}")
    print(f"Output Log: {output_log}")
    
    # Load knowledge base
    print("\n[1/3] Loading knowledge base (reusing Step 9 logic)...")
    knowledge_base = {
        'lof_genes': acmg.load_lof_genes(),
        'ps1_changes': acmg.load_ps1_aa_changes(),
        'pm1_domains': acmg.load_pm1_domains(),
        'pp2_genes': acmg.load_pp2_genes(),
        'bp1_genes': acmg.load_bp1_genes(),
        'bs2_variants': acmg.load_bs2_variants(),
        'ps4_variants': acmg.load_ps4_variants(),
    }
    
    # Read variants
    print(f"\n[2/3] Reading variants from {input_vcf}...")
    variants = []
    if not os.path.exists(input_vcf):
        print(f"Error: Input file {input_vcf} not found.")
        return

    with open(input_vcf, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    variant = {
                        'chrom': parts[0],
                        'pos': parts[1],
                        'id': parts[2],
                        'ref': parts[3],
                        'alt': parts[4],
                        'qual': parts[5],
                        'filter': parts[6],
                        'info': acmg.parse_info_field(parts[7]),
                        'raw_info': parts[7]
                    }
                    variants.append(variant)
    
    print(f"  Processing {len(variants)} variants...")
    
    print("\n[3/3] Generating detailed log...")
    with open(output_log, 'w') as f:
        f.write("HERITA Pipeline - Step 10: Detailed ACMG Analysis Log\n")
        f.write(f"Date: {os.popen('date').read().strip()}\n")
        f.write(f"Input: {input_vcf}\n")
        f.write("=" * 80 + "\n\n")
        
        for i, var in enumerate(variants, 1):
            # Re-evaluate evidence to get details
            path_ev, ben_ev, details = acmg.evaluate_acmg_evidence(var, knowledge_base)
            classification, reason = acmg.classify_acmg(path_ev, ben_ev)
            
            gene = acmg.get_str(var['info'], 'genename', '.')
            
            f.write(f"Variant #{i}: {var['chrom']}:{var['pos']} {var['ref']}>{var['alt']}\n")
            f.write(f"Gene: {gene}\n")
            f.write("-" * 80 + "\n")
            f.write(f"ACMG Classification: {classification}\n")
            f.write("-" * 80 + "\n")
            
            f.write("1. Basic Information:\n")
            f.write(f"   - Position: {var['chrom']}:{var['pos']}\n")
            f.write(f"   - Ref/Alt: {var['ref']}/{var['alt']}\n")
            f.write(f"   - Gene: {gene}\n")
            aa_change = f"{acmg.get_str(var['info'], 'aaref', '')}{acmg.get_str(var['info'], 'aapos', '')}{acmg.get_str(var['info'], 'aaalt', '')}"
            if aa_change:
                f.write(f"   - AA Change: {aa_change}\n")
            
            f.write("\n2. Pathogenic Evidence:\n")
            if path_ev:
                for ev in sorted(list(set(path_ev))): # Unique and sorted
                    # Find details for this evidence code
                    ev_details = [d for d in details if d.startswith(ev)]
                    for d in ev_details:
                        f.write(f"   [+] {d}\n")
            else:
                f.write("   (None)\n")
            
            f.write("\n3. Benign Evidence:\n")
            if ben_ev:
                for ev in sorted(list(set(ben_ev))):
                    ev_details = [d for d in details if d.startswith(ev)]
                    for d in ev_details:
                        f.write(f"   [-] {d}\n")
            else:
                f.write("   (None)\n")
            
            f.write("\n4. Classification Reason:\n")
            f.write(f"   {reason}\n")
            
            f.write("\n5. Key Annotations:\n")
            keys = [
                ('gnomAD AF', 'gnomAD4.1_joint_POPMAX_AF'),
                ('ClinVar', 'clinvar_clnsig'),
                ('REVEL', 'REVEL_score'),
                ('CADD', 'CADD_phred'),
                ('MetaRNN', 'MetaRNN_score'),
                ('dbscSNV ADA', 'DBSCSNV_ADA'),
                ('dbscSNV RF', 'DBSCSNV_RF')
            ]
            for label, k in keys:
                val = var['info'].get(k, '.')
                if val != '.':
                    f.write(f"   - {label}: {val}\n")
            
            f.write("=" * 80 + "\n\n")
            
    print(f"\nDetailed log generated: {output_log}")
    print("✓ Step 10 Completed!")

if __name__ == "__main__":
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_vcf = os.path.join(project_dir, "results", "final_intervar_classified.vcf")
    output_log = os.path.join(project_dir, "results", "STEP10_ACMG_DETAILED_LOG.txt")
    
    # Allow command line overrides
    if len(sys.argv) > 1:
        input_vcf = sys.argv[1]
    if len(sys.argv) > 2:
        output_log = sys.argv[2]
        
    generate_detailed_log(input_vcf, output_log)
