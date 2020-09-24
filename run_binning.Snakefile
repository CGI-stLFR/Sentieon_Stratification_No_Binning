
import os.path

configfile: "run_binning.config"

stratification_map = {}
with open(config["stratification_region_tsv"]) as fh:
  for line in fh:
    name, strat_bed = line.rstrip().split('\t')
    stratification_map[name] = config["stratification_bed_dir"] + '/' + strat_bed

def get_all_output_files():
  '''
  Find all of the expected outputs
  '''
  expected_files = []
  for sample in config["samples"]:
    for ref_name in ("hs37d5", ):
      expected_files.append("synthetic_stratified_counts/" + sample + '/' + ref_name + "/all/counts.txt")
      expected_files.append("synthetic_truth_vcf_annotated/" + sample + '/' + ref_name + "/calls.vcf.gz")
      expected_files.append("stratified_counts_table/" + sample + '/' + ref_name + "/counts_table.txt")
      for hap in ['0', '1']:
        expected_files.append("synthetic_ref_align/" + sample + '/' + hap + '/' + ref_name + "/align.paf.gz")
        expected_files.append("synthetic_aligned/" + sample + '/' + hap + '/' + ref_name + "/aligned.bam")
      for truth_hap in ["hapA", "hapB"]:
        expected_files.append("synthetic_aligned/truth/" + truth_hap + '/' + ref_name +  "/aligned.bam")
  print(expected_files)
  return expected_files

rule all:
  input:
    all_outputs = get_all_output_files()


##################################
# Align each contig to reference #
##################################

rule ref_synthetic_align:
  input:
    input_fa = config["input_dir"] + "/{sample}/LinkedRead.seed_{hap}.fasta",
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    minimap_binary = config["minimap2"]
  output:
    aligned = "synthetic_ref_align/{sample}/{hap}/{ref_name}/align.paf.gz"
  threads:
    32
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.aligned})
    mkdir -p $outdir
    {input.minimap_binary} -x asm5 -t {threads} --paf-no-hit --secondary=no --cs -y {input.ref} <(sed 's/=/:i:/g' {input.input_fa}) | gzip -c > {output.aligned}
    """

##########################################
# Align everything to a reference genome #
##########################################


rule align_truths:
  input:
    input_fa = config['truth_dir'] + 'mergeScaftig_normalized_{truth_hap}.fa',
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    minimap_binary = config["minimap2"],
    k8_binary = config["k8"],
    sam_flt = config["sam_flt"],
    samtools = config["samtools"]
  output:
    aligned = "synthetic_aligned/{sample}/{truth_hap}/{ref_name}/aligned.bam"
  threads:
    32
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.aligned})
    mkdir -p $outdir
    {input.minimap_binary} -a -x asm5 -t 16 --paf-no-hit --secondary=no -r2k {input.ref} {input.input_fa} | {input.k8_binary} {input.sam_flt} /dev/stdin | {input.samtools} sort -m4G -@4 -o {output.aligned}
    """


rule align_samples:
  input:
    input_fa = config["input_dir"] + "/{sample}/LinkedRead.seed_{hap}.fasta",
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    minimap_binary = config["minimap2"],
    k8_binary = config["k8"],
    sam_flt = config["sam_flt"],
    samtools = config["samtools"]
  output:
    aligned = "synthetic_aligned/{sample}/{hap}/{ref_name}/aligned.bam"
  threads:
    32
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.aligned})
    mkdir -p $outdir
    {input.minimap_binary} -a -x asm5 -t 16 --paf-no-hit --secondary=no -r2k {input.ref} {input.input_fa} | {input.k8_binary} {input.sam_flt} /dev/stdin | {input.samtools} sort -m4G -@4 -o {output.aligned}
    """

##################################################################################
# Call variants from the aligned BAMs. Merge multiple samples in to a single VCF #
##################################################################################
rule call_all_against_ref:
  input:
    hap0_sample = "synthetic_aligned/{sample}/0/{ref_name}/aligned.bam",
    hap1_sample = "synthetic_aligned/{sample}/1/{ref_name}/aligned.bam",
    hapA_truth = "synthetic_aligned/truth/hapA/{ref_name}/aligned.bam",
    hapB_truth = "synthetic_aligned/truth/hapB/{ref_name}/aligned.bam",
    ref = lambda wildcards: os.path.realpath(config["ref_files"][wildcards.ref_name]),
    htsbox = config["htsbox"],
    sentieon = config["sentieon"]
  output:
    sample_truth_vcf = "synthetic_truth_vcf/{sample}/{ref_name}/calls.vcf.gz"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.sample_truth_vcf})
    mkdir -p $outdir
    {input.htsbox} pileup -q5 -evcf {input.ref} {input.hap0_sample} {input.hap1_sample} {input.hapA_truth} {input.hapB_truth} | {input.sentieon} util vcfconvert - {output.sample_truth_vcf}
    """

################################
# Annotate the called variants #
################################
rule annotate_against_ref:
  input:
    sample_truth_vcf = "synthetic_truth_vcf/{sample}/{ref_name}/calls.vcf.gz",
    label_calls_script = config["label_calls"]
  output:
    annotated_truth_vcf = "synthetic_truth_vcf_annotated/{sample}/{ref_name}/calls.vcf.gz"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.annotated_truth_vcf})
    mkdir -p $outdir
    export PYTHONPATH=/opt/sentieon-genomics-201911/lib/python/sentieon
    python3 {input.label_calls_script} {input.sample_truth_vcf} {output.annotated_truth_vcf}
    """

#################################################
# Get variant counts for each stratified region #
#################################################

rule get_filter_list:
  input:
    annotated_truth_vcf = "synthetic_truth_vcf_annotated/{sample}/{ref_name}/calls.vcf.gz",
    bcftools = config["bcftools"]
  output:
    filter_list = "synthetic_filters/{sample}/{ref_name}/filters.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.filter_list})
    mkdir -p $outdir
    {input.bcftools} view -h {input.annotated_truth_vcf} | grep "^##FILTER" | sed -e 's/^##FILTER=<ID=//' -e 's/,Description=.*$//' > {output.filter_list}
    """

rule filter_stratified:
  input:
    stratification_bed = lambda wildcards: stratification_map[wildcards.stratificaiton_name],
    annotated_truth_vcf = "synthetic_truth_vcf_annotated/{sample}/{ref_name}/calls.vcf.gz",
    bcftools = config["bcftools"],
    sentieon = config["sentieon"]
  output:
    stratified_truth_vcf = "synthetic_stratified_vcf/{sample}/{ref_name}/{stratificaiton_name}/stratified.vcf.gz"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.stratified_truth_vcf})
    mkdir -p $outdir
    {input.bcftools} view -T {input.stratification_bed} {input.annotated_truth_vcf} | {input.sentieon} util vcfconvert - $outdir/stratified.vcf.gz
    """

def get_count_vcf(wildcards):
  '''
  Find the stratified VCF from the wildcards
  '''
  if wildcards.stratificaiton_name == "all":
    return "synthetic_truth_vcf_annotated/" + wildcards.sample + '/' + wildcards.ref_name + "/calls.vcf.gz"
  else:
    return "synthetic_stratified_vcf/" + wildcards.sample + '/' + wildcards.ref_name + '/' + wildcards.stratificaiton_name + "/stratified.vcf.gz"

rule count_filters:
  input:
    filter_list = "synthetic_filters/{sample}/{ref_name}/filters.txt",
    vcf = get_count_vcf,
    bcftools = config["bcftools"]
  output:
    stratified_truth_counts = "synthetic_stratified_counts/{sample}/{ref_name}/{stratificaiton_name}/counts.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.stratified_truth_counts})
    mkdir -p $outdir
    echo "Total variants:" > {output.stratified_truth_counts}
    {input.bcftools} view -H {input.vcf} | wc -l >> {output.stratified_truth_counts}
    while IFS= read -r filter_name; do
      echo "Counting filter $filter_name" >> {output.stratified_truth_counts}
      {input.bcftools} view -H -f "$filter_name" {input.vcf} | wc -l >> {output.stratified_truth_counts}
    done < "{input.filter_list}"
    """

###########################################
# Turn the counts data in to a nice table #
###########################################
rule counts_to_table:
  input:
    "synthetic_stratified_counts/{sample}/{ref_name}/all/counts.txt",
    all_counts_files = expand("synthetic_stratified_counts/{sample}/{ref_name}/{stratificaiton_name}/counts.txt", sample="{sample}", ref_name="{ref_name}", stratificaiton_name=stratification_map.keys()),
    stratification_to_table = config["stratification_to_table"]
  output:
    stratified_counts_table = "stratified_counts_table/{sample}/{ref_name}/counts_table.txt"
  threads:
    1
  shell:
    """
    set -exvo pipefail

    outdir=$(dirname {output.stratified_counts_table})
    mkdir -p $outdir
    python3 {input.stratification_to_table} $(dirname $(dirname {input.all_counts_files[0]})) > {output.stratified_counts_table}
    """
