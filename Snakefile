
import os
import sys
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.1.2")

configfile: "config.yaml"

def get_samples(wildcards):
    return config["samples"][wildcards.sample]

def get_all_samples(wildcards):
    return config["samples"].values()

rule target:
	input:
		expand("results/sniffles_genotypes/{sample}.vcf",
			sample=config["samples"]),
		expand("results/mosdepth/{sample}.mosdepth.global.dist.txt",
			sample=config["samples"]),
		"results/mosdepth/regions.combined.gz"
#		"results/mosdepth_global_plot/global.html"

rule find_fastq:
	input:
		get_samples
	output:
		"fq/{sample}.fq"
	shell:
		"""
		zcat {input}  > {output}
		"""

rule ngmlr_map:
	input:
		fq = get_samples,
		REF = config["genome"]
	output:
		"results/align/{sample}.bam"
	log:
		"logs/align/{sample}.log"
	threads: 24
	conda: "envs/ngmlr.yaml"
	shell:
		"zcat {input.fq} | ngmlr --presets ont -t {threads} -r {input.REF} | \
		 samtools sort -@ 8 -o {output} - 2> {log}"

rule samtools_index:
	input:
		"results/align/{sample}.bam"
	output:
		"results/align/{sample}.bam.bai"
	log:
		"logs/align/samtools_index/{sample}.log"
	threads: 4
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools index -@ {threads} {input} 2> {log}"


rule sniffles_call:
	"""Call SVs per sample. Min support changed from default 10. """
	input:
		bam = "results/align/{sample}.bam",
		bai = "results/align/{sample}.bam.bai"
	output:
		"results/sniffles_calls/{sample}.vcf"
	params:
		read_support = 3
	log:
		"logs/sniffles_calls/{sample}.log"
	threads: 12
	conda:
		"envs/sniffles.yaml"
	shell:
		"sniffles --mapped_reads {input.bam} --vcf {output} --threads {threads} -s {params.read_support}  2> {log}"

rule survivor:
	input:
		[f"results/sniffles_calls/{sample}.vcf" for sample in config["samples"]]
	output:
		vcf = "results/sniffles_combined/calls.vcf",
		fofn = "results/sniffles_combined/samples.fofn"
	params:
		distance = config["parameters"]["survivor_distance"],
		caller_support = 1,
		same_type = 1,
		same_strand = -1,
		estimate_distance = -1,
		minimum_size = -1,
	conda: "envs/survivor.yaml"
	log:
		"logs/sniffles_combined/calls.vcf"
	shell:
		"ls {input} > {output.fofn} ; \
		SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
		{params.same_type} {params.same_strand} {params.estimate_distance}  \
		{params.minimum_size} {output.vcf} 2> {log}"

rule sniffles_genotype:
	input:
		bam = "results/align/{sample}.bam",
		ivcf = "results/sniffles_combined/calls.vcf"
	output:
		"results/sniffles_genotypes_temp/{sample}.vcf"
	conda: "envs/sniffles.yaml"
	log:
		"logs/sniffles_genotypes_temp/{sample}.log"
	shell:
		"sniffles --mapped_reads {input.bam} \
			--vcf {output} \
			--threads 22 \
			--cluster \
			--Ivcf {input.ivcf} 2> {log}"

rule bcftools_reheader_sniffles:
	"""Rule to be deleted as soon as ngmlr uses read groups correctly"""
	input:
		"results/sniffles_genotypes_temp/{sample}.vcf"
	output:
		vcf = "results/sniffles_genotypes/{sample}.vcf",
		sample = "results/sniffles_genotypes/sample_{sample}.txt"
	threads: 12
	log:
		"logs/bcftools_reheader/{sample}.log"
	conda: "envs/bcftools.yaml"
	shell:
		"""
		echo {wildcards.sample} > {output.sample} &&
		bcftools reheader -s {output.sample} {input} -o {output.vcf} 2> {log}
		"""

rule mosdepth_get:
	input:
		bam = "results/align/{sample}.bam", 
		bai = "results/align/{sample}.bam.bai"
	threads: 12
	output:
		"results/mosdepth/{sample}.mosdepth.global.dist.txt",
		"results/mosdepth/{sample}.regions.bed.gz"
	params:
		windowsize = 500,
		outdir = "mosdepth/{sample}"
	log:
		"logs/mosdepth/mosdepth_{sample}.log"
	conda: "envs/mosdepth.yaml"
	shell:
		"mosdepth --threads {threads} \
			-n \
			--by {params.windowsize} \
			{params.outdir} {input.bam} 2> {log}"
# not tested
rule mosdepth_combine:
	input:
		[f"results/mosdepth/{sample}.regions.bed.gz" for sample in config["samples"]]
	output:
		"results/mosdepth/regions.combined.gz"
	threads: 12
	log:
		"logs/mosdepth/mosdepth_combine.log"
	shell:
		os.path.join(workflow.basedir, "scripts/combine_mosdepth.py") + \
			" {input} -o {output} 2> {log}"

# not tested
rule mosdepth_global_plot:
	input:
		[f"results/mosdepth/{sample}.global.dist.txt" for sample in config["samples"]]
	output:
		"results/mosdepth_global_plot/global.html"
	threads: 12
	log:
		"logs/mosdepth/mosdepth_global_plot.log"
	shell:
		os.path.join(workflow.basedir, "scripts/mosdepth_plot-dist.py") + \
			" {input} -o {output} 2> {log}"
 
