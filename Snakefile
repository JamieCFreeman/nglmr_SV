
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
		expand("align/{sample}.bam.bai",
               		sample=config["samples"]),
		expand("sniffles_genotypes_temp/{sample}.vcf",
			sample=config["samples"])		

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
		"align/{sample}.bam"
	log:
		"align/{sample}.log"
	conda: "envs/ngmlr.yaml"
	shell:
		"zcat {input.fq} | ngmlr --presets ont -t 22 -r {input.REF} | \
		 samtools sort -@ 8 -o {output} - 2> {log}"

rule samtools_index:
	input:
		"align/{sample}.bam"
	output:
		"align/{sample}.bam.bai"
	log:
		"align/samtools_index/{sample}.log"
	conda:
		"envs/samtools.yaml"
	shell:
		"samtools index -@ 8 {input} 2> {log}"


rule sniffles_call:
	input:
		bam = "align/{sample}.bam",
		bai = "align/{sample}.bam.bai"
	output:
		"sniffles_calls/{sample}.vcf"
	log:
		"sniffles_calls/{sample}.log"
	conda:
		"envs/sniffles.yaml"
	shell:
		"sniffles --mapped_reads {input.bam} --vcf {output} --threads 22  2> {log}"

rule survivor:
	input:
		[f"sniffles_calls/{sample}.vcf" for sample in config["samples"]]
#		 expand("sniffles_calls/{sample}.vcf",
#                       sample=config["samples"]),
	output:
		vcf = "sniffles_combined/calls.vcf",
		fofn = "sniffles_combined/samples.fofn"
	params:
		distance = config["parameters"]["survivor_distance"],
		caller_support = 1,
		same_type = 1,
		same_strand = -1,
		estimate_distance = -1,
		minimum_size = -1,
	conda: "envs/survivor.yaml"
	log:
		"sniffles_combined/calls.vcf"
	shell:
		"ls {input} > {output.fofn} ; \
		SURVIVOR merge {output.fofn} {params.distance} {params.caller_support} \
		{params.same_type} {params.same_strand} {params.estimate_distance}  \
		{params.minimum_size} {output.vcf} 2> {log}"

rule sniffles_genotype:
	input:
		bam = "align/{sample}.bam",
		ivcf = "sniffles_combined/calls.vcf"
	output:
		"sniffles_genotypes_temp/{sample}.vcf"
	conda: "envs/sniffles.yaml"
	log:
		"sniffles_genotype/{sample}.log"
	shell:
		"sniffles --mapped_reads {input.bam} \
			--vcf {output} \
			--threads 22 \
			--cluster \
			--Ivcf {input.ivcf} 2> {log}"

 
