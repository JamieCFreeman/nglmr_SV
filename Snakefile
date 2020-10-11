
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
	conda: "../envs/ngmlr.yaml"
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
		"../envs/samtools.yaml"
	shell:
		"samtools index -@ 8 {input} 2> {log}"






 
