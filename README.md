Using Nanopore reads for variant calling. 

Largely based on https://github.com/wdecoster/nano-snakemake and https://github.com/nanoporetech/ont_tutorial_sv

1. Map reads with NGLMR.
1. Call SV with Sniffles (min read support set low (3) for my low coverage data).
1. Use Survivor to merge per sample VCFs (paramters for when to merge variants are important).
1. Force Sniffles to call genotype at the site of each variant in the merged file.

Conda is used to manage environments on a per rule basis, so be sure to deploy Snakemake with the --use-conda flag.

![DAG](https://github.com/JamieCFreeman/nglmr_SV/blob/main/README_files/rulegraph.svg?raw=true)
 


