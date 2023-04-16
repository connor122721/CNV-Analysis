#!/usr/bin/env bash

# Modules to load
module load anaconda/2020.11-py3.8
source activate msprime_env

# Working & temp directory
wd=/project/berglandlab/connor/BACKUP_scratch/all_bam/Euro_bams
out=/project/berglandlab/connor/control_cnvs

# April_2017_DBunk_13_finalmap_mdup.bam; Scaffold_9201_HRSCAF_10758_70001_200000_CN1

# Run script
samplot plot \
-b ${wd}/April_2017_DBunk_13_finalmap_mdup.bam \
-c Scaffold_9201_HRSCAF_10758 \
-s 1 \
-e 2000000 \
--coverage_only \
--max_coverage 100 \
-o ${out}/April_2017_DBunk_13_CN1_chitin

# Finish
echo "Finish"
echo date
