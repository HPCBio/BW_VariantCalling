#!/bin/bash
#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/log.bcftools.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/log.bcftools.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
set -x
output=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/results/run8
goldenFile=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/reads/synthetic.Oct2016.WES30x/sample.1/Merged/synthetic.Oct2016.WES30x_merged_golden.vcf
workflowFile=$output/delivery/jointVCFs/jointVCFcalled.vcf
targeted_region=/home/groups/hpcbio_shared/azza/TargetedRegions-Azza-has-permission/ZachUses_Illumina_truseq_exome_targeted_regions.hg19.chr.bed

########################### Preparatory stages:
set +x
module load tabix

set -x
bgzip -c $goldenFile > $goldenFile.gz
tabix -p vcf $goldenFile.gz

bgzip -c $workflowFile > $workflowFile.gz
tabix -p vcf $workflowFile.gz

########################### Comparison stage:

set -x

bcftools stats $workflowFile.gz $goldenFile.gz -R $targeted_region > $output/variant_comparison_bcftools

plot-vcfstats -p comps $output/variant_comparison_bcftools

