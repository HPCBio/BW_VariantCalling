#!/bin/bash
#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/log.neat_comp.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/log.neat_comp.e 
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
module load python

set -x
neatdir=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/

########################### Comparison stage:
python $neatdir/vcf_compare_OLD.py \
  -r $referencedir/HG19_GATKbundle2.8_noDecoys.fa \
  -g $goldenFile \
  -w $workflowFile \
  -o $output/variant_comparison_neat_vcf_compare \
  -t $targeted_region \
  --vcf-out \
  --no-plot 
