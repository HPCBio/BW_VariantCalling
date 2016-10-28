#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/log.hap.py.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/log.hap.py.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
output=/home/groups/hpcbio_shared/azza/GIAB/results/run8/
goldenFile=/home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf
workflowFile=$output/delivery/jointVCFs/jointVCFcalled.vcf

referencedir=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome
########################### Preparatory stages:
module load python
set -x
neatdir=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/

########################### Comparison stage:
python $neatdir/vcf_compare_OLD.py \
  -r $referencedir/ucsc.hg19.fasta \
  -g $goldenFile \
  -w $workflowFile \
  -o $output/variant_comparison_neat_vcf_compare \
  --vcf-out \
  --no-plot \
  -T 90
