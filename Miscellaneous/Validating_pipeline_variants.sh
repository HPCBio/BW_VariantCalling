#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/log.hap.py.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/log.hap.py.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
goldenFile=/home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf
workflowFile=/home/groups/hpcbio_shared/azza/GIAB/results/run8/delivery/jointVCFs/jointVCFcalled.vcf

export HGREF=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome/ucsc.hg19.fasta 

########################### Preparatory stages:
module load tabix
bgzip -c $goldenFile > $goldenFile.gz
tabix -p vcf $goldenFile.gz

bgzip -c $workflowFile > $workflowFile.gz
tabix -p vcf $workflowFile.gz

########################### Comparison stage:
module load hap.py/0.3.0
hap.py $goldenFile $workflowFile -o /home/groups/hpcbio_shared/azza/GIAB/results/run8/variant_compare_hap.py \
 --logfile /home/groups/hpcbio_shared/azza/GIAB/results/run8/hap.py.log \
 --no-internal-leftshift \
  --no-internal-preprocessing \
  --include-nonpass-truth -V \
  --roc Q_GQ \
  --roc-filter LowGQX
