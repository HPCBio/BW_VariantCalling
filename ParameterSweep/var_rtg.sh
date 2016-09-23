#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/src/rtg.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/src/rtg.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
goldenFile=/home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf
workflowFile=/home/groups/hpcbio_shared/azza/GIAB/results/run8/delivery/jointVCFs/jointVCFcalled.vcf

export HGREF=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome/ucsc.hg19.fasta 

########################### Needed tools and preps:
module load rtg-tools/3.6.2  
rtg format -o $reference/ucsc.hg19.sdf $reference/ucsc.hg19.fasta

###########################
rtg vcfeval -b $goldenFile -c $workflowFile -o /home/groups/hpcbio_shared/azza/GIAB/results/run8/variant_compare_rtg -t $reference/ucsc.hg19.sdf

