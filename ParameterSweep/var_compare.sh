#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o /home/groups/hpcbio_shared/azza/GIAB/neat.ou
#PBS -e /home/groups/hpcbio_shared/azza/GIAB/neat.e 
#PBS -M aeahmed@illinois.edu
#PBS -m abe

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
goldenFile=/home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf
workflowFile=/home/groups/hpcbio_shared/azza/GIAB/results/run8/delivery/jointVCFs/jointVCFcalled.vcf

module load python/2.7.9
vcf_compare_dir=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/builds/NEAT/neat-genreads/utilities

python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g ${goldenFile} -w ${workflowFile}  -o /home/groups/hpcbio_shared/azza/GIAB/results/run8/variant_compare_neat  --incl-homs --incl-fail --vcf-out --no-plot

