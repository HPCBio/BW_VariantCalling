#!/bin/bash
#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

module load samtools/1.3.1
# Extract the required reference seqenceL
samtools faidx /home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome//HG19_GATKbundle2.8_noDecoys.fa chr1 chr2 chr3 > /home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/small_ref/HG19_GATKbundle2.8_chr1_chr2_chr3.fa

module load bwa/0.7.15
# Index the new reference
bwa index -a bwtsw /home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/small_ref/HG19_GATKbundle2.8_chr1_chr2_chr3.fa

samtools faidx /home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/small_ref/HG19_GATKbundle2.8_chr1_chr2_chr3.fa

module load picard-tools/2.4.1
# Create the dictionary
picard=/home/apps/picard-tools/picard-tools-2.4.1/picard.jar
java -jar $picard CreateSequenceDictionary \
	REFERENCE=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/small_ref/HG19_GATKbundle2.8_chr1_chr2_chr3.fa \
	OUTPUT=/home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/data/genome/small_ref/HG19_GATKbundle2.8_chr1_chr2_chr3.dict



