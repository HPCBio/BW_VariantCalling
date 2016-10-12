#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs-hap.py_VariantEval.log.txt
#PBS -e localhost:$HOME/errors-hap.py_VariantEval.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
reads="/home/groups/hpcbio_shared/azza/GIAB/reads"
results="/home/groups/hpcbio_shared/azza/GIAB/results"

TargetedRegions="/home/groups/hpcbio_shared/azza/TargetedRegions"

#### Needed modules:

module load hap.py/0.3.0

### hap.py prefers vcf.gz files, so, if needed:
module load tabix
bgzip -c file.vcf > file.vcf.gz
tabix -p vcf file.vcf.gz

# The objectives of this exercise:
## 1) See the effect of different preprocessig steps in the called variants
## 2) See how the different programs compare variants, and understand which is approapriate in this study. What is special here is that were are particularily interested in rare variants, and so the comparison needs to be very exact. In and ideal world, the comparison would further tell me where the 2 files agree and disagree (besides the stats on sensitivity and specifity)

# ######################################################################################################################## #

if [ ! -d $results/p5_Haplotypecaller/comparison_stages/hap.py ]; then
        mkdir -p $results/p5_Haplotypecaller/comparison_stages/hap.py
fi

cd $results/p5_Haplotypecaller/comparison_stages/hap.py

##
#2. Using hap.py 

### 1.
START=$(date +%s)
hap.py \
 $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf \
 $results/p5_Haplotypecaller/vars_default_alinged_deduped.vcf \
 -r $reference/ucsc.hg19.fasta \
 -o default_aligned_deduped_hap.py
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from illumina, hap.py = $exitcode; Runtime =$DIFF (aligned deduped)" >> "../variants_comparison_summary_hap.py.txt"

### 2.
START=$(date +%s)
hap.py \
 $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf \
 $results/p5_Haplotypecaller/vars_ics7_best_bqsr.vcf \
 -r $reference/ucsc.hg19.fasta \
 -o best_bqsr_ics7_hap.py 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from illumina, hap.py = $exitcode; Runtime =$DIFF (best BQSR ics7)" >> "../variants_comparison_summary_hap.py.txt"

## 3.
START=$(date +%s)
hap.py \
 $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf \
 $results/p5_Haplotypecaller/vars_bqsrBAQGOP20_worst_bqsr.vcf \ 
 -r $reference/ucsc.hg19.fasta \
 -o worst_bqsr_bqsrBAQGOP20_hap.py 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from illumina, hap.py = $exitcode; Runtime =$DIFF (worst BQSR bqsrBAQGOP20)" >> "../variants_comparison_summary_hap.py.txt"

## 3b.
START=$(date +%s)
hap.py \
 $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf \
 $results/p5_Haplotypecaller/vars_bqsr_default_bqsr.vcf \
 -r $reference/ucsc.hg19.fasta \
 -o default_bqsr_hap.py 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from illumina, hap.py = $exitcode; Runtime =$DIFF (default bqsr)" >> "../variants_comparison_summary_hap.py.txt"

## 4.
START=$(date +%s)
hap.py \
 $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf \
 $results/p5_Haplotypecaller/vars_default_aligned.vcf \
 -r $reference/ucsc.hg19.fasta \
 -o default_aligned_hap.py 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from illumina, hap.py = $exitcode; Runtime =$DIFF (default aligment)" >> "../variants_comparison_summary_hap.py.txt"

## 5.
START=$(date +%s)
hap.py \
 $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf \
 $results/p5_Haplotypecaller/vars_worst_aligned_U4.vcf \
 -r $reference/ucsc.hg19.fasta \
 -o worst_aligned_U4_hap.py 
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from illumina, hap.py = $exitcode; Runtime =$DIFF (worst aligment U4)" >> "../variants_comparison_summary_hap.py.txt"


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
