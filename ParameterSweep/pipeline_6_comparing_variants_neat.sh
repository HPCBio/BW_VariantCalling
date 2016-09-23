#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs-neat_var_comp.log.txt
#PBS -e localhost:$HOME/errors-neat_var_comp.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
reads="/home/groups/hpcbio_shared/azza/GIAB/reads"
results="/home/groups/hpcbio_shared/azza/GIAB/results"

#### Needed modules:

module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

module load R/3.2.3

# HC was called for 5 different files, and it would be nice to compare the variants called in each case. These files are:
# 1. $results/p5_Haplotypecaller 
# 2. $results/p5_Haplotypecaller
# 3. $results/p5_Haplotypecaller
# 4. $results/p5_Haplotypecaller
# 5. $results/p5_Haplotypecaller
# The objectives of this exercise:
## 1) See the effect of different preprocessig steps in the called variants
## 2) See how the different programs compare variants, and understand which is approapriate in this study. What is special here is that were are particularily interested in rare variants, and so the comparison needs to be very exact. In and ideal world, the comparison would further tell me where the 2 files agree and disagree (besides the stats on sensitivity and specifity)

# ######################################################################################################################## #
# The following tools should be tested:
# 1. vcf_compare from UIUC
# 2. varEval from GATK
# 3. hap.py from Illumina (https://github.com/Illumina/hap.py/tree/ee33640b218b5c575e999be619884992781ab9ac)
# 4. rtg-tools from RealTimeGenomics (https://github.com/RealTimeGenomics/rtg-tools/tree/fdeeb3554badf1286fffc413d57913a9b7dfaa46) 
#               note that the latest updates on hap.py and rtg-tools date back to 2015 (Sept, and Oct respectively)
# ######################################################################################################################## #

if [ ! -d $results/p5_Haplotypecaller/comparison_stages/ ]; then
	mkdir -p $results/p5_Haplotypecaller/comparison_stages/
fi

cd $results/p5_Haplotypecaller/comparison_stages/

##
#1. Using vcf_compare.py from UIUC

module load python/2.7.9
vcf_compare_dir=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/builds/NEAT/neat-genreads/utilities

# 1.
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf -w $results/p5_Haplotypecaller/vars_default_alinged_deduped.vcf -o default_aligned_deduped_neat $results/tmp --incl-homs --incl-fail  -t $reads/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed -T 90 --vcf-out --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (aligned deduped)" >> "../variants_comparison.summary.txt"

### 2.
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf -w $results/p5_Haplotypecaller/vars_ics7_best_bqsr.vcf -o best_bqsr_ics7_neat $results/tmp --incl-homs --incl-fail  -t $reads/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed -T 90 --vcf-out --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (best BQSR ics7)" >> "../variants_comparison.summary.txt"

## 3.
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf -w $results/p5_Haplotypecaller/vars_bqsrBAQGOP20_worst_bqsr.vcf -o worst_bqsr_bqsrBAQGOP20_neat $results/tmp --incl-homs --incl-fail  -t $reads/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed -T 90 --vcf-out --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (worst BQSR bqsrBAQGOP20)" >> "../variants_comparison.summary.txt"

## 3b.
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf -w $results/p5_Haplotypecaller/vars_bqsr_default_bqsr.vcf -o default_bqsr_neat $results/tmp --incl-homs --incl-fail  -t $reads/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed -T 90 --vcf-out --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (default bqsr)" >> "../variants_comparison.summary.txt"

## 4.
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf -w $results/p5_Haplotypecaller/vars_default_aligned.vcf -o default_aligned_neat $results/tmp --incl-homs --incl-fail  -t $reads/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed -T 90 --vcf-out --no-plot

END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (default aligment)" >> "../variants_comparison.summary.txt"

## 5.
START=$(date +%s)
python $vcf_compare_dir/vcf_compare_OLD.py -r $reference/ucsc.hg19.fasta -g $reads/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf -w $results/p5_Haplotypecaller/vars_worst_aligned_U4.vcf -o worst_aligned_U4_neat $results/tmp --incl-homs --incl-fail  -t $reads/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed -T 90 --vcf-out --no-plot
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from Neat, var_comp.py = $exitcode; Runtime =$DIFF (worst aligment U4)" >> "../variants_comparison.summary.txt"



####################################################################### generated vcfs
# 1. $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam
#$results/p5_Haplotypecaller/vars_default_bqsr.vcf

# 2. $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.5.bam
#$results/p5_Haplotypecaller/vars_ics5_bqsr.vcf

# 3. $results/p4_indelRealign_bqsr/bqsr/default/recal.ics.7.bam
#$results/p5_Haplotypecaller/vars_ics7_bqsr.vcf

# 4. $results/p4_indelRealign_bqsr/realignedBam.default.0.bam
#$results/p5_Haplotypecaller/vars_default_realigned.vcf

# 5. $results/p4_indelRealign_bqsr/markduplicates_RG_added2.bam
#$results/p5_Haplotypecaller/vars_default_markDuplicates.vcf


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
