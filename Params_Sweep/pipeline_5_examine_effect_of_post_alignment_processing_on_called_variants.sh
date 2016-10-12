#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs-hc.different.files.log.txt
#PBS -e localhost:$HOME/errors-hc.different.files.log.txt		
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

module load novocraft

# HC needs to be called for 3 different files. While at it, let's compare the run time in each case. For fun, let's do it also for the file that is just aligned, just realigned, in addition to the best scoring files from BQSR. These are: 
# 1. Aligned deduped bam 
# 2. 'Best' bqsr file
# 3. 'Worst' bqsr file
# 3.b default bqsr
# 4. 'default' aligned file
# 5. 'Worst' aligned file
# Note that the last couple of files will need their read group info added first!

if [ ! -d $results/p5_Haplotypecaller ]; then
	mkdir -p $results/p5_Haplotypecaller
fi

cd $results/p5_Haplotypecaller

if [ $done ]; then
# 1. Aligned Deduped file
START=$(date +%s)
java -jar $GenomeAnalysisTK\
	-T HaplotypeCaller\
	-R $reference/ucsc.hg19.fasta\
	-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
	--dbsnp $reference/dbsnp_138.hg19.vcf\
	-o vars_default_alinged_deduped.g.vcf\
	-ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (alnded default)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_default_alinged_deduped.g.vcf\
        -o vars_default_alinged_deduped.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (alnded default)" >> "hc.summary-various.stages.txt"

# 2. Best bqsr file 
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/bqsr/ics/recal.ics.7.bam\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_ics7_best_bqsr.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (best BQSR ics7)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_ics7_best_bqsr.g.vcf\
        -o vars_ics7_best_bqsr.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (best BQSR ics7)" >> "hc.summary-various.stages.txt"

# 3. Worst bqsr file 
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/bqsr/bqsrBAQGOP/recal.bqsrBAQGOP.20.bam\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_bqsrBAQGOP20_worst_bqsr.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (worst BQSR bqsrBAGOP20)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_bqsrBAQGOP20_worst_bqsr.g.vcf\
        -o vars_bqsrBAQGOP20_worst_bqsr.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (worst BQSR bqsrBAGOP20)" >> "hc.summary-various.stages.txt"

####
# 3b. default bqsr file 
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p4_indelRealign_bqsr/bqsr/default/recal.default.0.bam\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_bqsr_default_bqsr.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (default BQSR )" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_bqsr_default_bqsr.g.vcf\
        -o vars_bqsr_default_bqsr.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (default BQSR )" >> "hc.summary-various.stages.txt"

fi #### these are processed successfully, let's address the remaining ones:
####

# rememeber? GATK only works with indexed bams!
novosort $results/p2_alignment/default/a.default.0.bam -i -o $results/p2_alignment/default/a.default.0.indexed.bam


# 4. Default alignment file
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I a.default.0.indexed.bam\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_default_aligned.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (aligned default)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_default_aligned.g.vcf\
        -o vars_default_aligned.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (aligned default)" >> "hc.summary-various.stages.txt"

### indexing this as well:
novosort $results/p2_alignment/U/a.U.4.bam -i -o $results/p2_alignment/U/a.U.4.indexed.bam

# 5. Worst alignment file
START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T HaplotypeCaller\
        -R $reference/ucsc.hg19.fasta\
        -I $results/p2_alignment/U/a.U.4.indexed.bam\
        --dbsnp $reference/dbsnp_138.hg19.vcf\
        -o vars_worst_aligned_U4.g.vcf\
        -ERC GVCF
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from HaplotypeCalle = $exitcode; Runtime =$DIFF  (worst aligned U4)" >> "hc.summary-various.stages.txt"

START=$(date +%s)
java -jar $GenomeAnalysisTK\
        -T GenotypeGVCFs\
        -R $reference/ucsc.hg19.fasta\
        -V vars_worst_aligned_U4.g.vcf\
        -o vars_worst_aligned_U4.vcf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from GenotypeGVCFs = $exitcode; Runtime =$DIFF (worst aligned U4)" >> "hc.summary-various.stages.txt"

echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
