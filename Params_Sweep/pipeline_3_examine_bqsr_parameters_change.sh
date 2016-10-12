#!/bin/bash

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs-bqsr.log.txt
#PBS -e localhost:$HOME/errors-bqsr.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m abe

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"
	set -x
######### Paths defintions:
reference="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"
reads="/home/groups/hpcbio_shared/azza/GIAB/reads"
results="/home/groups/hpcbio_shared/azza/GIAB/results"

# Default run
if [ ! -d  $results/p4_indelRealign_bqsr/default ]; then
	mkdir -p  $results/p4_indelRealign_bqsr/bqsr/default
fi

cd $results/p4_indelRealign_bqsr/bqsr/default

module load gatk/3.6
GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar

module load R/3.2.3

echo Run using the default parameters: > "../../bqsr.summary.txt"
echo "The parameters are:     ( ics maxCycle mcs bqsrBAQGOP ddq idq lqt mdq )">> "../../bqsr.summary.txt"
echo "The default values are: ( 3     500     2      40     45   45  2  -1  )">> "../../bqsr.summary.txt"
echo  "------------------------------------------------------------------------------------------------------------" >> "../../bqsr.summary.txt"

#if [ $run_default ]; then
	START=$(date +%s)
	java -jar $GenomeAnalysisTK\
		-T BaseRecalibrator\
		-R $reference/ucsc.hg19.fasta\
		-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
		-L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		-knownSites $reference/dbsnp_138.hg19.vcf\
	  	-knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
	 	-knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
		-o recal.table.default.0
	END=$(date +%s)
	DIFF=$(( $END - $START ))
	exitcode=$?
	echo "exit code from BaseRecalibrator using: defaults 0 = $exitcode ; execution time = $DIFF" >> "../../bqsr.summary.txt"
	
	START=$(date +%s)
	java -jar $GenomeAnalysisTK\
	       	-T BaseRecalibrator\
	        -R $reference/ucsc.hg19.fasta\
		-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
	        -L $reads/TruSeq_exome_targeted_regions.hg19.bed\
	        -knownSites $reference/dbsnp_138.hg19.vcf\
	        -knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
	        -knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
	        -BQSR recal.table.default.0\
	        -o after_recal.table.default.0
	END=$(date +%s)
	DIFF=$(( $END - $START ))
	exitcode=$?
	echo "exit code from the rerun of BaseRecalibrator using: defaults 0 = $exitcode  ; execution time = $DIFF" >> "../../bqsr.summary.txt"

	START=$(date +%s)
	java -jar $GenomeAnalysisTK\
		-T AnalyzeCovariates\
		-R $reference/ucsc.hg19.fasta\
		-L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		-before recal.table.default.0\
		-after after_recal.table.default.0\
		-plots recalibration_plots.default.0.pdf
	END=$(date +%s)
	DIFF=$(( $END - $START ))
	exitcode=$?
	echo "exit code from AnalyzeCovariates using: defaults 0 = $exitcode  ; execution time = $DIFF" >> "../../bqsr.summary.txt"
	
	
	START=$(date +%s)
	java -jar $GenomeAnalysisTK\
		-T PrintReads\
		-R $reference/ucsc.hg19.fasta\
		-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
		-L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		-BQSR recal.table.default.0\
		-o recal.default.0.bam
	END=$(date +%s)
	DIFF=$(( $END - $START ))
	exitcode=$?
	echo "exit code from PrintReads using: defaults 0 = $exitcode  ; execution time = $DIFF" >> "../../bqsr.summary.txt"
	
	echo >> "../../bqsr.summary.txt"
	echo "####################################################################################################" >> ../../bqsr.summary.txt
#fi

declare -a parameters=(ics maxCycle mcs bqsrBAQGOP ddq idq lqt mdq)
declare -a min=(1 250 1 10 10 10 1 2)
declare -a step=(2 150 2 10 10 10 2 2)
declare -a max=(13 1000 13 70 70 70 12 12) 


cd $results/p4_indelRealign_bqsr/bqsr
mkdir ${parameters[@]}

echo The parameters being tested and their ranges are given below: >> "../bqsr.summary.txt"
echo paramters: ${parameters[@]}, >> "../bqsr.summary.txt"
echo minimum  : ${min[@]} >> "../bqsr.summary.txt"
echo maximum  : ${max[@]} >> "../bqsr.summary.txt"
echo >> "../bqsr.summary.txt"
echo  "------------------------------------------------------------------------------------------------------------" >> ../"bqsr.summary.txt"

pos=0
echo start while loop ; 

while [ $pos -lt ${#parameters[@]} ]
do
	echo inside while
	module load gatk/3.6
	GenomeAnalysisTK=/home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar
	module load R/3.2.3
	par=${parameters[pos]}
	cd $results/p4_indelRealign_bqsr/bqsr/$par
	for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
		START=$(date +%s)
		java -jar $GenomeAnalysisTK\
		        -T BaseRecalibrator\
		        -R $reference/ucsc.hg19.fasta\
			-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
		        -L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		        -knownSites $reference/dbsnp_138.hg19.vcf\
		        -knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
		        -knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
			-$par "$i"\
		        -o recal.table.$par.$i
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
		echo "exit code from BaseRecalibrator using: -$par $i = $exitcode ; execution time = $DIFF" >> "../bqsr.summary.txt"

		if [ -s "recal.table.$par.$i" ]; then
			echo 'successful parameter test $par=$i'
		else
			echo "BQSR failed using -$par $i">> "../bqsr.summary.txt"
			echo "Execution time is :$DIFF: seconds" >> "../bqsr.summary.txt"
			echo
			continue
		fi
			
		START=$(date +%s)
		java -jar $GenomeAnalysisTK\
		        -T BaseRecalibrator\
		        -R $reference/ucsc.hg19.fasta\
			-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
		        -L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		        -knownSites $reference/dbsnp_138.hg19.vcf\
		        -knownSites $reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf\
		        -knownSites $reference/1000G_phase1.indels.hg19.sites.vcf\
		        -BQSR recal.table.$par.$i\
		        -o after_recal.table.$par.$i
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
	        echo "exit code from the rerun of  BaseRecalibrator with: -$par $i = $exitcode ; execution time = $DIFF" >> "../bqsr.summary.txt"
		
		START=$(date +%s)
		java -jar $GenomeAnalysisTK\
		        -T AnalyzeCovariates\
		        -R $reference/ucsc.hg19.fasta\
		        -L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		        -before recal.table.$par.$i\
		        -after after_recal.table.$par.$i\
		        -plots recalibration_plots.$par.$i.pdf
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
	        echo "exit code from AnalyzeCovariates using $par $i = $exitcode; execution time = $DIFF" >> "../bqsr.summary.txt"

		START=$(date +%s)
		java -jar $GenomeAnalysisTK\
		        -T PrintReads\
		        -R $reference/ucsc.hg19.fasta\
			-I $reads/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam\
		        -L $reads/TruSeq_exome_targeted_regions.hg19.bed\
		        -BQSR recal.table.$par.$i\
		        -o recal.$par.$i.bam
		END=$(date +%s)
		DIFF=$(( $END -$START + $DIFF ))
		exitcode=$?
		echo "exit code from PrintReads with: -$par $i = $exitcode  ; execution time = $DIFF" >> "../bqsr.summary.txt"
		echo  >> "../bqsr.summary.txt"
	done
let pos+=1
done



echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
