#!/bin/bash 
#	
#  script to perform variant calling with HaplotypeCaller ONLY on a region of the sample
#  This module is called from within the realign module
######################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu

set -x
echo `date`
ulimit -s unlimited
umask 0027
scriptfile=$0
inputfile=$1
outputfile=$2
outputdir=$3
target=$4
chr=$5
runfile=$6
failedlog=$7
email=$8
qsubfile=$9

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$failedlog\noutputlog=$failedlog"


set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;
       
if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
   exit 1;
fi

rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
snvcaller=$( cat $runfile | grep -w SNV_CALLER | cut -d '=' -f2 )
skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
dbsnp=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
genderinfo=$( cat $runfile | grep -w GENDERINFORMATION | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
variantcmd=$( cat $runfile | grep -w VARIANT_CMD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
variantAnalysis=$( cat $runfile | grep -w VARIANT_ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
tabixdir=$( cat $runfile | grep -w TABIXDIR | cut -d '=' -f2 )


set +x; echo -e "\n\n########### should we skip variant calling     #############\n\n" >&2; set -x;

if [ $skipvcall == "1" -o $skipvcall == "YES" ]
then
	echo "skipping the execution of this variant calling module"
	exit 0;
fi

set +x; echo -e "\n\n########### checking input type     #############\n\n" >&2; set -x;

if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
    input_type="WGS"
fi

if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
    input_type="WES"
fi

set +x; echo -e "\n\n########### checking variant calling analysis     #############\n\n" >&2; set -x;

if [ $variantAnalysis == "GENOTYPE" -o $variantAnalysis == "GENOTYPING" ]
then
	variantAnalysis="GENOTYPING"
fi

set +x; echo -e "\n\n########### checking variant calling command     #############\n\n" >&2; set -x;

if [ $variantcmd != "HC" -a $variantcmd != "HAPLOTYPECALLER"  ]
then
	MSG="VARIANT_CMD=$variantcmd this case is not analyzed in this script"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
elif [ $variantcmd == "HC" -o $variantcmd == "HAPLOTYPECALLER"  ]
then
	variantcmd="HC"
fi

set +x; echo -e "\n\n########### checking variant calling toolkit    #############\n\n" >&2; set -x;

if [ $snvcaller != "GATK" ]
then
	MSG="$snvcaller snvcaller tool is not used in this script"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi


set +x; echo -e "\n\n########### checking inputfile     #############\n\n" >&2; set -x;

if [ ! -s $inputfile ]
then
	MSG="$inputfile realigned bam file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

set +x; echo -e "\n\n########### checking callsets     #############\n\n" >&2; set -x;

if [ ! -d $refdir ]
then
	MSG="$refdir reference genome directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi      
if [ ! -s $refdir/$ref ]
then
	MSG="$ref reference genome not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi
if [ ! -s $refdir/$dbsnp ]
then
	MSG="$refdir/$dbsnp dbSNP for reference genome not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

set +x; echo -e "\n\n########### checking tool directories     #############\n\n" >&2; set -x;

if [ ! -d $samdir ]
then
	MSG="$samdir samtools directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi
if [ -z $javadir ]
then
	MSG="A value must be specified for JAVADIR in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi      
if [ ! -d $gatk ]
then
	MSG="$gatk GATK directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

if [ ! -d $tabixdir ]
then
	MSG="$tabixdir tabix directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

if [ ! -d $outputdir ]
then
    mkdir -p $outputdir
fi

echo `date`

#tabix needs to be on your path for GATK to produce *tbi files
export PATH=${PATH}:$tabixdir

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  PREPARATORY WORK                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $outputdir/..
sample=`basename $PWD`
cd $$outputdir
intervals=" "                                    # for file with WES intervals list, if it is provided

set +x; echo -e "\n\n########### defining region     #############\n\n" >&2; set -x;

if [ $input_type == "WES" -a $target != "NOTARGET" -a -s $target ]
then
	intervals=" -L $target "
else
	intervals=" -L $chr "
fi


set +x; echo "############################################################################################" >&2
echo "############################################################################################" >&2
echo "#############  HaplotypeCaller cases: with or without gender information  ###############" >&2; 
echo "############################################################################################" >&2
set -x;       

if [ $input_type == "WGS" -a  -s $genderinfo ]
then
	set +x; echo -e "\n\n" >&2;
	echo "############################################################################################" >&2       
	echo "############################################################################################" >&2
	echo "##########   HaplotypeCaller   with GENDERINFO file AND WGS input data           ###########" >&2
	echo "############################################################################################" >&2        
	echo -e "\n\n" >&2; set -x;

	echo `date`

	set +x; echo -e "\n\n#############    gathering gender and sample info  ###############\n\n" >&2; set -x;


	sample_id=$sample
	gender=`cat $genderinfo | grep -v "#" | grep $sample_id | awk -F$'\t' '{print $2}'` # Can take this line out later
	ploidy=2
	echo "Sample: $sample_id, gender: $gender, chr: $chr, default ploidy:2" # Can take this line out later                                                                 


	set +x; echo -e "\n\n#############    var-calling conditions start here  ###############\n\n" >&2; set -x;        


	set +x; echo -e "\n\n#############    Calling chr X   ###############\n\n" >&2; set -x;  
       
	if [ $chr == "X" ];
	then
       
		set +x; echo -e "\n\n#############   calling chr=X with for TWO genders and THREE regions X_PAR1 X_PAR2 X_nonPAR ###############\n\n" >&2; set -x;        

		# The X PAR coordinates can be found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt
		x_par1="X:60001-2699520"
		x_par2="X:154931044-155260560"

		### Calling male X
		if [ $gender == "Male" ];
		then
 
			set +x; echo -e "\n\n#############   calling chr=X gender=male X_PAR1  ###############\n\n" >&2; set -x;        

			site="-L "$x_par1
			site_name="X_PAR1"
			outfile=$sample.$site_name.raw.g.vcf.gz

			java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $refdir/$ref \
			-I ${inputfile} \
			--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
			-gt_mode DISCOVERY \
			-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
			-stand_call_conf 30 \
			-stand_emit_conf 30 \
			--sample_ploidy $ploidy \
			-nt 1 -nct $thr  \
			--dbsnp $refdir/$dbsnp  \
			$site \
			-o $outfile

			exitcode=$?

			echo `date`		
			if [ $exitcode -ne 0 ]
			then
				MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
				exit $exitcode;
			fi
			if [ ! -s $outfile ] 
			then
				MSG="$outfile HaplotypeCaller file not created. vcall failed."
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
				exit 1;
			fi           

			set +x; echo -e "\n\n#############   calling chr=X gender=male X_PAR2  ###############\n\n" >&2; set -x;        


			site="-L "$x_par2
			site_name="X_PAR2"
			outfile=$sample.$site_name.raw.g.vcf.gz

			java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $refdir/$ref \
			-I ${inputfile} \
			--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
			-gt_mode DISCOVERY \
			-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
			-stand_call_conf 30 \
			-stand_emit_conf 30 \
			--sample_ploidy $ploidy \
			-nt 1 -nct $thr  \
			--dbsnp $refdir/$dbsnp  \
			$site \
			-o $outfile

			exitcode=$?

			echo `date`		
			if [ $exitcode -ne 0 ]
			then
				MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
				exit $exitcode;
			fi
			if [ ! -s $outfile ] 
			then
				MSG="$outfile HaplotypeCaller file not created. vcall failed."
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
				exit 1;
			fi           

			set +x; echo -e "\n\n#############   calling chr=X gender=male X_nonPAR  ###############\n\n" >&2; set -x;        


			ploidy=1
			site="-L X -XL "$x_par1" -XL "$x_par2
			site_name="X_nonPAR"
			outfile=$sample.$site_name.raw.g.vcf.gz

			java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $refdir/$ref \
			-I ${inputfile} \
			--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
			-gt_mode DISCOVERY \
			-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
			-stand_call_conf 30 \
			-stand_emit_conf 30 \
			--sample_ploidy $ploidy \
			-nt 1 -nct $thr  \
			--dbsnp $refdir/$dbsnp  \
			$site \
			-o $outfile

			exitcode=$?

			echo `date`		
			if [ $exitcode -ne 0 ]
			then
				MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
				exit $exitcode;
			fi
			if [ ! -s $outfile ] 
			then
				MSG="$outfile HaplotypeCaller file not created. vcall failed."
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
				exit 1;
			fi           



		fi ### End calling male on chr X  
		echo `date`

		### Calling female X
		if [ $gender == "Female" ];
		then
 
			set +x; echo -e "\n\n#############   calling chr=X gender=female   ###############\n\n" >&2; set -x;        

			ploidy=2
			outfile=$sample.$chr.raw.g.vcf.gz

			java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $refdir/$ref \
			-I ${inputfile} \
			--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
			-gt_mode DISCOVERY \
			-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
			-stand_call_conf 30 \
			-stand_emit_conf 30 \
			--sample_ploidy $ploidy \
			-nt 1 -nct $thr  \
			--dbsnp $refdir/$dbsnp  \
			-o $outfile

			exitcode=$?

			echo `date`		
			if [ $exitcode -ne 0 ]
			then
				MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
				exit $exitcode;
			fi
			if [ ! -s $outfile ] 
			then
				MSG="$outfile HaplotypeCaller file not created. vcall failed."
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
				exit 1;
			fi           
		fi  ### End calling female on chr X

 
		set +x; echo -e "\n\n#############   end calling chr=X   ###############\n\n" >&2; set -x;        
		echo `date`

	set +x; echo -e "\n\n#############    Calling chr Y   ###############\n\n" >&2; set -x;
	elif [ $chr == "Y" -a $gender == "Male" ]
	then
       
		 set +x; echo -e "\n\n#############   calling chr=Y with for ONE gender and THREE regions Y_PAR1 Y_PAR2 Y_nonPAR ###############\n\n" >&2; set -x;        

		 # The X PAR coordinates can be found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt
		 y_par1="Y:10001-2649520"
		 y_par2="Y:59034050-59363566"

 
		set +x; echo -e "\n\n#############   calling chr=Y gender=male Y_PAR1  ###############\n\n" >&2; set -x;        

		site="-L "$y_par1
		site_name="Y_PAR1"
		outfile=$sample.$site_name.raw.g.vcf.gz

		java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $refdir/$ref \
		-I ${inputfile} \
		--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		-gt_mode DISCOVERY \
		-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		-stand_call_conf 30 \
		-stand_emit_conf 30 \
		--sample_ploidy $ploidy \
		-nt 1 -nct $thr  \
		--dbsnp $refdir/$dbsnp  \
		$site \
		-o $outfile


		exitcode=$?

		echo `date`		
		if [ $exitcode -ne 0 ]
		then
			MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
			exit $exitcode;
		fi
		if [ ! -s $outfile ] 
		then
			MSG="$outfile HaplotypeCaller file not created. vcall failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
			exit 1;
		fi           


		set +x; echo -e "\n\n#############   calling chr=Y gender=male Y_PAR2  ###############\n\n" >&2; set -x;        

		site="-L "$y_par2
		site_name="Y_PAR2"
		outfile=$sample.$site_name.raw.g.vcf.gz

		java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $refdir/$ref \
		-I ${inputfile} \
		--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		-gt_mode DISCOVERY \
		-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		-stand_call_conf 30 \
		-stand_emit_conf 30 \
		--sample_ploidy $ploidy \
		-nt 1 -nct $thr  \
		--dbsnp $refdir/$dbsnp  \
		$site \
		-o $outfile

		exitcode=$?

		echo `date`		
		if [ $exitcode -ne 0 ]
		then
			MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
			exit $exitcode;
		fi
		if [ ! -s $outfile ] 
		then
			MSG="$outfile HaplotypeCaller file not created. vcall failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
			exit 1;
		fi           


		set +x; echo -e "\n\n#############   calling chr=Y gender=male Y_nonPAR  ###############\n\n" >&2; set -x;        

		ploidy=1
		site="-L Y -XL "$y_par1" -XL "$y_par2
		site_name="Y_nonPAR"
		outfile=$isample.$site_name.raw.g.vcf.gz

		java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $refdir/$ref \
		-I ${inputfile} \
		--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		-gt_mode DISCOVERY \
		-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		-stand_call_conf 30 \
		-stand_emit_conf 30 \
		--sample_ploidy $ploidy \
		-nt 1 -nct $thr  \
		--dbsnp $refdir/$dbsnp  \
		$site \
		-o $outfile

		exitcode=$?

		echo `date`		
		if [ $exitcode -ne 0 ]
		then
			MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
			exit $exitcode;
		fi
		if [ ! -s $outfile ] 
		then
			MSG="$outfile HaplotypeCaller file not created. vcall failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
			exit 1;
		fi   

		set +x; echo -e "\n\n#############   end calling chr=X   ###############\n\n" >&2; set -x;        
		echo `date`


	set +x; echo -e "\n\n#############    Calling chr Y   ###############\n\n" >&2; set -x;


	elif [ $chr == "MT" ];
	then
       
		set +x; echo -e "\n\n#############   calling chr=MT   ###############\n\n" >&2; set -x;        

		ploidy=1
		site_name=$chr
		site="-L $chr"
		outfile=$sample.$site_name.raw.g.vcf.gz

		java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $refdir/$ref \
		-I ${inputfile} \
		--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		-gt_mode DISCOVERY \
		-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		-stand_call_conf 30 \
		-stand_emit_conf 30 \
		--sample_ploidy $ploidy \
		-nt 1 -nct $thr  \
		--dbsnp $refdir/$dbsnp  \
		$site \
		-o $outfile

		exitcode=$?

		echo `date`		
		if [ $exitcode -ne 0 ]
		then
			MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
			exit $exitcode;
		fi
		if [ ! -s $outfile ] 
		then
			MSG="$outfile HaplotypeCaller file not created. vcall failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
			exit 1;
		fi           
	 
	 
		set +x; echo -e "\n\n#############   end calling chr=MT   ###############\n\n" >&2; set -x;        
		echo `date`

       else
       
		set +x; echo -e "\n\n#############  Calling the autosomes, chr-1-chr22    ###############\n\n" >&2; set -x;        

		echo "Call 1->22" # Can take this line out later
		ploidy=2
		outfile=$outputfile

		java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R $refdir/$ref \
		-I ${inputfile} \
		--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		-gt_mode DISCOVERY \
		-A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		-stand_call_conf 30 \
		-stand_emit_conf 30 \
		--sample_ploidy $ploidy \
		-nt 1 -nct $thr  \
		--dbsnp $refdir/$dbsnp $intervals \
		-o $outfile
		
		exitcode=$?

		echo `date`		
		if [ $exitcode -ne 0 ]
		then
			MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
			exit $exitcode;
		fi
		if [ ! -s $outfile ] 
		then
			MSG="$outfile HaplotypeCaller file not created. vcall failed."
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
			exit 1;
		fi           
       fi   ## End calling the autosomes

       set +x; echo -e "\n\n" >&2; 
       echo "############################################################################################" >&2
       echo "############################################################################################" >&2
       echo "####################     End HaplotypeCaller   with GENDERINFO file     ####################" >&2
       echo "############################################################################################" >&2
       echo -e "\n\n" >&2; set -x; 

else
       set +x; echo -e "\n\n" >&2; 
       echo "############################################################################################" >&2       
       echo "############################################################################################" >&2
       echo "##########   HaplotypeCaller   without GENDERINFO file                           ###########" >&2
       echo "############################################################################################" >&2        
       echo -e "\n\n" >&2; set -x;
       
	if [ $chr == "MT" ]
	then
		ploidy=1
	else
		ploidy=2
	fi
   
	outfile=$outputfile

	java -Xmx8g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
         -T HaplotypeCaller \
         -R $refdir/$ref \
         -I ${inputfile} \
         --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
         -gt_mode DISCOVERY \
         -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
         -stand_call_conf 30 \
         -stand_emit_conf 30 \
         -nt 1 -nct $thr  \
         --dbsnp $refdir/$dbsnp $intervals \
	--sample_ploidy $ploidy \
         -o $outfile

	exitcode=$?

	echo `date`		
	if [ $exitcode -ne 0 ]
	then
		MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi
	if [ ! -s $outfile ] 
	then
		MSG="$outfile HaplotypeCaller file not created. vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		exit 1;
	fi           

       set +x; echo -e "\n\n" >&2; 
       echo "############################################################################################" >&2
       echo "############################################################################################" >&2
       echo "####################     End HaplotypeCaller   without GENDERINFO file     #################" >&2
       echo "############################################################################################" >&2
       echo -e "\n\n" >&2; set -x; 
         

fi

echo `date`

set +x; echo "############################################################################################" >&2
echo "############################################################################################" >&2
echo "#############  DONE. EXITING NOW                                             ###############" >&2; 
echo "############################################################################################" >&2
set -x;


