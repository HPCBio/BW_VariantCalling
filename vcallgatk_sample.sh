#!/bin/bash 
#	
#  script to perform variant calling with HaplotypeCaller ONLY
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
runfile=$5
failedlog=$6
email=$7
qsubfile=$8

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
deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
ped=$( cat $runfile | grep -w PEDIGREE | cut -d '=' -f2 )
allsites=$( cat $runfile | grep -w EMIT_ALL_SITES | cut -d '=' -f2 )
snvcaller=$( cat $runfile | grep -w SNV_CALLER | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
dbsnp=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
indelfile=$( cat $runfile | grep -w INDELFILE | cut -d '=' -f2 )
genderinfo=$( cat $runfile | grep -w GENDERINFORMATION | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
variantcmd=$( cat $runfile | grep -w VARIANT_CMD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
variantAnalysis=$( cat $runfile | grep -w VARIANT_ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
tabixdir=$( cat $runfile | grep -w TABIXDIR | cut -d '=' -f2 )
qc_result=$rootdir/QC_Results.txt

if [ $variantAnalysis == "GENOTYPE" -o $variantAnalysis == "GENOTYPING" ]
then
	MSG="VARIANT_ANALYSIS=$variantAnalysis This case is not analyzed in this script; it needs to be done by region rather than by sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
else

if [ $variantcmd != "HC" -o $variantcmd != "HAPLOTYPECALLER"  ]
then
	MSG="VARIANT_CMD=$variantcmd this case is not analyzed in this script"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

if [ `expr ${#deliveryfolder}` -lt 2 ]
then
    deliverydir=$rootdir/delivery/Vcfs
else
    deliverydir=$rootdir/$deliveryfolder/Vcfs
fi

if [ ! -d $deliverydir ]
then
    `mkdir -p $deliverydir`
fi

if [ $skipvcall != "1" -a $skipvcall != "0" -a $skipvcall != "YES" -a $skipvcall != "NO" ]
then
	MSG="Invalid value for SKIPVCALL=$skipvcall"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
elif [ $skipvcall == "1" -o $skipvcall == "YES" ]
then
	echo "skipping the execution of this variant calling module"
	exit 0;
fi

if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
    input_type="WGS"
fi
if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
    input_type="WES"
fi

if [ $variantcmd == "HC" -o $variantcmd == "HAPLOTYPECALLER"  ]
then
    variantcmd="HAPLOTYPECALLER"
fi

if [ ! -s $inputfile ]
then
	MSG="$inputfile realigned bam file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

if [ ! -d $outputdir ]
then
    mkdir -p $outputdir
fi

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
if [ -z $snvcaller ]
then
	MSG="$snvcaller snvcaller tool was not specified in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

if [ $snvcaller == "GATK" -a $variantcmd != "HAPLOTYPECALLER" ]
then
	MSG="VARIANT_CMD=$variantcmd specified in configuration file. This script only runs HaplotypeCaller"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

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
echo `date`

#tabix needs to be on your path for GATK to produce *tbi files
export PATH=${PATH}:$tabixdir

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  PREPARATORY WORK                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $$outputdir/..
sample=`basename $PWD`
cd $$outputdir
intervals=" "                                    # for file with WES intervals list, if it is provided
qctestflag=""                                    # for qc test

if [ $input_type == "WES" -a $target != "NOTARGET" -a -s $target ]
then
	intervals=" -L $target "
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP1: run HaplotypeCaller                  ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

java -Xmx8g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
         -T HaplotypeCaller \
         -R $refdir/$ref \
         -I $inputfile \
         --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
         -gt_mode DISCOVERY \
         -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
         -stand_call_conf 30 \
         -stand_emit_conf 30 \
         -nt 1 -nct $thr \
         --dbsnp $refdir/$dbsnp $intervals \
         -o $outputfile

exitcode=$?
echo `date`		
if [ $exitcode -ne 0 ]
then
	MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

if [ ! -s $outputfile ] 
then
	MSG="$outputfile HaplotypeCaller file not created. varcalling stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi 

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP2: simple QC step                       ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

# total variants in outputfile is approx. equal to the number of NON-HEADER lines in that file

numvariants=$( grep -v "^#" -c $outputfile )

if [ `expr ${#numvariants}` -lt 1 ]
then
	# when the output file has only header lines, the command to calculate numvariants returns nothing
	qctestflag="FAILED"
elfi [ ${numvariants} -lt 1 ]
	# when the output file has only header lines, the command to calculate numvariants returns 0
	qctestflag="FAILED"
else
	# otherwise, the command to calculate numvariants returns positive integer
	qctestflag="PASSED"
fi

# writing result to QC_Results
detail=$( echo -e "$sample\t$fqctestflag\tTotalVariants=$numvariants\tMinVariant_Count=1" )
echo -e "$detail" >> $qc_result


set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP3: COPY OUTPUT FILES TO DELIVERY FOLDER      ##############" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


echo `date`

cp ${outputfile $deliverydir
cp ${outputfile}.tbi $deliverydir

echo `date`	 

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  DONE. EXITING NOW                                ##############" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;
