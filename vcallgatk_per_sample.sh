#!/bin/bash 
#	
#  script to perform variant calling with HaplotypeCaller ONLY on an entire sample
#  This module is called from within the realign module
######################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu

set -x
echo `date`
ulimit -s unlimited
umask 0027
scriptfile=$0
realignedbam=$1
gvcfFile=$2
plainVcfFile=$3
outputdir=$4
target=$5
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
deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
snvcaller=$( cat $runfile | grep -w SNV_CALLER | cut -d '=' -f2 )
skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
dbsnp=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
indelfile=$( cat $runfile | grep -w INDELFILE | cut -d '=' -f2 )
genderinfo=$( cat $runfile | grep -w GENDERINFORMATION | cut -d '=' -f2 )
variantcmd=$( cat $runfile | grep -w VARIANT_CMD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
variantAnalysis=$( cat $runfile | grep -w VARIANT_ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
gvcf2vcf=$( cat $runfile | grep -w CONVERT_GVCF2VCF | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
tabixdir=$( cat $runfile | grep -w TABIXDIR | cut -d '=' -f2 )
qc_result=$rootdir/QC_Results.txt

if [ $variantAnalysis == "GENOTYPE" -o $variantAnalysis == "GENOTYPING" ]
then
	MSG="VARIANT_ANALYSIS=$variantAnalysis This case is not analyzed in this script; it needs to be done by region rather than by sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi

if [ $variantcmd != "HC" -a $variantcmd != "HAPLOTYPECALLER"  ]
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

if [ $skipvcall == "1" -o $skipvcall == "YES" ]
then
	echo "skipping the execution of this variant calling module"
	exit 0;
fi
if [ -z $snvcaller ]
then
	MSG="$snvcaller snvcaller tool was not specified in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi
if [ $variantcmd == "HC" -o $variantcmd == "HAPLOTYPECALLER"  ]
then
    variantcmd="HC"
fi

if [ $snvcaller == "GATK" -a $variantcmd != "HC" ]
then
	MSG="VARIANT_CMD=$variantcmd specified in configuration file. This script only runs HaplotypeCaller"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi
if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
    input_type="WGS"
fi

if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
    input_type="WES"
fi


if [ ! -s $realignedbam ]
then
	MSG="$realignedbam realigned bam file not found"
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

set +x; echo -e "\n\n########### checking GVCF 2 VCF conversion     #############\n\n" >&2; set -x;

if [ `expr ${#gvcf2vcf}` -lt 1 ]
then
	gvcf2vcf="NO"
fi
if [ $gvcf2vcf == "1" ]
then
	gvcf2vcf="YES"
fi
if [ $gvcf2vcf == "0" ]
then
	gvcf2vcf="NO"
fi 

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  PREPARATORY WORK                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $outputdir/..
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
         -I $realignedbam \
         --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
         -gt_mode DISCOVERY \
         -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
         -stand_call_conf 30 \
         -stand_emit_conf 30 \
         -nt 1 -nct $thr \
         --dbsnp $refdir/$dbsnp $intervals \
         -o $gvcfFile

exitcode=$?

echo `date`		
if [ $exitcode -ne 0 ]
then
	MSG="HaplotypeCaller command failed  exitcode=$exitcode varcalling stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

if [ ! -s $gvcfFile ] 
then
	MSG="$ogvcfFile HaplotypeCaller file not created. varcalling stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	exit 1;
fi 

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP3: COPY OUTPUT FILES TO DELIVERY FOLDER      ##############" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


echo `date`

cp $gvcfFile $deliverydir

echo `date`	 

set +x; echo -e "\n\n" >&2;
echo -e "##################################################################################"  >&2
echo -e "########### STEP4: convert GVCF to PLAIN    VCF                              #####" >&2
echo -e "##################################################################################" >&2
echo -e "\n\n" >&2; set -x;

if [ $gvcf2vcf == "NO" ]
then

	set +x; echo -e "\n\n" >&2;
	echo -e "##################################################################################"  >&2
	echo -e "########### DONE. Exiting now                                                #####" >&2
	echo -e "##################################################################################" >&2
	echo -e "\n\n" >&2; set -x;

	exit 0;
fi

set +x; echo -e "\n\n########### variables      #############\n\n" >&2; set -x;

plainTmpVcf=tmp_${sample}.regular.raw.vcf


set +x; echo -e "\n\n########### run GenotypeGVCFs to convert GVCF to VCF     #############\n\n" >&2; set -x;

java -Xmx50g  -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	 -R $refdir/$ref \
	 --dbsnp $refdir/$dbsnp \
         -T GenotypeGVCFs \
         -o $plainTmpVcf \
         -nt $thr \
         --variant $gvcfFile



exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="GATK GenotypeGVCFs  command failed exitcode=$exitcode. sample $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

if [ ! -s $plainTmpVcf ]
then
	MSG="GATK GenotypeGVCFs did not generate file exitcode=$exitcode. sample $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi 

set +x; echo -e "\n\n########### run VariantAnnotator to add missing tags that UnifiedGenotyper always includes  #############\n\n" >&2; set -x;


java -Xmx50g  -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	 -R $refdir/$ref \
	 --dbsnp $refdir/$dbsnp \
	 -T VariantAnnotator \
	 -I $realignedbam \
         -V $plainTmpVcf \
         --disable_auto_index_creation_and_locking_when_reading_rods \
         -A VariantType \
         -A AlleleBalance -A BaseCounts -A BaseQualityRankSumTest -A ChromosomeCounts \
         -A Coverage -A FisherStrand -A GCContent -A HaplotypeScore \
         -A HomopolymerRun -A InbreedingCoeff -A LowMQ -A MappingQualityRankSumTest \
         -A MappingQualityZero -A NBaseCount -A QualByDepth -A RMSMappingQuality \
         -A ReadPosRankSumTest -A SpanningDeletions -A TandemRepeatAnnotator \
         -nt $thr \
         -o $plainVcfFile 

exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="GATK VariantAnnotator  command failed exitcode=$exitcode. sample $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

if [ ! -s $plainVcfFile ]
then
	MSG="GATK VariantAnnotator did not generate file exitcode=$exitcode. sample $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi 

set +x; echo -e "\n\n########### copy file to delivery folder  #############\n\n" >&2; set -x;

echo `date`

cp $plainVcfFile $deliverydir

echo `date`


set +x; echo -e "\n\n" >&2;
echo -e "##################################################################################"  >&2
echo -e "########### DONE. Exiting now                                                #####" >&2
echo -e "##################################################################################" >&2
echo -e "\n\n" >&2; set -x;
