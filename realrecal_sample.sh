#!/bin/bash
#	
#  script to realign and recalibrate a sample without splitting by region
########################################################
redmine=hpcbio-redmine@igb.illinois.edu
set -x
#if [ $# != 12 ];
#then
#	MSG="parameter mismatch."
#        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
#        exit 1;
#fi					

echo `date`
umask 0027
scriptfile=$0
realrecaldir=$1
inputfile=$2
realignOutputfile=$3
recalibratedOutputfile=$4
target=$5
runfile=$6
failedlog=$7
email=$8
qsubfile=$9

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$failedlog\noutputlog=$failedlog"




if [ ! -s $runfile ]
then
	MSG="$runfile configuration file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
indelfile=$( cat $runfile | grep -w INDELFILE | cut -d '=' -f2 )
realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
thr=$threads
runverify=$( cat $runfile | grep -w RUNVERIFYBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
verifydir=$( cat $runfile | grep -w VERIFYBAMDIR | cut -d '=' -f2 )    
freemix_cutoff=$( cat $runfile | grep -w FREEMIX_CUTOFF | cut -d '=' -f2 )
sites=$( cat $runfile | grep -w OMNISITES | cut -d '=' -f2 ) 
deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
RGparms=$( echo $RGparms | tr "::" ":" | sed "s/:/\tRG/g" | sed "s/ID\=/RGID\=/" )
RealignOutputLogs=`dirname $failedlog`
qc_result=$rootdir/QC_Results.txt




if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
	input_type="WES"
fi

if [ ! -d $picardir ]
then
	MSG="$picardir picard directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -d $samdir ]
then
	MSG="$samdir samtools directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -d $sambambadir ]
then
	MSG="$sambambadir sambamba directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -d $gatk ]
then
	MSG="$gatk GATK directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi


if [ -z $javadir ]
then
	MSG="Value for JAVADIR must be specified in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
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

if [ ! -s $refdir/$dbSNP ]
then
	MSG="$dbSNP DBSNP file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -s $indelfile ]
then
	MSG="$$indelfile indels file  not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi


if [ ! -s $inputfile ]
then
	MSG="$inputfile aligned bam file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi


if [ ! -d $realrecaldir ]
then
	MSG="$realrecaldir realign directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi


if [ $runverify != "1" -a $runverify != "0" -a $runverify != "YES" -a $runverify != "NO" ]
then
	MSG="Invalid value for RUNVERIFYBAM=$runverify"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ $runverify == "1" ]
then
	runverify="YES"
fi

if [ $runverify == "0" ]
then
	runverify="NO"
fi   

if [ `expr ${#deliveryfolder}` -lt 2 ]
then
	delivery=$rootdir/delivery
else
	delivery=$rootdir/$deliveryfolder
fi

if [ ! -d $delivery/Cleaned_BAMS ]
then
	mkdir -p $delivery/Cleaned_BAMS
fi


set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  PREPARATORY WORK                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $realrecaldir/..
sample=`basename $PWD`
cd $realrecaldir

inputFilename=$sample
realtargetfile=${inputFilename}.list             # temp file created by TargetCreator command
recalReport=${inputFilename}.report.grp          # temp file created by BQRS command
precalmdBam=${inputFilename}.precalmd.bam        # temp file created by PrintReads command
recalCsv=${inputFilename}.data.csv               # temp file created by CountCov command
verifiedOut=${inputFilename}.verifyOutput        # output prefix for verifybam command
skipRealign="NO"                                 # flag to skip realignment if  $realtargetfile is empty
intervals=" "                                    # for file with WES intervals list, if it is provided

if [ $input_type == "WES" -a $target != "NOTARGET" -a -s $target ]
then
	intervals=" -L $target "
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP1: realignment                          ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
    -R $refdir/$ref \
    -I $inputfile \
    -T RealignerTargetCreator  $intervals \
    -nt $thr \
    -known $indelfile \
    -o $realtargetfile


exitcode=$?
echo `date`

echo -e "########### checking to see if command failed                   #############"
if [ $exitcode -ne 0 ]
then
	MSG="RealignerTargetCreator command failed exitcode=$exitcode realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

echo -e "########### checking to see if the target list is empty                   #############"	
if [ -z $realtargetfile ]
then
	## the target list is empty, we need to skip IndelAligner command
	skipRealign="YES"
fi

if [ $skipRealign != "YES" ]
then
	     echo -e "########### the target list is NOT empty. Proceed with IndelRealigner #############"
     java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
    -R $refdir/$ref \
    -I $inputfile \
    -T IndelRealigner \
    -o $realignOutputfile \
    -targetIntervals $realtargetfile

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
		MSG="indelrealigner command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
else
	echo -e "########### the target list is empty. Skip IndelRealigner #############"
	ln  $inputfile  $realignOutputfile
fi

if [ ! -s $realignOutputfile ]
then
	MSG="$realignOutputfile indelrealigner file not created. realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP2: recalibration                        ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


if [ "$recalibrator" == "BQSR" ]
then
	set +x; echo -e "\n\n" >&2;       
	echo "#################################################################################" >&2
	echo "                            recalibrator == BQSR" >&2
	echo "#################################################################################" >&2
	echo -e "\n\n" >&2; set -x;         

	java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	--knownSites $refdir/$dbSNP \
	--knownSites $indelfile $intervals \
	-I $realignOutputfile \
	-T BaseRecalibrator \
	--out $recalReport \
	-nct $thr

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
		MSG="BQSR command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi


	if [ ! -s $recalReport ]
	then
		MSG="$recalReport recalibration report file not created, realignment-recalibration stopped for for $infile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
	
	echo `date`


	java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	-I $realignOutputfile \
	-T PrintReads \
	-BQSR $recalReport \
	--out $precalmdBam \
	-nct $thr

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
		MSG="BQSR PrintReads command failed exitcode=$exitcode  recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi


	if [ ! -s $precalmdBam ]
	then
		MSG="$precalmdBam recalibrated file not created, realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
else
	set +x; echo -e "\n\n" >&2;       
	echo "#################################################################################" >&2
	echo "                            recalibrator IS NOT BQSR" >&2
	echo "#################################################################################" >&2
	echo -e "\n\n" >&2; set -x;         


	java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	--knownSites $refdir/$dbSNP \
	--knownSites $indelfile $intervals \
	-I $realignOutputfile \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-recalFile  $recalCsv

	exitcode=$?
	echo `date`		
	if [ $exitcode -ne 0 ]
	then
		MSG="countcovariates command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi

	if [ ! -s $recalCsv ]
	then
		MSG="$recalCsv countcovariates file not created realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
	echo `date`

	java -Xmx8g -Xms1024m  -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	-I $realignOutputfile \
	-T TableRecalibration \
	--out $precalmdBam \
	-recalFile $recalCsv 

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
		MSG="tablerecalibration command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi
	echo `date`		


	if [ ! -s $precalmdBam ]
	then
		MSG="$precalmdBam tablerecalibration file not created realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi

fi # end choosing between BQSR and CountCovariates/TableRecalibration

set +x; echo -e "\n\n" >&2;                
echo "#################################################################################" >&2
echo "#########  STEP3: run samtools calmd to Generate the MD tag            ##########" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

$samdir/samtools calmd -Erb $precalmdBam $refdir/$ref > $recalibratedOutputfile
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="samtools calmd command failed with exitcode=$exitcode. realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi
if [ ! -s $recalibratedOutputfile ]
then
	MSG="$recalibratedOutputfile output file not created. realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n" >&2;                
echo "#################################################################################" >&2
echo "#########  STEP4: create index file to make GATK happy                 ##########" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;
 

$samdir/samtools index $recalibratedOutputfile
exitcode=$?
echo `date`		
if [ $exitcode -ne 0 ]
then
	MSG="samtools index command failed with  exitcode=$exitcode realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

echo `date`    


set +x; echo -e "\n\n" >&2;                
echo "#################################################################################" >&2
echo "#########  STEP5: copy outputfile to delivery folder                   ##########" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cp $recalibratedOutputfile $delivery/Cleaned_BAMS

exitcode=$?
echo `date`		
if [ $exitcode -ne 0 ]
then
	MSG="copy of results to delivery folder  failed with  exitcode=$exitcode realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

echo `date`
    
set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP6: QC with VerifyBam - optional         ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


if [ $runverify == "YES" -a ! -s $sites ]
then
	MSG="$sites file not found cannot perform QC exitcode=$exitcode abort QC for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 0;	
fi


if [ $runverify == "YES" ]
then	

	$verifydir/verifyBamID --vcf $sites --bam $recalibratedOutputfile --out $verifiedOut --verbose --ignoreRB --chip-none 
	exitcode=$?
	`echo date`

	if [ $exitcode -ne 0 ]
	then
		MSG="verifyBamID command failed exitcode=$exitcode abort QC for $inputfile"
		exit 0;
	fi

	if [ ! -s $verifiedOut* ]
	then
		MSG="verifyBamID did not produced output files  exitcode=$exitcode abort QC for $inputfile"
		exit 0;
	fi

	set +x; echo -e "\n\n" >&2; 
	echo -e "##########################################################################" >&2
	echo -e "#############  parse output file                         #################" >&2
	echo -e "##########################################################################" >&2
	echo -e "\n\n" >&2; set -x;


	verifyfile=`find ./ -name "*.selfSM"`
	filterflag=""
	freemix=$( cat $verifyfile | sed '2q;d' | cut -f 7 ) 

	set +x; echo -e "\n\n" >&2; 
	echo -e "##########################################################################" >&2
	echo -e "checking that we actually grabbed something that is a number" >&2
	echo -e "using a weird trick to do the comparison in bash between noninteger variables" >&2
	echo -e "##########################################################################" >&2
	echo -e "\n\n" >&2; set -x;    

	if [ $(echo "$freemix < $freemix_cutoff"|bc) -eq 1 ]
	then
		set +x; echo -e "\n\n" >&2;    
		echo -e "##########################################################################" >&2
		echo -e "sample-$sample PASSED verifyBamID filter" >&2
		echo -e "##########################################################################" >&2
		echo -e "\n\n" >&2; set -x;	

		filterflag="PASSED"	
	else

		set +x; echo -e "\n\n" >&2;  
		echo -e "##########################################################################" >&2
		echo -e "sample-$sample FAILED verifyBamID filter" >&2
		echo -e "##########################################################################" >&2
		echo -e "\n\n" >&2; set -x;	

		filterflag="FAILED"
	fi
	
	# writing result to QC_Results
	detail=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
	echo -e "$detail" >> $qc_result
	


fi # end if runverify

echo `date`

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  DONE. Exiting now                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


