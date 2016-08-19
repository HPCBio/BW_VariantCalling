#!/bin/bash
#	
#  script to merge all clean bams for a sample that were realigned-recalibrated by region
#  it also copies outputfile to delivery folder and runs QC with VerifyBamID
########################################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu

#if [ $# != 12 ];
#then
#	MSG="parameter mismatch."
#        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
#        exit 1;
#fi					


set -x
echo `date`
umask 0027
scriptfile=$0
outputfile=$1
realrecaldir=$2
runfile=$3
failedlog=$4
email=$5
qsubfile=$6

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$failedlog\noutputlog=$failedlog"

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
runverify=$( cat $runfile | grep -w RUNVERIFYBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
verifydir=$( cat $runfile | grep -w VERIFYBAMDIR | cut -d '=' -f2 )    
freemix_cutoff=$( cat $runfile | grep -w FREEMIX_CUTOFF | cut -d '=' -f2 )
sites=$( cat $runfile | grep -w OMNISITES | cut -d '=' -f2 ) 
deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
qc_result=$rootdir/QC_Results.txt
 
set +x; echo -e "\n\n########### checking tool directories     #############\n\n" >&2; set -x;

if [ ! -d $samdir ]
then
	MSG="$samdir samtools directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -d $novodir ]
then
	MSG="$novodir novosort directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -d $verifydir ]
then
	MSG="$verifydir verifydir directory not found"
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

set +x; echo -e "\n\n########### checking output folder     #############\n\n" >&2; set -x;

if [ ! -d $realrecaldir ]
then
	MSG="$realrecaldir realign directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n########### checking option to run QC with VerifyBamID     #############\n\n" >&2; set -x;


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

set +x; echo -e "\n\n########### checking delivery folder     #############\n\n" >&2; set -x;

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
verifiedOut=${inputFilename}.verifyOutput        # output prefix for verifybam command

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP1: gather all files to be merged        ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

`echo date`

chrbamfiles=`find ./ -name "${sample}*.calmd.bam"`

if [ `expr ${#chrbamfiles}` -lt 1 ]
then
	MSG="no bam files to merge for sample: $sample in folder $realrecaldir  mergeCleanedBams failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
	
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP2: merge files with novosort            ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

`echo date`

$novodir/novosort --index --tmpdir $realrecaldir --threads $thr -m 16g -o $outputfile ${sample}*.calmd.bam

exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="novosort  failed exitcode=$exitcode on $sample in folder $realrecaldir  mergeCleanedBams failed"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;

fi

if [ ! -s $outputfile ]
then
    MSG="novosort did not create outputfile file.  $sample in folder $realrecaldir  mergeCleanedBams failed"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi


set +x; echo -e "\n\n" >&2; 
echo -e "##########################################################################" >&2
echo -e "##########################################################################" >&2
echo -e "#############  STEP3: populating delivery folder         #################" >&2
echo -e "##########################################################################" >&2
echo -e "\n\n" >&2; set -x;

`echo date`

cp $outputfile $delivery/Cleaned_BAMS

`echo date`

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP4: QC with VerifyBam - optional         ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


if [ $runverify == "YES" -a ! -s $sites ]
then
	MSG="$sites file not found cannot perform QC $sample in folder $realrecaldir  mergeCleanedBams failed"
	echo -e "program=$scriptfile warning at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 0;	
fi


if [ $runverify == "YES" ]
then	

	$verifydir/verifyBamID --vcf $sites --bam $outputfile --out $verifiedOut --verbose --ignoreRB --chip-none 
	exitcode=$?
	`echo date`

	if [ $exitcode -ne 0 ]
	then
		MSG="verifyBamID command failed exitcode=$exitcode abort QC for $outputfile"
		echo -e "program=$scriptfile warning at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog

		exit 0;
	fi

	if [ ! -s $verifiedOut* ]
	then
		MSG="verifyBamID did not produced output files  exitcode=$exitcode abort QC for $inputfile"
		echo -e "program=$scriptfile warning at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
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


