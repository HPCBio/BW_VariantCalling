#!/bin/bash
#	
#  script to realign and recalibrate the aligned file(s) by region ONLY.
#  Use mergeCleanBams.sh  to perform QC with VerifyBamID and to copy results to delivery folder
########################################################
#redmine=hpcbio-redmine@igb.illinois.edu
remine=grendon@illinois.edu

set -x

#if [ $# != 11 ];
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
chr=$5  
target=$6
indelparam=$7
runfile=$8
failedlog=$9
email=${10}
qsubfile=${11}

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"



set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )


set +x; echo -e "\n\n########### checking tool directories     #############\n\n" >&2; set -x;

if [ ! -d $samdir ]
then
	MSG="$samdir samtools directory not found"
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

if [ ! -s $refdir/$dbSNP ]
then
	MSG="$dbSNP DBSNP file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi


set +x; echo -e "\n\n########### checking inputfile     #############\n\n" >&2; set -x;

if [ ! -s $inputfile ]
then
	MSG="$inputfile aligned bam file not found"
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


set +x; echo -e "\n\n########### checking input type     #############\n\n" >&2; set -x;

if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
	input_type="WES"
fi


set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  PREPARATORY WORK                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

set +x; echo -e "\n\n########### define sample     #############\n\n" >&2; set -x;

cd $realrecaldir/..
sample=`basename $PWD`
cd $realrecaldir

set +x; echo -e "\n\n########### define temp files and variables     #############\n\n" >&2; set -x;

inputFilename=${sample}.${chr}
alignedbam=${inputFilename}.alignedDedup.bam
alignedSortedbam=${inputFilename}.alignedDedup.sorted.bam
realtargetfile=${inputFilename}.list             # temp file created by TargetCreator command
recalReport=${inputFilename}.report.grp          # temp file created by BQRS command
precalmdBam=${inputFilename}.precalmd.bam        # temp file created by PrintReads command
recalCsv=${inputFilename}.data.csv               # temp file created by CountCov command
skipRealign="NO"                                 # flag to skip realignment if  $realtargetfile is empty
intervals=" "                                    # for file with WES intervals list, if it is provided

set +x; echo -e "\n\n########### define parameters for GATK commands     #############\n\n" >&2; set -x;

if [ $input_type == "WES" -a $target != "NOTARGET" -a -s $target ]
then
	intervals=" -L $target "
else
	intervals=" -L $chr "
fi

indelfile=$( echo $indelparams | tr ":" " " )

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP1: extract chromosome from aligned bam  ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

$samdir/samtools view -bu -@ $thr -h $inputfile $chr > $alignedbam  
 
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="split by chromosome samtools command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi
if [ ! -s $alignedbam ]
then
	MSG="split by chromosome samtools command failed it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi


set +x; echo -e "\n\n###### sort the new aligned.bam file     ##############\n\n" >&2; set -x;

$novodir/novosort --index --threads $thr --tmpdir $realrecaldir -m 16g -o  $alignedSortedbam $alignedbam
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="novosort command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi
if [ ! -s $alignedSortedbam ]
then
	MSG="novosort command failed. it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi


set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP2: realignment                          ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
    -R $refdir/$ref \
    -I $alignedSortedbam \
    -T RealignerTargetCreator \
    -nt $thr \
    -known $indelfile  $intervals \
    -o $realtargetfile


exitcode=$?
echo `date`

set +x; echo -e "\n\n########### checking to see if RealignerTargetCreator command failed     #############\n\n" >&2; set -x;

if [ $exitcode -ne 0 ]
then
	MSG="RealignerTargetCreator command failed exitcode=$exitcode realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n########### checking to see if the target list is empty    #############\n\n" >&2; set -x;
	
if [ -z $realtargetfile ]
then
	## the target list is empty, we need to skip IndelAligner command
	skipRealign="YES"
fi

if [ $skipRealign != "YES" ]
then
	set +x; echo -e "\n\n########### the target list is NOT empty. Proceed with IndelRealigner #############\n\n" >&2; set -x;
	     
     java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
    -R $refdir/$ref \
    -I $alignedSortedbam \
    -T IndelRealigner \
    -o $realignOutputfile \
    -targetIntervals $realtargetfile

	exitcode=$?
	echo `date`
	set +x; echo -e "\n\n########### check to see if IndelRealigner command failed #############\n\n" >&2; set -x;
	
	if [ $exitcode -ne 0 ]
	then
		MSG="indelrealigner command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
else
	set +x; echo -e "\n\n########### the target list is empty. Skip IndelRealigner #############\n\n" >&2; set -x;
	ln  $alignedSortedbam  $realignOutputfile
fi

set +x; echo -e "\n\n########### check to see if IndelRealigner created an output file #############\n\n" >&2; set -x;

if [ ! -s $realignOutputfile ]
then
	MSG="$realignOutputfile indelrealigner file not created. realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP3: recalibration. It could be BQSR or CountCovariates #####" >&2
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
	
	set +x; echo -e "\n\n########### check to see if BaseRecalibrator command failed #############\n\n" >&2; set -x;
	
	if [ $exitcode -ne 0 ]
	then
		MSG="BQSR command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi

	set +x; echo -e "\n\n########### check to see if BaseRecalibrator did not create a file #############\n\n" >&2; set -x;
	
	if [ ! -s $recalReport ]
	then
		MSG="$recalReport recalibration report file not created, realignment-recalibration stopped for for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
	
	echo `date`

	set +x; echo -e "\n\n########### Run PrintReads command #############\n\n" >&2; set -x;
	
	java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	-I $realignOutputfile \
	-T PrintReads \
	-BQSR $recalReport \
	--out $precalmdBam \
	-nct $thr

	exitcode=$?
	echo `date`
	
	set +x; echo -e "\n\n########### check to see if PrintReads command failed #############\n\n" >&2; set -x;
	
	if [ $exitcode -ne 0 ]
	then
		MSG="BQSR PrintReads command failed exitcode=$exitcode  recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi

	set +x; echo -e "\n\n########### check to see if PrintReads produced a file #############\n\n" >&2; set -x;
	
	if [ ! -s $precalmdBam ]
	then
		MSG="$precalmdBam recalibrated file not created, realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
else
	set +x; echo -e "\n\n" >&2;       
	echo "#################################################################################" >&2
	echo "                            recalibrator == CountCovariantes + TableRecalibration" >&2
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
	
	set +x; echo -e "\n\n########### check to see if CountCovariates command failed #############\n\n" >&2; set -x;
	
	if [ $exitcode -ne 0 ]
	then
		MSG="countcovariates command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi

	set +x; echo -e "\n\n########### check to see if CountCovariates produced a file #############\n\n" >&2; set -x;
	
	if [ ! -s $recalCsv ]
	then
		MSG="$recalCsv countcovariates file not created realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi
	echo `date`

	set +x; echo -e "\n\n########### run TableRecalibration command #############\n\n" >&2; set -x;
	
	java -Xmx8g -Xms1024m  -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	-I $realignOutputfile \
	-T TableRecalibration \
	--out $precalmdBam \
	-recalFile $recalCsv 

	exitcode=$?
	echo `date`
	
	set +x; echo -e "\n\n########### check to see if TableRecalibration command failed #############\n\n" >&2; set -x;
	
	if [ $exitcode -ne 0 ]
	then
		MSG="tablerecalibration command failed exitcode=$exitcode  realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit $exitcode;
	fi
	echo `date`		

	set +x; echo -e "\n\n########### check to see if TableRecalibration produced a file #############\n\n" >&2; set -x;
	
	if [ ! -s $precalmdBam ]
	then
		MSG="$precalmdBam tablerecalibration file not created realignment-recalibration stopped for $inputfile"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
		exit 1;
	fi

fi # end choosing between BQSR and CountCovariates/TableRecalibration

set +x; echo -e "\n\n" >&2;                
echo "#################################################################################" >&2
echo "#########  STEP4: run samtools calmd to Generate the MD tag            ##########" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

$samdir/samtools calmd -Erb $precalmdBam $refdir/$ref > $recalibratedOutputfile
exitcode=$?
echo `date`

set +x; echo -e "\n\n########### check to see if samtools calmd  command failed #############\n\n" >&2; set -x;

if [ $exitcode -ne 0 ]
then
	MSG="samtools calmd command failed with exitcode=$exitcode. realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi

set +x; echo -e "\n\n########### check to see if samtools calmd  produced a file #############\n\n" >&2; set -x;

if [ ! -s $recalibratedOutputfile ]
then
	MSG="$recalibratedOutputfile output file not created. realignment-recalibration stopped for $inputfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

set +x; echo -e "\n\n" >&2;                
echo "#################################################################################" >&2
echo "#########  STEP5: create index file to make GATK happy                 ##########" >&2
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
echo "################  DONE. Exiting now                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;


