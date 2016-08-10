#!/bin/bash
#	
#  script to realign by region the aligned-dedup file(s) 
########################################################
redmine=hpcbio-redmine@igb.illinois.edu
set -x
if [ $# != 13 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi					

echo `date`
umask 0027
scriptfile=$0
realrecaldir=$1
outputfile=$2
inputfile=$3
chr=$4       
RGparms=$5
target=$6
realparms=$7
runfile=$8
elog=$9
olog=${10}
email=${11}
qsubfile=${12}
failedlog=${13}
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
RealignOutputLogs=`dirname $elog`
infile=`basename $inputfile`

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

input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )        
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
thr=`expr $threads "-" 1`
chrinfiles=$( echo $chrinfiles | tr ":" " " )
RGparms=$( echo $RGparms | tr "::" ":" | sed "s/:/\tRG/g" | sed "s/ID\=/RGID\=/" )
real2parms=$( echo $realparms | sed "s/:\n//" | tr ":" "\n" | sed "s/^/-known /g" )
realparms=$( echo $realparms | sed "s/:\n//" | tr ":" "\n" | sed "s/^/--known /g" )

if [ ! -s $inputfile ]
then
    MSG="$inputfile input bam file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog    
    exit 1;
fi


if [ ! -d $realrecaldir ]
then
    MSG="$realrecaldir realign directory not found"
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
    echo -e "program=$scriptfile WARNING at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog    
    #exit 1;
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


if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
     input_type="WGS"
fi
if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
     input_type="WES"
fi

if [ $input_type == "WES" -a $target == "NOTARGET" ]
then
    target=$chr
fi


set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP1: extract chromosome from aligned bam  ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $realrecaldir

$samdir/samtools view -bu -@ $thr -h $inputfile $chr > presorted_wrg.${chr}.$infile  
  
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="split by chromosome samtools command failed exitcode=$exitcode  realignment stopped for $infile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi
if [ ! -s presorted_wrg.${chr}.$infile ]
then
    MSG="split by chromosome samtools command failed it produced an empty file exitcode=$exitcode  realignment stopped for $infile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi

# sorting and creating an index file with novosort

$novodir/novosort --index --threads $thr --tmpdir $realrecaldir -m 8g -o sorted_wrg.${chr}.$infile  presorted_wrg.${chr}.$infile 
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="split by chromosome novosort command failed exitcode=$exitcode  realignment stopped for $infile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi
if [ ! -s sorted_wrg.${chr}.$infile ]
then
    MSG="split by chromosome  novosort command failed. it produced an empty file exitcode=$exitcode  realignment stopped for $infile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi
 
# generating the header file

$samdir/samtools view -H  sorted_wrg.${chr}.$infile > sorted_wrg.${chr}.${infile}.header  

echo `date`

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP2: RealignerTargetCreator               ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

$javadir/java -Xmx8g -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
    -R $refdir/$ref \
    -I sorted_wrg.$chr.$infile \
    -T RealignerTargetCreator \
    -nt $thr \
    -L $target \
    -o realign.$chr.$infile.list $realparms


exitcode=$?
echo `date`

if [ $exitcode -ne 0 ]
then
    MSG="GATK  RealignerTargetCreator command failed exitcode=$exitcode  realignment stopped for $infile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi

echo -e "########### checking to see if the target list is empty                   #############"	
if [ -z realign.$chr.$infile.list ]
then
    ## the target list is empty, we need to skip IndelAligner command
    skipRealign="YES"
    MSG="realign.$chr.$infile.list realignertargetcreator created an empty list. Skipping IndelRealigner cmd"
    echo -e "program=$scriptfile WARNING at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
fi


set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP3: IndelRealigner                       ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

if [ $skipRealign != "YES" ]
then
     echo -e "########### the target list is NOT empty. Proceed with IndelRealigner #############"
     $javadir/java -Xmx8g -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
    -R $refdir/$ref \
    -I sorted_wrg.$chr.$infile \
    -T IndelRealigner \
    -L $target \
    -o $outputfile \
    -targetIntervals realign.$chr.$infile.list $real2parms

    exitcode=$?
    echo `date`
    if [ $exitcode -ne 0 ]
    then
	MSG="indelrealigner command failed exitcode=$exitcode  realignmen stopped"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
    fi
else
     echo -e "########### the target list is empty. Skip IndelRealigner #############"
     mv sorted_wrg.$chr.$infile   $outputfile
fi

if [ ! -s $outputfile ]
then
    MSG="realigned.bam indelrealigner file not created. realignment stopped for $infile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi


