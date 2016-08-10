#!/bin/bash
##redmine=hpcbio-redmine@igb.illinois.edu
#if [ $# != 11 ]
#then
#        MSG="parameter mismatch"
#        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
#        exit 1;
#fi

umask 0027	
set -x
echo `date`
scriptfile=$0
outputdir=$1
inputfiles=$2
alignedbam=$3
dedupbam=$4
RGparms=$5
runfile=$6
elog=$7
olog=$8
email=$9
qsubfile=${10}
failedlog=${11}

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


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
        
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
samprocessing=$( cat $runfile | grep -w SAMPROCESSING | cut -d "=" -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
samblasterdir=$( cat $runfile | grep -w SAMBLASTERDIR | cut -d "=" -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d "=" -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d "=" -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d "=" -f2 )
memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
alignerdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
refindex=$( cat $runfile | grep -w BWAMEMINDEX | cut -d '=' -f2 )
alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 | tr " " "_" )

if [ ! -d $alignerdir ]
then
    MSG="$alignerdir directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $refindexed ]
then
    MSG="$refindexed directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $samdir ]
then
    MSG="$samdir samtools directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $novodir ]
then
    MSG="$novodir novocraft directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi

header=$( echo $RGparms  | tr ":" "\t" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${header}"  | tr "=" ":" )
threads=`expr $thr "-" 1`
fqfiles=$( echo $inputfiles | tr ":" " " )
prefix=${alignedbam%.bam}
samfile=${prefix}.sam
unsortedbam=${prefix}.unsorted.bam

set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    Run bwa mem                                                                         #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $outputdir
echo `date`


$alignerdir/bwa mem  $alignparms -t $threads -R "${rgheader}" $refindex $fqfiles >  $samfile
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="bwa mem command failed on $fqfiles   exitcode=$exitcode. alignment failed"
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi        
if [ ! -s $samfile ]
then
    MSG="$$samfile aligned sam file not created. alignment failed"
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi            

set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    SAM to BAM conversion                                                               #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


$samdir/samtools view -bSu -@ $thr $samfile > $unsortedbam
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="samtools view command failed on $$unsortedbam   exitcode=$exitcode. alignment failed" 
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi        
if [ ! -s $unsortedbam ]
then
    MSG="$$unsortedbam aligned bam file not created. alignment failed"
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi        

set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    QC: count alignments                                                               #######" >&2
echo -e "#######################################################################################################" >&2

echo -e "\n\n" >&2; set -x;
#numAlignments=$( $sambambadir/sambamba view -c -t $thr $unsortedbam ) 
numAlignments=$( $samdir/samtools view -c -@ $thr $unsortedbam )
echo `date`

if [ `expr ${#numAlignments}` -lt 1 ]
then
    MSG="bwa mem command did not produce alignments. alignment failed"
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi

set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    Sort output and generate header                                                     #######" >&2
echo -e "#######################################################################################################" >&2

$novodir/novosort --index --tmpdir $outputdir --threads $threads -m 16g --compression 1 -o $alignedbam $unsortedbam
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="novosort  failed on $unsortedbam  exitcode=$exitcode. alignment failed"
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi


if [ ! -s $alignedbam ]
then
    MSG="$alignedbam aligned-sorted bam file not created. alignment failed"
    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi

$samdir/samtools view -H $alignedbam > ${alignedbam}.header


echo `date`



set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    Done with bwa mem alignment                                                         #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x;