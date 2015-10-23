#!/bin/bash
##redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 14 ]
then
        MSG="parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
umask 0037	
set -x
echo `date`
scriptfile=$0
aligndir=$1
parms=$2
ref=$3
outputdir=$4
bamprefix=$5
Rone=$6
Rtwo=$7
runfile=$8
elog=${9}
olog=${10}
email=${11}
qsubfile=${12}
RGparms=${13}
AlignOutputLogs=${14}

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"



if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

echo "####################################################################################################"
echo "#####################################                       ########################################"
echo "##################################### PARSING RUN INFO FILE ########################################"
echo "##################################### AND SANITY CHECK      ########################################"
echo "####################################################################################################"

rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d "=" -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
markduptool=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d "=" -f2 | tr '[a-z]' '[A-Z]' )
samprocessing=$( cat $runfile | grep -w SAMPROCESSING | cut -d "=" -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
samblasterdir=$( cat $runfile | grep -w SAMBLASTERDIR | cut -d "=" -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d "=" -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d "=" -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d "=" -f2 )
alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
header=$( echo $RGparms  | tr ":" "\t" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${header}"  | tr "=" ":" )
dup_cutoff=$( cat $runfile | grep -w DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w MAP_CUTOFF | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )

if [ ! -d $picardir ]
then
    MSG="$picardir picard directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

if [ ! -d $samdir ]
then
    MSG="$samdir samtools directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ ! -d $samblasterdir ]
then
    MSG="$samblasterdir samblaster directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ ! -d $sambambadir ]
then
    MSG="$sambambadir sambamba directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ ! -d $novodir ]
then
    MSG="$novodir novocraft directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ `expr ${#markduptool}` -lt 1 ]
then
    MSG="MARKDUPLICATESTOOL=$markduptool invalid value"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

if [ ! -s $rootdir/SAMPLENAMES_multiplexed.list ]
then
    MSG="$rootdir/SAMPLENAMES_multiplexed.list file not found. This is not a Baylor sample"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

threads=`expr $thr "-" 1`

all_exitcodes=0


if [ $markduptool != "PICARD" ]
then
    MSG="$markduptool IS NOT PICARD. This is not a Baylor sample"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

        echo -e "####################################### CASE2:   MARKDUPLICATESTOOL == PICARD  OF A BAYLOR SAMPLE              #################"
        echo -e "####################################### part 1: align, then convert sam to bam, then sort, then QC test        #################"
        echo -e "####################################### part 2: mark duplicates, then sort                                     #################"
        echo -e "################################################################################################################################"
        echo -e "####################################### This script performs part 1                                            #################"
        
        cd $outputdir
        echo `date`
        #AlignCase="CHPC" #testing alignment Gerrit's way

        $aligndir/bwa mem  $alignparms -R "${rgheader}" $ref $Rone $Rtwo |  $samdir/samtools view -bSu -> ${bamprefix}.tmp.bam
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem command failed on $Rone $Rtwo   exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit $exitcode;
        fi        
        if [ ! -s ${bamprefix}.tmp.bam ]
        then
            MSG="${bamprefix}.tmp.bam aligned bam file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit 1;
        fi        

        ### sometimes we may have a BAM file with NO alignmnets, just the header
#        numAlignments=$( $samdir/samtools view -c ${bamprefix}.tmp.bam ) 
        numAlignments=$( $sambambadir/sambamba view -c -t $thr ${bamprefix}.tmp.bam ) 
        echo `date`
        if [ $numAlignments -eq 0 ]
        then
            MSG="bwa mem command did not produce alignments. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            exit 1;
        fi
        #if [ $AlignCase != "CHPC" ]
        #then
            echo -e "sorting aligned.bam file using novosort"

            $novodir/novosort --index --tmpdir $outputdir --threads $threads -m 16g --kt --compression 1 -o ${bamprefix}.sorted.bam ${bamprefix}.tmp.bam
            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="novosort  failed on ${bamprefix}.tmp.bam  exitcode=$exitcode. alignment failed"
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
                cp $qsubfile $AlignOutputLogs/FAILEDjobs/
                exit $exitcode;
            fi
        
        #else 
        #    echo -e "sorting aligned.bam file  Gerrit's way, using samtools sort"
        #    $samdir/samtools sort -T sort -@ $threads -m 1G -l 1 -o ${bamprefix}.sorted.bam  ${bamprefix}.tmp.bam
        #    exitcode=$?
        #    echo `date`
        #    if [ $exitcode -ne 0 ]
        #    then
        #        MSG="samtools sort  failed on   ${bamprefix}.tmp.bam exitcode=$exitcode. alignment failed"
	#        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
        #        echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
        #        cp $qsubfile $AlignOutputLogs/FAILEDjobs/
        #        exit $exitcode;
        #    fi  
            
        #    $samdir/samtools index ${bamprefix}.sorted.bam
        #    exitcode=$?
        #    echo `date`
        #    if [ $exitcode -ne 0 ]
        #    then
        #        MSG="samtools index  failed on   ${bamprefix}.sorted.bam exitcode=$exitcode. alignment failed"
	#        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
        #        echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
        #        cp $qsubfile $AlignOutputLogs/FAILEDjobs/
        #        exit $exitcode;
        #    fi  
             
        #fi
        
        if [ ! -s ${bamprefix}.sorted.bam ]
        then
            MSG="${bamprefix}.sorted.bam aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit 1;
        fi

        $samdir/samtools view -H ${bamprefix}.sorted.bam > ${bamprefix}.sorted.bam.header


        echo `date`

        echo -e "#######################################################################################################"
        echo -e "########    Done with part 1 of CASE 2                                                        #######"
        echo -e "#######################################################################################################"
