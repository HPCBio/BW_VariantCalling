#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 14 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        aligndir=$1
        parms=$2
        ref=$3
        outputdir=$4
	samfile=$5
	bamfile=$6
        Rone=$7
        Rtwo=$8
        scriptdir=$9
        samdir=${10}
        elog=${11}
        olog=${12}
        email=${13}
        qsubfile=${14}
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        parameters=$( echo $parms | tr "_" " " )


        cd $outputdir
        $aligndir/bwa mem $parameters $ref $Rone $Rtwo > $samfile

        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa aln command failed on $R.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $samfile ]
        then
            MSG="$samfile aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        echo `date`

        ## sam2bam conversion
        $samdir/samtools view -bS -o $bamfile $samfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="samtools view command failed.  exitcode=$exitcode. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" |  myproxy-logon -s tfca.ncsa.illinois.edu -p 7512 "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $bamfile ]
        then
            MSG="$bamfile bam file not created. sam2bam step failed during alignment."
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" |  myproxy-logon -s tfca.ncsa.illinois.edu -p 7512 "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        echo `date`

fi
