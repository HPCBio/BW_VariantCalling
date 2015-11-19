#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 13 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        umask 0027
        scriptfile=$0
        aligndir=$1
        parms=$2
        ref=$3
        outputdir=$4
	samfile=$5
	bamfile=$6
        R=$7
        scriptdir=$8
        samdir=$9
        elog=${10}
        olog=${11}
        email=${12}
        qsubfile=${13}
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        parameters=$( echo $parms | tr "_" " " )


        cd $outputdir
        $aligndir/bwa mem $parameters $ref $R > $samfile

        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa aln command failed on $R.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $samfile ]
        then
            MSG="$samfile aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
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
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" |  myproxy-logon -s tfca.ncsa.illinois.edu -p 7512 "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $bamfile ]
        then
            MSG="$bamfile bam file not created. sam2bam step failed during alignment."
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" |  myproxy-logon -s tfca.ncsa.illinois.edu -p 7512 "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        echo `date`

fi
