#!/bin/sh
######################################
#  script to calculate quality information of fastq file
######################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=lmainzer@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 8 ];
then
	MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        fastqcdir=$1
        outputdir=$2
        fastqcparms=$3
        fastqfile=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"        
       
        if [ ! -d $fastqcdir ]
        then
              MSG="fastqcdir not found.  exitcode=$exitcode. Execution of the pipeline continues."
              #echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
              exit 0;
        fi
        if [ ! -d $outputdir ]
        then
              MSG="directory for depositing results of fastqc was not found.  exitcode=$exitcode. Execution of the pipeline continues."
              #echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
              exit 0;
        fi
        if [ ! -s $fastqfile ]
        then
              MSG="input reads file for fastqc command was not found.  exitcode=$exitcode. Execution of the pipeline continues."
              #echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
              exit 0;
        fi



        parameters=$( echo $fastqcparms | tr "_" " " )
        cd $outputdir
        $fastqcdir/fastqc -o $outputdir $parameters $fastqfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
              MSG="fastqc command failed.  exitcode=$exitcode. Execution of the pipeline continues."
              #echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
              exit 0;
        fi
        totlines=`ls -1 *.zip | wc -l | cut -d ' ' -f 1`
        if [ $totlines -lt 1 ]
        then
              MSG="fastqc file for $fastqfile was not created. Execution of the pipeline continues."
              #echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
              exit 0;
        fi
fi
