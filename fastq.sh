#!/bin/bash
######################################
#  script to calculate quality information of fastq file
######################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=lmainzer@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 8 ]
then
	MSG="parameter mismatch"
        echo -e "Program $0 stopped. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
        exit 1;
else					

        echo -e "\n\n############# CHECK FASTQ QUALITY  ###############\n\n" >&2

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
       
        set +x; echo -e "\n\n############# preliminary parameter check  ###############\n\n" >&2
 
        # get folder path to logs
        LogsPath=`dirname $elog`

        if [ ! -d $fastqcdir ]
        then
              MSG="fastqcdir not found.  Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi
        if [ ! -d $outputdir ]
        then
              MSG="directory for depositing results of fastqc was not found.  Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi
        if [ ! -s $fastqfile ]
        then
              MSG="input reads file for fastqc command was not found.  Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi



        parameters=$( echo $fastqcparms | tr "_" " " )
        cd $outputdir
        $fastqcdir/fastqc -o $outputdir $parameters $fastqfile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
              MSG="fastqc command failed.  exitcode=$exitcode. Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi
        totlines=`ls -1 *.zip | wc -l | cut -d ' ' -f 1`
        if [ $totlines -lt 1 ]
        then
              MSG="fastqc file for $fastqfile was not created. Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi
fi
