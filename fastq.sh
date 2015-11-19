#!/bin/bash
######################################
#  script to calculate quality information of fastq file
######################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=lmainzer@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 9 ]
then
	MSG="parameter mismatch"
        echo -e "Program $0 stopped. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
        exit 1;
else					

        echo -e "\n\n############# CHECK FASTQ QUALITY  ###############\n\n" >&2

	set -x
	echo `date`
        umask 0027
        scriptfile=$0
        runfile=$1
        fastqcdir=$2
        outputdir=$3
        fastqcparms=$4
        inputfiles=$5
        elog=$6
        olog=$7
        email=$8
        qsubfile=$9
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"        
       
        #set +x; echo -e "\n\n############# preliminary parameter check  ###############\n\n" >&2
 
        # get folder path to logs
        LogsPath=`dirname $elog`

        if [ ! -s $runfile ]
        then
              MSG="$runfile file not found.  This script stops, but execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi
        if [ ! -d $fastqcdir ]
        then
              MSG="fastqcdir not found.  This script stops, but execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi
        if [ ! -d $outputdir ]
        then
              MSG="directory for depositing results of fastqc was not found.  Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi

        echo -e "preparing delivery folder"
        deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
        rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        if [ `expr ${#deliveryfolder}` -lt 2 ]
        then
             delivery=$rootdir/delivery
        else
             delivery=$rootdir/$deliveryfolder
        fi
        if [ ! -d $delivery/fastqc ]
        then
            mkdir -p $delivery/fastqc
        fi

        echo -e "gathering the pieces of the fastqc command"

        fastqfile=$( echo $inputfiles | tr ":" " " )
        parameters=$( echo $fastqcparms | tr "_" " " )

        echo -e "running the fastqc command"

        cd $outputdir
        $fastqcdir/fastqc -o $outputdir $parameters $fastqfile
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
              MSG="fastqc command failed.  exitcode=$exitcode. Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi

        echo -e "checking that output file was created"

        totlines=`ls -1 *.zip | wc -l | cut -d ' ' -f 1`
        if [ $totlines -lt 1 ]
        then
              MSG="fastqc file for $fastqfile was not created. Execution of the pipeline continues."
              echo -e "program=$scriptfile failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> ${LogsPath}/FAILED_fastq
              exit 1;
        fi

        echo -e "transfer results to delivery folder"

        cp *.zip $delivery/fastqc/
	echo `date`

fi
