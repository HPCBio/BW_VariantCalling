#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 4 ]
then
        MSG="Parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        
        exit 1;
else
	set -x
	echo `date`	
        scriptfile=$0
        outputdir=$1
        email=$2
        exitstatus=$3
        reportticket=$4
        numdays=30

	listjobids=$( cat $outputdir/logs/*pbs cat $outputdir/logs/*/*pbs | sort | uniq | tr "\n" "\t" )
	pipeid=$( cat $outputdir/logs/CONFIGUREpbs )

        if [ $exitstatus == "exitok" ]
        then
            MSG="Variant calling workflow with id: [$pipeid] by username: $USER finished successfully at: "$( echo `date` )
        else
            MSG="Variant calling workflow with id: [$pipeid] by username: $USER finished on iforge at: "$( echo `date` )
        fi
        LOGS="Results and execution logs can be found at \n$outputdir\n\nJOBIDS\n\n$listjobids\n\nThis jobid:${PBS_JOBID}\n\n"
        detjobids=""
        nl="\n"
#        for jobid in $listjobids
#        do
#           report=`tracejob -q -n $numdays $jobid`
#           report=$( echo $report | tr ";" "\t" )
#           detjobids=${detjobids}${nl}$report
#        done
        echo -e "$MSG\n\nDetails:\n\n$LOGS\n$detjobids\n\nPlease view $outputdir/logs/Summary.Report" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        echo -e "$MSG\n\nDetails:\n\n$LOGS\n$detjobids" > $outputdir/logs/Summary.Report

fi
