#!/bin/sh
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
        runfile=$1
        email=$2
        exitstatus=$3
        reportticket=$4
        numdays=30
        outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )

        if [ `expr ${#deliveryfolder}` -lt 2 ]
        then
            delivery=$outputdir/delivery
        else
	    delivery=$outputdir/$deliveryfolder
	fi

        echo -e "populating the delivery folder with documents runfiles etc"
        cp $outputdir/*.txt $delivery
        cp $outputdir/*.list $delivery

	listjobids=$( cat $outputdir/logs/pbs.* cat $outputdir/logs/*/pbs.* | sort | uniq | tr "\n" "\t" )
	pipeid=$( cat $outputdir/logs/pbs.CONFIGURE )

        if [ $exitstatus == "exitok" ]
        then
            MSG="Variant calling workflow with id: [$pipeid] by username: $USER finished with ALL  jobs with exit code 0 at: "$( echo `date` )
        else
            MSG="Variant calling workflow with id: [$pipeid] by username: $USER finished with SOME jobs with exit code 0  at: "$( echo `date` )
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
        echo -e "$MSG\n\nDetails:\n\n$LOGS\n$detjobids" > $delivery/Summary.Report

fi
