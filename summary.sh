#!/bin/bash
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
        genderinfo=$( cat $runfile | grep -w GENDERINFORMATION | cut -d '=' -f2 )
        sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )

        echo -e "the delivery folder should be populated already with BAMs and VCFs"
        if [ `expr ${#deliveryfolder}` -lt 2 ]
        then
            delivery=$outputdir/delivery
        else
	    delivery=$outputdir/$deliveryfolder
	fi

        echo -e "populating the delivery/docs folder with documents runfiles etc"
        
        mkdir -p ${delivery}/docs
        cp $outputdir/*.txt ${delivery}/docs
        cp $outputdir/*.list ${delivery}/docs
        cp $sampleinfo  ${delivery}/docs
        cp $genderinfo  ${delivery}/docs

        chmod -R 660 ${delivery}/docs

        echo -e "now putting together the second part of the Summary.Report file with the list of jobs executed inside this pipeline"
        
	listjobids=$( cat $outputdir/logs/pbs.* cat $outputdir/logs/*/pbs.* | sort | uniq | tr "\n" "\t" )
	pipeid=$( cat $outputdir/logs/pbs.CONFIGURE )

        if [ $exitstatus == "exitok" ]
        then
            MSG="Variant calling workflow with id: [$pipeid] by username: $USER finished with ALL  jobs with exit code 0 at: "$( echo `date` )
        else
            MSG="Variant calling workflow with id: [$pipeid] by username: $USER finished with SOME jobs with exit code 0  at: "$( echo `date` )
        fi
        
        LOGS="Results and execution logs can be found at \n$outputdir\n\nJOBIDS\n\n$listjobids\n\nThis jobid:${PBS_JOBID}\n\n"
        echo -e "$MSG\n\nDetails:\n\n$LOGS\n$detjobids\n\nPlease view $outputdir/logs/Summary.Report" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        echo -e "$MSG\n\nDetails:\n\n$LOGS\n$detjobids" >> $outputdir/logs/Summary.Report
        cp  $outputdir/logs/Summary.Report ${delivery}/docs/Summary.Report      
fi
