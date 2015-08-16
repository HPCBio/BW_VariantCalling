#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 3 ]
        MSG="Parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`	
        scriptfile=$0
        outputdir=$1
        email=$2
        exitstatus=$3
        
	listjobids=$( cat $outputdir/logs/pbs.* | sed "s/\.[a-z]*//g" | tr "\n" "\t" )        
        if [ $exitstatus == "exitok" ]
        then
            MSG="GGPS pipeline finished successfully on iforge by username:$USER at: "$( echo `date` )
        else
            MSG="GGPS pipeline finished on iforge by username:$USER at: "$( echo `date` )
        fi
        LOGS="Results and execution logs can be found on iforge at $outputdir\n\nJOBIDS\n\n$listjobids\n\nThis jobid:${PBS_JOBID}\n\n"
        detjobids=""
        nl="\n"
        for jobid in $listjobids
        do

           cputime=`tracejob -q -n 10 $jobid | grep -m 1 "resources_used.walltime=" | sed 's/^.*resources_used.walltime=//' |  awk 'BEGIN { FS=":"} { print $1*3600+$2*60+$3}'`
           price=$( echo "scale=2; 16*$cputime/3600*0.216" | bc )

           report=`checkjob -A $jobid`
           report=$( echo "$report\tAPPROX_PRICE=\$$price" | tr ";" "\t" )
           detjobids=${detjobids}${nl}$report


           echo $report

        done
        echo -e "$MSG\n\nDetails:\n\n$LOGS\n$detjobids" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""

fi
