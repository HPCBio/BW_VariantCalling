#!/bin/sh

########################### 
#		$1		=	       run info file
###########################
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 1 ]
then
        MSG="Parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
        exit 1;
else
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found."
           exit 1;
        fi

        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "Program $0 stopped.\n\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi
        if [ -z $email ]
        then
           MSG="Invalid value for parameter PBSEMAIL=$email in configuration file"
           echo -e "Program $0 stopped.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi

        if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for parameter INPUTTYPE=$input_type  in configuration file."
                echo -e "Program $0 stopped.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
            fi
        fi

        if [ ! -d $outputdir ]
        then
            mkdir -p $outputdir
            mkdir -p $outputdir/logs
        else 
            echo "resetting directory"
	    `rm -r $outputdir/*`
            mkdir -p $outputdir/logs
        fi
        `chmod -R 770 $outputdir/`
        `chmod 750 $epilogue`
	`cp $runfile $outputdir/runfile.txt`
        runfile=$outputdir/runfile.txt


        # initialize autodocumentation script
        truncate -s 0 $outputdit/WorkflowAutodocumentationScript.sh
        echo "#!/bin/bash" > $outputdit/WorkflowAutodocumentationScript.sh
        WorkflowName=`basename $outputdir`  
        echo "# @begin $WorkflowName" > $outputdit/WorkflowAutodocumentationScript.sh


        outputlogs=$outputdir/logs
	echo "launching the pipeline configuration script"
        qsub1=$outputlogs/qsub.configure
        echo "#PBS -V" > $qsub1
        echo "#PBS -A $pbsprj" >> $qsub1
        echo "#PBS -N CONFIGURE" >> $qsub1
        echo "#PBS -l epilogue=$epilogue" >> $qsub1
	echo "#PBS -l walltime=00:03:00" >> $qsub1
	echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	echo "#PBS -o $outputlogs/CONFIGURE.ou" >> $qsub1
	echo "#PBS -e $outputlogs/CONFIGURE.in" >> $qsub1
        echo "#PBS -q debug" >> $qsub1
        echo "#PBS -m ae" >> $qsub1
        echo "#PBS -M $email" >> $qsub1
        echo "$scriptdir/configure.sh $runfile batch $outputlogs/CONFIGURE.in $outputlogs/CONFIGURE.ou $email $outputlogs/qsub.configure" >> $qsub1
        `chmod a+r $qsub1`               
        jobid=`qsub $qsub1`
        pipeid=$( echo $jobid | sed "s/\.[a-z]*[0-9]*//g" )
        echo $pipeid >> $outputlogs/CONFIGUREpbs
        echo `date`

        MSG="Variant calling workflow with id:[${pipeid}] started by username:$USER at: "$( echo `date` )
        LOGS="jobid=${jobid}\nqsubfile=$outputlogs/qsub.configure\nrunfile=$outputdir/runfile.txt\nerrorlog=$outputlogs/CONFIGURE.in\noutputlog=$outputlogs/CONFIGURE.ou"
        echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"


fi
