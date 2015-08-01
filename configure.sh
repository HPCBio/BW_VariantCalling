#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 6 ]
then
        MSG="Parameter mismatch."
        echo -e "Program $0 stopped. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
        exit 1;
else
        
        echo -e "\n\n############# CONFIGURE VARIANT CALLING WORKFLOW ###############\n\n" >&2
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        runmode=$2
        elog=$3
        olog=$4
        email=$5
        qsubfile=$6
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"





        set +x; echo -e "\n\n############# CHECKING PARAMETERS ###############\n\n" >&2; set -x;
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found."
           exit 1;
        fi

 
        set +x; echo -e "\nadditional variable assignment from runfile and sanity check\n" >&2; set -x;

        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        resortbam=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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
		MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
                echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
            fi
        fi

        if [ -z $thr -o -z $outputdir -o -z $pbsprj -o -z $epilogue ]
        then
 		MSG="Invalid value specified for any of these paramaters in configuration file:\nPBSTHREADS=$thr\nOUTPUTDIR=$outputdir\nPBSPROJECTID=$pbsprj\nEPILOGUE=$epilogue"
                echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
        fi

        if [ $inputformat != "FASTQ" -a $inputformat != "BAM" ]
        then
            MSG="Incorrect value for INPUTFORMAT=$inputformat in the configuration file."
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi


        if [ $resortbam != "1" -a $resortbam != "0" -a $resortbam != "YES" -a $resortbam != "NO" ]
        then
           MSG="Invalid value for RESORTBAM=$resortbam"
           echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
            if [ $resortbam == "1" ]
            then
                $resortbam="YES"
            fi
            if [ $resortbam == "0" ]
            then
                $resortbam="NO"
            fi
        fi

        if [ $bamtofastqflag != "YES" -a $bamtofastqflag != "NO" -a $bamtofastqflag != "1" -a $bamtofastqflag != "0" ]
        then
            MSG="BM2FASTQFLAG=$bamtofastqflag  invalid value"
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        else
            if [ $bamtofastqflag == "1" ]
            then
                $bamtofastqflag="YES"
            fi
            if [ $bamtofastqflag == "0" ]
            then
                $bamtofastqflag="NO"
            fi
        fi
        if [ $multisample != "YES" -a $multisample != "NO" -a $multisample != "1" -a $multisample != "0" ]
        then
            MSG="MULTISAMPLE=$multisample  invalid value"
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        else
            if [ $multisample == "1" ]
            then
                $multisample="YES"
            fi
            if [ $multisample == "0" ]
            then
                $multisample="NO"
            fi
        fi


        if [ $resortbam == "YES" -a $bamtofastqflag == "YES" ]
        then
            MSG="Incompatible values for the pair RESORTBAM=$resortbam and BAM2FASTQFLAG=$bam2fqflag in the configuration file."
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi

        if [ $multisample == "NO" -a $analysis == "MULTIPLEXED" ]
        then
            MSG="Incompatible values for the pair MULTISAMPLE=$multisample and ANALYSIS=$analysis in the configuration file."
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        elif [ $multisample == "YES" -a $analysis == "MULTIPLEXED" ]
        then 
            set +x; echo -e "\nChecking that a tab delimited file with information about samples and lanes must exist for the MULTIPLEXED analysis pipeline to start.\n" >&2; set -x;
            if [ ! -s $sampleinfo ]
            then
		MSG="SAMPLEINFORMATION=$sampleinfo invalid value. A tab delimited file with lanes and samples must be specified. "
                echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
            fi
        fi

        if [ ! -d $scriptdir ]
        then
            MSG="SCRIPTDIR=$scriptdir directory not found"
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
        if [ ! -d $outputdir ]
        then
            mkdir -p $outputdir
            mkdir -p $outputdir/logs
	elif [ $runmode != "batch" ]
        then
	    set +x; echo -e "\n resetting logs\n" >&2; set -x; 
	    `rm -r $outputdir/logs/*`
        fi
	`chmod -R 770 $outputdir`
        `chmod 740 $epilogue`
        TopOutputLogs=$outputdir/logs
        pipeid=$( cat $TopOutputLogs/pbs.CONFIGURE )





        #############################
	set +x; echo -e "\n\nconstructing files with list(s) of input files to analyze in this run of the pipeline" >&2;
        echo -e "IF input is NOT multiplexed, then ONE file will be created, otherwise THREE files will be created\n\n" >&2; set -x;
 
	if [ $analysis == "MULTIPLEXED" ]
	then
            set +x; echo -e "\n\n ############# ANALYSIS IS MULTIPLEXED ! ############\n\n" >&2; 
            echo -e "produce SAMPLENAMES.list SAMPLENAMES_multiplexed.list SAMPLEGROUPS.list " >&2
            echo -e "from Baylor's info sheet specified in runfile in line SAMPLEINFORMATION=$sampleinfo." >&2; set -x;
	    perl $scriptdir/Baylor2SAMPLENAMES.pl $outputdir $sampleinfo SAMPLENAMES.list SAMPLENAMES_multiplexed.list SAMPLEGROUPS.list
            
            set +x; echo -e "\n ### testing that the perl script actually worked ### \n"  >&2; set -x;
	    if [ ! -s $outputdir/SAMPLENAMES.list ]
	    then
		MSG="$outputdir/SAMPLENAMES.list is empty"
                 echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    fi
	    if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
	    then
		MSG="$outputdir/SAMPLENAMES_multiplexed.list is empty"
                echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    fi
	    if [ ! -s $outputdir/SAMPLEGROUPS.list ]
	    then
		MSG="$outputdir/SAMPLEGROUPS.list is empty"
                echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    fi
	elif [ ! -e $outputdir/SAMPLENAMES.list ]
	then
            set +x; echo -e "\n ###### ANALYSIS = $analysis, IS NOT = MULTIPLEXED ########" >&2;
            echo -e "\n ###### produce just one file with list of files in directory specified in run file in line INPUTDIR=$sampledir ########" >&2; set -x;
	    truncate -s 0 $outputdir/SAMPLENAMES.tmp.list
	    if  [ $inputformat == "FASTQ" ]
	    then
		for inputfile in $sampledir/*
		do
                   # strip path, which read (left/right), and extension from input files
                   # and put that info into the SampleNames file
		   SampleName=$( basename $inputfile | sed 's/_read.\?\..*$//' )
                   echo -e "$SampleName" >> $outputdir/SAMPLENAMES.tmp.list
		done
	    elif [ $inputformat == "BAM" ]
	    then
		for inputfile in $sampledir/*
		do
                   # strip path, which read (left/right), and extension from input files
                   # and put that info into the SampleNames file
		   SampleName=$( basename $inputfile .bam )
                   echo -e "$SampleName" >> $outputdir/SAMPLENAMES.tmp.list
		done
	    fi
            # paired-ended fastq will produce duplicate lines in the SampleNames file, so remove the duplicates
	    sort $outputdir/SAMPLENAMES.tmp.list | uniq > $outputdir/SAMPLENAMES.list
	    sed -i '/^\s*$/d' $outputdir/SAMPLENAMES.list # remove blank lines
	    rm  $outputdir/SAMPLENAMES.tmp.list
	fi
        set +x; echo -e "\n### check that this actually worked, " >&2
        echo -e "### because otherwise the bash script will just go on, as if there is no problem \n"  >&2; set -x;
	if [ ! -s $outputdir/SAMPLENAMES.list ]
	then
	    MSG="$outputdir/SAMPLENAMES.list is empty"
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit 1;
	fi

	numsamples=`wc -l $outputdir/SAMPLENAMES.list | cut -d ' ' -f 1`
	if [ $numsamples -lt 1 ]
	then
	    MSG="No samples found in INPUTDIR=$sampledir."
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit 1;
	fi




        ###################################################
	set +x; echo -e "\n\nOne more configuration task: generate a qsub header so we would not have to repeat the same lines\n\n" >&2; set -x;
	generic_qsub_header=$outputdir/qsubGenericHeader
	truncate -s 0 $generic_qsub_header
	echo "#!/bin/bash" > $generic_qsub_header
	echo "#PBS -V" >> $generic_qsub_header
	echo "#PBS -A $pbsprj" >> $generic_qsub_header
	echo "#PBS -q $pbsqueue" >> $generic_qsub_header
	echo "#PBS -m ae" >> $generic_qsub_header
	echo "#PBS -M $email" >> $generic_qsub_header
        set +x; echo -e "\n### check that this actually worked, " >&2
        echo -e "### because otherwise the bash script will just go on, as if there is no problem \n"  >&2; set -x;
        if [ ! -s $generic_qsub_header ]
	then 
	    MSG="$generic_qsub_header is empty"
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit 1;
	fi





        ###################################################
        set +x; 
        echo -e "\n\n"; 
        echo -e "done with preprocessing and configuration steps" >&2
        echo -e "now we select analysis or analyses to run" >&2
        echo -e "based on the value specified in runfile in line ANALYSIS"  >&2 
        echo -e "\n\n" >&2; set -x;


	case=""
	    if [ $analysis == "ALIGN" -o $analysis == "ALIGNMENT" ]
	    then
		echo "Type of analysis to run: ALIGNMENT only"      
		qsub1=$TopOutputLogs/qsub.start_align_block
                echo "#!/bin/bash" > $qsub1
		echo "#PBS -A $pbsprj" >> $qsub1
		echo "#PBS -N ${pipeid}_START_ALIGN_BLOCK" >> $qsub1
		echo "#pbs -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=00:30:00" >> $qsub1
		echo "#PBS -l nodes=1:ppn=1" >> $qsub1
		echo "#PBS -o $TopOutputLogs/log.start_align_block.ou" >> $qsub1
		echo "#PBS -e $TopOutputLogs/log.start_align_block.in" >> $qsub1
		echo "#PBS -q $pbsqueue" >> $qsub1
		echo "#PBS -m ae" >> $qsub1
		echo "#PBS -M $email" >> $qsub1
		echo "$scriptdir/start_align_block.sh $runfile $TopOutputLogs/log.start_align_block.in $TopOutputLogs/log.start_align_block.ou $email $TopOutputLogs/qsub.start_align_block" >> $qsub1
		`chmod a+r $qsub1`               
		`qsub $qsub1 >> $TopOutputLogs/pbs.ALIGN`
		echo `date`
		case="alignonly"
	    fi
	    if [ $analysis == "REALIGNONLY" -o $analysis == "REALIGN_ONLY" ]
	    then
		echo "Type of analysis to run: REALIGNMENT only. bams provided"
		qsub2=$TopOutputLogs/qsub.start_realrecal_block
                echo "#!/bin/bash" > $qsub2
		echo "#PBS -A $pbsprj" >> $qsub2
		echo "#PBS -N ${pipeid}_START_REALRECAL_BLOCK"
		echo "#pbs -l epilogue=$epilogue" >> $qsub2
		echo "#PBS -l walltime=00:30:00" >> $qsub2
		echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		echo "#PBS -o $TopOutputLogs/log.start_realrecal_block.ou" >> $qsub2
		echo "#PBS -e $TopOutputLogs/log.start_realrecal_block.in" >> $qsub2
		echo "#PBS -q $pbsqueue" >> $qsub2
		echo "#PBS -m ae" >> $qsub2
		echo "#PBS -M $email" >> $qsub2
		echo "$scriptdir/start_realrecal_block.sh $runfile $TopOutputLogs/log.start_realrecal_block.in $TopOutputLogs/log.start_realrecal_block.ou $email $TopOutputLogs/qsub.start_realrecal_block" >> $qsub2
		`chmod a+r $qsub2` 
		`qsub $qsub2 >> $TopOutputLogs/pbs.REALRECAL`
		echo `date`
		case="realignonly" 
	    fi
	    if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" -o $analysis == "MULTIPLEXED" ]
	    then
		echo "Type of analysis to run: ALIGNMENT and REALIGNMENT"
		qsub1=$TopOutputLogs/qsub.start_align_block
                echo "#!/bin/bash" > $qsub1
		echo "#PBS -A $pbsprj" >> $qsub1
		echo "#PBS -N ${pipeid}_START_ALIGN_BLOCK" >> $qsub1
		echo "#pbs -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=00:30:00" >> $qsub1
		echo "#PBS -l nodes=1:ppn=1" >> $qsub1
		echo "#PBS -o $TopOutputLogs/log.start_align_block.ou" >> $qsub1
		echo "#PBS -e $TopOutputLogs/log.start_align_block.in" >> $qsub1
		echo "#PBS -q $pbsqueue" >> $qsub1
		echo "#PBS -m ae" >> $qsub1
		echo "#PBS -M $email" >> $qsub1
		echo "$scriptdir/start_align_block.sh $runfile $TopOutputLogs/log.start_align_block.in $TopOutputLogs/log.start_align_block.ou $email $TopOutputLogs/qsub.main.aln" >> $qsub1
		`chmod a+r $qsub1`               
		`qsub $qsub1 >> $TopOutputLogs/pbs.ALIGN`
		echo `date`
		echo "Note: realign module will be scheduled after align module ends"
		case="align and realign"  
	    fi
	    if [ $analysis == "VCALL_ONLY" -o $analysis == "VCALL" ]
	    then
		echo "variant calling only"
		qsub3=$TopOutputLogs/qsub.start_varcall_block
		echo "#!/bin/bash" > $qsub3
		echo "#PBS -A $pbsprj" >> $qsub3
		echo "#PBS -N ${pipeid}_START_VARCALL_BLOCK" >> $qsub3
		echo "#PBS -l epilogue=$epilogue" >> $qsub3
		echo "#PBS -l walltime=00:30:00" >> $qsub3
		echo "#PBS -l nodes=1:ppn=1" >> $qsub3
		echo "#PBS -o $TopOutputLogs/log.start_varcall_block.ou" >> $qsub3
		echo "#PBS -e $TopOutputLogs/log.start_varcall_block.in" >> $qsub3
		echo "#PBS -q $pbsqueue" >> $qsub3
		echo "#PBS -m ae" >> $qsub3
		echo "#PBS -M $email" >> $qsub3
		echo "$scriptdir/start_varcall_block.sh $runfile $TopOutputLogs/log.start_varcall_block.in $TopOutputLogs/log.start_varcall_block.ou $email $TopOutputLogs/qsub.start_varcall_block" >> $qsub3
		`chmod a+r $qsub3`
		vcalljobid=`qsub $qsub3`
		echo $vcalljobid >> $TopOutputLogs/pbs.VARCALL
		case="vcall_only"  
	    fi
	    if [ `expr ${#case}` -lt 1 ]
	    then
		MSG="Invalid value for parameter ANALYSIS=$analysis in configuration file."
                echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1; 
	    fi
fi
