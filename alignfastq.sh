#!/bin/bash
#
# alignfastq.sh
# align module to be used for input files in fastq format
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

        echo -e "\n\n############# BEGIN ALIGNFASTQ PROCEDURE: schedule fastqc, parse sample information and create alignment jobs  ###############\n\n" >&2

	set -x
	echo `date`
        scriptfile=$0
        runfile=$1
        elog=$2
        olog=$3
        email=$4
        qsubfile=$5
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        set +x; echo -e "\n\n############# CHECKING PARAMETERS ###############\n\n" >&2; set -x;
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found."
           echo -e "Program $0 stopped. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
           exit 1;
        fi


set +x; echo -e "\n\n" >&2; 
# wrapping commends in echoes, so that the output logs would be easier to read: they will have more structure
echo "####################################################################################################" >&2
echo "##################################### PARSING RUN INFO FILE ########################################" >&2
echo "##################################### AND SANITY CHECK      ########################################" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
        run_method=$( cat $runfile | grep -w RUNMETHOD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        chunkfastq=$( cat $runfile | grep -w CHUNKFASTQ | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        aligner=$( cat $runfile | grep -w ALIGNER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        fastqcdir=$( cat $runfile | grep -w FASTQCDIR | cut -d '=' -f2 )
        fastqcflag=$( cat $runfile | grep -w FASTQCFLAG | cut -d '=' -f2 )
        fastqcparms=$( cat $runfile | grep -w FASTQCPARMS | cut -d '=' -f2 | tr " " "_" )_-t_${thr}
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        samblasterdir=$( cat $runfile | grep -w SAMBLASTERDIR | cut -d '=' -f2 )
	dup=$( cat $runfile | grep -w MARKDUP  | cut -d '=' -f2 )
        dupflag=$( cat $runfile | grep -w REMOVE_DUP  | cut -d '=' -f2 )
        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        markduplicatestool=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        profiling=$( cat $runfile | grep -w PROFILING | cut -d '=' -f2 )
################## MUST PUT IN A CHECK FOR CORRECT SETTING ON THE PROFILING VARIABLE
        profiler=$( cat $runfile | grep -w PROFILER | cut -d '=' -f2 )
        cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )


        set +x; echo -e "\n\n\n############ checking workflow scripts directory\n" >&2; set -x;

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi
        if [ ! -d $outputdir ]
        then
           mkdir -p $outputdir
        fi


        set +x; echo -e "\n\n\n############ checking input type: WGS or WES\n" >&2; set -x

        if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
                echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
            fi
        fi


        set +x; echo -e "\n\n\n# check that sampledir exists\n" >&2; set -x
        if [ ! -d $sampledir ]
        then
           MSG="SAMPLEDIR=$sampledir directory not found"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi


        set +x; echo -e "\n\n\n# check number of samples\n" >&2; set -x
        numsamples=`wc -l $outputdir/SAMPLENAMES.list | cut -d ' ' -f 1`
        if [ $numsamples -gt 1 -a $multisample == "YES" ]
        then
            echo "multiple samples to be aligned."
        else
           if [ $numsamples -eq 1 -a $multisample == "NO" ]
           then
              echo "single sample to be aligned."
           fi
        fi


        set +x; echo -e "\n\n\n# check whether fastq will be chunked\n" >&2; set -x
        if [ $chunkfastq != "YES" -a $chunkfastq != "NO" ]
        then
            MSG="CHUNKFASTQ variable must be binary YES/NO; incorrect value encountered"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi


        set +x; echo -e "\n\n\n############ checking settings for marking of duplicates\n" >&2; set -x
        if [ $dup != "1" -a $dup != "0" -a $dup != "YES" -a $dup != "NO" ]
        then
           MSG="Invalid value for MARKDUP=$dup"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
            if [ $dup == "1" ]
            then
                $dup="YES"
            fi
            if [ $dup == "0" ]
            then
                $dup="NO"
            fi
        fi
        if [ $dupflag != "1" -a $dupflag != "0" -a $dupflag != "YES" -a $dupflag != "NO" ]
        then
           MSG="Invalid value for REMOVE_DUP=$dupflag"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
            if [ $dupflag == "1" ]
            then
                $dupflag="YES"
            fi
            if [ $dupflag == "0" ]
            then
                $dupflag="NO"
            fi
        fi
	dupparms=$( echo "dup=${dup}_flag=${dupflag}" )


        set +x; echo -e "\n\n\n############ checking FastQC settings\n" >&2; set -x 
        if [ $fastqcflag != "1" -a $fastqcflag != "0" -a $fastqcflag != "YES" -a $fastqcflag != "NO" ]
        then
           MSG="Invalid value for FASTQCFLAG=$fastqcflag"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi
        if [ $fastqcflag == "1" ]
        then
            $fastqcflag="YES"
        fi
        if [ $fastqcflag == "0" ]
        then
            $fastqcflag="NO"
        fi


        set +x; echo -e "\n\n\n############ checking Cleanup settings\n" >&2; set -x
        if [ $cleanupflag != "1" -a $cleanupflag != "0" -a $cleanupflag != "YES" -a $cleanupflag != "NO" ]
        then
           MSG="Invalid value for REMOVETEMPFILES=$cleanupflag"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
            if [ $cleanupflag == "1" ]
            then
                $cleanupflag="YES"
            fi
            if [ $cleanupflag == "0" ]
            then
                $cleanupflag="NO"
            fi
        fi


        set +x; echo -e "\n\n\n############ checking computational tools\n" >&2; set -x
        if [ $aligner != "NOVOALIGN" -a $aligner != "BWA_ALN" -a $aligner != "BWA_MEM"]
        then
            MSG="ALIGNER=$aligner  is not available at this site"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi
        if [ $aligner == "NOVOALIGN" ]
        then
            alignerdir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w NOVOINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w NOVOPARAMS | cut -d '=' -f2 | tr " " "_" )
        fi
        if [ $aligner == "BWA_ALN" ]
        then
            alignerdir=$( cat $runfile | grep -w BWAALNDIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w BWAALNINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w BWAALNPARAMS | cut -d '=' -f2 | tr " " "_" )
        fi
        if [ $aligner == "BWA_MEM" ]
        then
            alignerdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w BWAMEMINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 | tr " " "_" )
        fi
        if [ -z $sortool ]
        then
           MSG="Value for SORTOOL must be specified in configuration file"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
           if [ $sortool != "NOVOSORT" -a $sortool != "PICARD" ]
           then
               MSG="Invalid value for SORTOOL=$sortool in configuration file"
               echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit 1;
           fi
        fi      
        if [ -z $markduplicatestool ]
        then
           MSG="Value for MARKDUPLICATESTOOL must be specified in configuration file"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
           if [ $markduplicatestool != "PICARD" -a $markduplicatestool != "SAMBLASTER" ]
            then
               MSG="Invalid value for SORTOOL=$sortool in configuration file"
               echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit 1;
           fi
        fi

        if [ ! -d $alignerdir ]
        then
           MSG="$alignerdir aligner directory not found"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi

        if [ ! -d $picardir ]
        then
           MSG="PICARDIR=$picardir directory not found"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi

        if [ ! -d $samdir ]
        then
           MSG="SAMDIR=$samdir directory not found"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi

        if [ ! -e $profiler ]
        then
           MSG="PROFILER=$profiler not found"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi


        set +x; echo -e "\n\n\n############ checking presence of references\n" >&2; set -x
        if [ ! -s $refdir/$ref ]
        then
           MSG="$refdir/$ref reference genome not found"
           echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi


       set +x; echo -e "\n\n\n############ checking run method\n" >&2; set -x
       if [ $run_method != "LAUNCHER" -a $run_method != "QSUB" -a $run_method != "APRUN" -a $run_method != "SERVER" ]
       then
          MSG="Invalid value for RUNMETHOD=$run_method"
          echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
          exit 1;
       fi


set +x; echo -e "\n\n">&2
echo "############################################################################################################" >&2
echo "#####################################  CREATE  DIRECTORY STRUCTURE               ###########################" >&2
echo "############################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


        #AlignOutputDir=$outputdir/align

        # initialize output file name and path template for YesWorkflow
        AlignedFastqPathTemplate="align/"


        TopOutputLogs=$outputdir/logs
        if [ -d $TopOutputLogs ]
        then
           # not sure if we want to reset this: useful for debugging ....
           #set +x; echo -e "$TopOutputLogs is there; resetting it" >&2; set -x;
           #`rm -r $TopOutputLogs/*`
           pbsids=""
        else
           mkdir -p $TopOutputLogs
        fi

        AlignOutputLogs=$TopOutputLogs/align
        if [ ! -d $AlignOutputLogs ]
        then
            mkdir $AlignOutputLogs
        fi
        `chmod -R 770 $AlignOutputLogs/`
        # where messages about failures will go
        truncate -s 0 $AlignOutputLogs/FAILEDmessages
        if [ ! -d $AlignOutputLogs/FAILEDjobs ]
        then
            mkdir $AlignOutputLogs/FAILEDjobs
        else
            rm -r $AlignOutputLogs/FAILEDjobs/*
        fi
        `chmod -R 770 $AlignOutputLogs/FAILEDjobs`




        pipeid=$( cat $TopOutputLogs/pbs.CONFIGURE )

        #chunks=`expr $nodes "-" 1`
        #if [ $chunks -lt 1 ]
        #then
        #    chunks=$nodes
        #fi
        nthreads=`expr $thr "-" 1`
        if [ $nthreads -lt 1 ]
        then
            nthreads=$thr
        fi


set +x; echo -e "\n\n">&2
echo "############################################################################################################" >&2
echo "#####################################                               ########################################" >&2
echo "##################################### ALIGNMENT: LOOP1 OVER SAMPLES ########################################" >&2
echo "#####################################                               ########################################" >&2
echo "############################################################################################################" >&2
echo "######           select the file to read sample info from. This file was created in configuration.sh  ######" >&2
echo -e "\n\n" >&2; set -x;



        # clear the joblist
        truncate -s 0 $AlignOutputLogs/AlignAnisimov.joblist
        
        # select the file to read
        if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
        then 
             TheInputFile=$outputdir/SAMPLENAMES_multiplexed.list
        else
             TheInputFile=$outputdir/SAMPLENAMES.list
        fi


        while read SampleLine
        do
          set +x; echo -e "\n\n processing next sample \n" >&2; set -x;
          # this will evaluate the length of string
          if [ `expr ${#SampleLine}` -lt 1 ]
          then
             set +x; echo -e "\n\n skipping empty line \n" >&2; set -x;
          else

            
            set +x; echo -e "\n\n">&2
	    echo -e "#################################### PREP WORK FOR $SampleLine            ########################################\n" >&2
	    echo -e "#################################### PARSING READS and CREATING RGLINE    ########################################\n" >&2
            echo -e "\n\n" >&2; set -x;


            if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
	    then
                set +x; echo -e "\n\nSAMPLENAMES_multiplexed.list DOES NOT exist." >&2
		echo -e "\n ###### Current line has ONE field for SAMPLENAME." >&2;
		echo -e "\n ###### The code parses PL CN and LB values from runfile\n" >&2; set -x;
                # form and check the left reads file
                SampleName=$( echo $SampleLine )
		LeftReadsFastq=$( ls $sampledir/${SampleName}_read1.* )
		RightReadsFastq=$( ls $sampledir/${SampleName}_read2.* )

                # form the read group values for the aligner 
		sID=$SampleName
		sPU=$SampleName
		sSM=$SampleName
		sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
		sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
		sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
	    else
                set +x; echo -e "\n" >&2
                echo -e "\n ###### SAMPLENAMES_multiplexed.list DOES exist." >&2
                echo -e "\n ###### Current line has FIVE fields:" >&2
                echo -e "\n ###### col1=sampleid col2=flowcell_and_lane_name col3=lib col4=read1 col5=read2." >&2
		echo -e "\n ###### The code parses PL and CN from runfile\n" >&2; 
                echo -w "\n" >&2; set -x;

                SampleName=$( echo -e "$SampleLine" | cut -f 2 )
                LeftReadsFastq=${sampledir}/$( echo -e "$SampleLine" | cut -f 4 )
                RightReadsFastq=${sampledir}/$( echo -e "$SampleLine" | cut -f 5 )
		sID=$( echo -e "$SampleLine" | cut -f 2 )
		sPU=$( echo -e "$SampleLine" | cut -f 2 )
		sSM=$( echo -e "$SampleLine" | cut -f 1 )
		sLB=$( echo -e "$SampleLine" | cut -f 3 )
		sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
		sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
	    fi

            set +x; echo -e "\n\n checking that it actually worked: we must have RG line and file names of the reads... \n" >&2; set -x;

            if [ `expr ${#sID}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ] 
	    then
		MSG="ID=$sID PL=$sPL CN=$sCN invalid values. The RG line cannot be formed"
                echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    fi


	    RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )

	    if [ ! -s $LeftReadsFastq ]
	    then
                MSG="$LeftReadsFastq left reads file not found"
                echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
	    fi
	    echo `date`


	    if [ $paired -eq 1 ]
	    then
		if [ ! -s $RightReadsFastq ]
		then
		    MSG="$RightReadsFastq right reads file not found"
                    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		    exit 1;
		fi
	    fi


            set +x; echo -e "\n\n ################    CREATING OUTPUT FOLDERS / OUTPUT FILENAMES   ################\n" >&2; set -x;

            AlignOutputDir=$outputdir/$SampleName/align
            if [ -d $AlignOutputDir ]
            then
                # perhaps we should stop when a sample is seen more than once. But for now, we simply reset the folder
		set +x; echo -e "\n\n $AlignOutputDir is there; resetting it" >&2; set -x;
		`rm -r $AlignOutputDir/*`
            else
		mkdir -p $AlignOutputDir/logs
            fi
            `chmod -R 770 $AlignOutputDir/`
            `chmod -R 770 $AlignOutputDir/logs`

            # filenames of temporary files go here 
            outputsamfileprefix=$AlignOutputDir/$SampleName
            sortedplain=$outputsamfileprefix.wrg.sorted.bam
            outsortnodup=$outputsamfileprefix.nodups.sorted.bam
            outsortwdup=$outputsamfileprefix.wdups.sorted.bam



            set +x; echo -e "\n\n ###################################   END PREP WORK  ###################################### " >&2
	    echo -e "\n\n #######                       PRE-ALIGNMENT BLOCK BEGINS                            ####### " >&2
            echo -e "\n\n #######          CHECKING IF WE NEED TO LAUNCH PRE-ALIGNMENT TASK(S)                ####### " >&2; set -x;

            if [ $fastqcflag == "YES" ]
            then
                set +x; echo -e "\n\n" >&2; 
		echo "#################################### FASTQFLAG == YES RUNNING FASTQC ON UNALIGNED READS of $SampleName #####################" >&2
                echo -e "\n\n" >&2; set -x;

                # check that fastqc output folder for this sample is there
		FastqcOutputFolder=$outputdir/$SampleName/fastqc
		if [ -d $FastqcOutputFolder ]
		then
		    set +x; echo -e "\n\n $FastqcOutputFolder  is there; resetting it" >&2; set -x;
		    `rm -r $FastqcOutputFolder/*`
		else
		    mkdir -p $FastqcOutputFolder
		fi
		`chmod -R 770 $FastqcOutputFolder`


                # check that fastqc tool is there
		if [ ! -d $fastqcdir ]
		then
		    MSG="FASTQCDIR=$fastqcdir directory not found"
                    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		    exit 1;
		fi

                # create necessary log directories
                if [ ! -d $TopOutputLogs/fastqc ]
                then
                    mkdir $TopOutputLogs/fastqc
                    `chmod -R 770 $TopOutputLogs/fastqc`
                fi



                set +x; echo -e "\n############# Now that we know fastq will be quality-checked, " >&2
                echo -e "############# we need to initialize Workflow autodocumentation\n" >&2; set -x

                set +x; echo -e "### update autodocumentation script ###"; set -x;
                echo -e "# @begin FastQC" >> $outputdir/WorkflowAutodocumentationScript.sh
                # select how samples will be found
                #if [[ -s $outputdir/SAMPLENAMES_multiplexed.list ] && [ -s $outputdir/SAMPLEGROUPS.list ]]
                #then
                #    numinputs=`wc -l $outputdir/SAMPLENAMES_multiplexed.list | cut -d ' ' -f 1`
                #    numsamplegroups=`wc -l $outputdir/SAMPLEGROUPS.list | cut -d ' ' -f 1`
                #    echo -e "   # @in datafilestoalign @as fastq_inputs @URI SAMPLENAMES_multiplexed.list=${numinputs}_inputs_for_${numsamplegroups}_samples" >> $outputdir/WorkflowAutodocumentationScript.sh
                #else
                #    numinputs=`wc -l $outputdir/SAMPLENAMES.list | cut -d ' ' -f 1`
                #    echo -e "   # @in datafilestoalign @as fastq_inputs @URI SAMPLENAMES.list=${numinputs}_inputs_for_${numinputs}_samples" >> $outputdir/WorkflowAutodocumentationScript.sh
                #fi
                FastqcOutputFolder_basename=`(basename $FastqcOutputFolder)`
                echo -e "   # @in fastq_inputs @as fastq_inputs  @URI ${SampleName}" >> $outputdir/WorkflowAutodocumentationScript.sh
                echo -e "   # @out fastq_output @as fastqc_output_folder @URI ${FastqcOutputFolder_basename}/" >> $outputdir/WorkflowAutodocumentationScript.sh
                echo -e "# @end FastQC\n" >> $outputdir/WorkflowAutodocumentationScript.sh

                # create and launch the qsub jobs with fastqc runs
                
                left_fastqc_input=$LeftReadsFastq                

                qsub_fastqcR1=$TopOutputLogs/fastqc/qsub.fastqcR1.$SampleName
		echo "#PBS -V" > $qsub_fastqcR1
		echo "#PBS -A $pbsprj" >> $qsub_fastqcR1
		echo "#PBS -N ${pipeid}_fastqc_R1_${SampleName}" >> $qsub_fastqcR1
		echo "#PBS -l walltime=$pbscpu" >> $qsub_fastqcR1
		echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_fastqcR1
		echo "#PBS -o $TopOutputLogs/fastqc/log.fastqc_R1_${SampleName}.ou" >> $qsub_fastqcR1
		echo "#PBS -e $TopOutputLogs/fastqc/log.fastqc_R1_${SampleName}.in" >> $qsub_fastqcR1
		echo "#PBS -q $pbsqueue" >> $qsub_fastqcR1
		echo "#PBS -m ae" >> $qsub_fastqcR1
		echo "#PBS -M $email" >> $qsub_fastqcR1
		echo "aprun -n 1 -d $thr $scriptdir/fastq.sh $fastqcdir $FastqcOutputFolder $fastqcparms $left_fastqc_input $TopOutputLogs/fastqc/log.fastqc_R1_${SampleName}.in $TopOutputLogs/log.fastqc_R1_${SampleName}.ou $email $TopOutputLogs/fastqc/qsub.fastqc_R1_$SampleName" >> $qsub_fastqcR1

                echo -e "\n\n" >> $qsub_fastqcR1
                echo "exitcode=\$?" >> $qsub_fastqcR1
                echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub_fastqcR1
                echo "   echo -e \"\n\n fastq.sh failed with exit code = \$exitcode \n logfile=$TopOutputLogs/fastqc/log.fastqc_R1_${SampleName}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub_fastqcR1
                echo -e "   exit 1" >> $qsub_fastqcR1
                echo "fi" >> $qsub_fastqcR1


		`chmod a+r $qsub_fastqcR1`
                `qsub $qsub_fastqcR1 >> $TopOutputLogs/pbs.FASTQC`

		if [ $paired -eq 1 ]
		then
                    right_fastqc_input=$RightReadsFastq


                    qsub_fastqcR2=$TopOutputLogs/fastqc/qsub.fastqc_R2_$SampleName
		    echo "#PBS -V" > $qsub_fastqcR2
		    echo "#PBS -A $pbsprj" >> $qsub_fastqcR2
		    echo "#PBS -N ${pipeid}_fastqc_R2_${SampleName}" >> $qsub_fastqcR2
		    echo "#PBS -l walltime=$pbscpu" >> $qsub_fastqcR2
		    echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_fastqcR2
		    echo "#PBS -o $TopOutputLogs/fastqc/log.fastqc_R2_${SampleName}.ou" >> $qsub_fastqcR2
		    echo "#PBS -e $TopOutputLogs/fastqc/log.fastqc_R2_${SampleName}.in" >> $qsub_fastqcR2
		    echo "#PBS -q $pbsqueue" >> $qsub_fastqcR2
		    echo "#PBS -m ae" >> $qsub_fastqcR2
		    echo "#PBS -M $email" >> $qsub_fastqcR2
		    echo "aprun -n 1 -d $thr $scriptdir/fastq.sh $fastqcdir $FastqcOutputFolder $fastqcparms $right_fastqc_input $TopOutputLogs/fastqc/log.fastqc_R2_${SampleName}.in $TopOutputLogs/fastqc/log.fastqc_R2_$SampleName.ou $email $TopOutputLogs/qsub.fastqc_R2_$SampleName" >> $qsub_fastqcR2

                    echo -e "\n\n" >> $qsub_fastqcR2
                    echo "exitcode=\$?" >> $qsub_fastqcR2
                    echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub_fastqcR2
                    echo "   echo -e \"\n\n fastq.sh failed with exit code = \$exitcode \n logfile=$TopOutputLogs/fastqc/log.fastqc_R2_${SampleName}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub_fastqcR2
                    echo -e "  exit 1" >> $qsub_fastqcR2
                    echo "fi" >> $qsub_fastqcR2

		    `chmod a+r $qsub_fastqcR2`
                    `qsub $qsub_fastqcR2 >> $TopOutputLogs/pbs.FASTQC`
		fi
            else
		set +x; echo -e "\n\n ############ FASTQCFLAG == NO. Quality information for fastq files will NOT be calculated." >&2; set -x;
            fi
            
            ## done with generating quality info for each read file


            set +x; echo -e "\n\n" >&2;
            echo -e "#################################### DONE WITH PRE-ALIGNMENT TASK(S)  ########################################\n\n" >&2
	    echo -e "#################################### SELECT TO CHUNK/Not2CHUNK READS for $SampleName  ########################################\n\n" >&2; set -x;

            # create new names for chunks of fastq
            # doing this outside the chunkfastq conditional, because
            # non-chunking is handled by creating a symbolic link from chunk 0 to original fastq
            
            cd $AlignOutputDir
            LeftReadsChunkNamePrefix=leftreads_chunk
            RightReadsChunkNamePrefix=rightreads_chunk
            if [ $chunkfastq == "YES" ]
            then

               ## splitting files into chunks before aligning;
               ## remember that one fastq read is made up of four lines
               NumChunks=`expr $nodes "-" 1`
               if [ $NumChunks -lt 1 ]
               then
	   	NumChunks=$nodes
               fi
               #if [ $NumChunks -lt 1 ]
               #then
               #    NumChunks=1
               #fi

               NumLinesInLeftFastq=`wc -l $LeftReadsFastq | cut -d ' ' -f 1`
               NumReadsInOriginalFastq=`expr $NumLinesInLeftFastq "/" 4`
               NumReadsPerChunk=`expr $NumReadsInOriginalFastq "/" $NumChunks`
               RemainderReads=`expr $NumReadsInOriginalFastq "%" $NumChunks`
               NumLinesPerChunk=`expr $NumReadsPerChunk "*" 4`
               if [ $RemainderReads -eq 0  ]
               then
	   	set +x; echo -e "\n # mod is 0; no reads for last chunk file, one idle node \n" >&2; set -x;
	   	let NumChunks-=1
               fi

               set +x; echo -e "\n # splitting read file1=$LeftReadsFastq \n" >&2; set -x;
               `split -l $NumLinesPerChunk -a 2 -d $LeftReadsFastq $LeftReadsChunkNamePrefix`

               exitcode=$?
               if [ $exitcode -ne 0 ]
               then
                      MSG="splitting of read file $LeftReadsFastq failed. exitcode=$exitcode"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                      exit $exitcode;
               fi
               if [ $paired -eq 1 ]
               then
                   set +x; echo -e "\n # splitting read file2=$RightReadsFastq \n" >&2; set -x;
 	     	   `split -l $NumLinesPerChunk -a 2 -d $RightReadsFastq $RightReadsChunkNamePrefix`
	   	   exitcode=$?
	   	   if [ $exitcode -ne 0 ]
                   then
                       MSG="splitting of read file $RightReadsFastq failed.  exitcode=$exitcode"
                       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                       exit $exitcode;
	   	   fi
               fi
            else
               NumChunks=0
               set +x; echo -e "\n\n" >&2; 
               echo "# copying original fastq into a single chunk takes too long for whole genome" >&2;
               echo "# instead, we will create symbolic link to the original file" >&2;
               echo "# adding double-0 to the chunk name, because expecting to chunk files into tens of pieces" >&2;
               echo "# need consistency in naming: two digits for chunk number" >&2;
               echo -e "\n\n" >&2; set -x; 
               #   `cp $LeftReadsFastq ${LeftReadsChunkNamePrefix}0`
               ln -s $LeftReadsFastq leftreads_chunk00

               if [ $paired -eq 1 ]
               then
               #   `cp $RightReadsFastq ${RightReadsChunkNamePrefix}0`
                  ln -s $RightReadsFastq rightreads_chunk00
               fi
            fi

            ## done chunking input fastq

            set +x; echo -e "\n\n\n" >&2;
            echo -e "#################################### DONE WITH CHUNKING DATA                 ########################################" >&2
	    echo -e "#################################### CREATING ALL QSUBS FOR ALIGMENT         ########################################" >&2
	    echo -e "#################################### SELECTING CASE --BASED ON ALIGNER       ########################################" >&2
            echo -e "\n\n" >&2; set -x

           
############################################################################################################
#####################################                               ########################################
##################################### ALIGNMENT: LOOP2 OVER CHUNKS   ########################################
#####################################                               ########################################
############################################################################################################

            allfiles=""
            for i in $(seq 0 $NumChunks)
            do
                set +x; echo -e "\n # step 1: grab  chunk $i " >&2; set -x
		echo `date`
                if (( $i < 10 ))
                then
                   Rone=${LeftReadsChunkNamePrefix}0$i
                   OutputFileSuffix=0${i}
                else
                   Rone=${LeftReadsChunkNamePrefix}$i
                   OutputFileSuffix=${i}
                fi
                if [ ! -s $Rone ]
                then
                   MSG="chunk $i of read file $LeftReadsFastq file not found"
                   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                   exit 1;
                fi
		if [ $paired -eq 1 ]
		then
                    if (( $i < 10 ))
                    then
                       Rtwo=${RightReadsChunkNamePrefix}0$i
                    else
                       Rtwo=${RightReadsChunkNamePrefix}$i
                    fi
                    if [ ! -s $Rtwo ]
                    then
			MSG="chunk $i of read file $RightReadsFastq file not found"
                        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                        #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
                fi
                
                set +x; echo -e "\n # step 2: generate the alignment cmd depending on the aligner " >&2; set -x
                
                if [ $aligner == "NOVOALIGN"  ]
		then
                    set +x; echo -e "\n###############   novoalign is used as aligner. input file in fastq format ################\n" >&2; set -x

                    if [ $paired -eq 1 ]
                    then
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $scriptdir $runfile $paired $AlignOutputDir/$Rone $AlignOutputDir/$Rtwo $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.novoaln.$SampleName.node$OutputFileSuffix" >> $AlignOutputDir/logs/novosplit.${SampleName}.node$OutputFileSuffix
                    else
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $scriptdir $runfile $paired $AlignOutputDir/$Rone $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.novoaln.$SampleName.node$OutputFileSuffix" >> $AlignOutputDir/logs/novosplit.${SampleName}.node$OutputFileSuffix
                    fi
		elif [ $aligner == "BWA_ALN" ] 
                then
                    set +x; echo -e "\n############# bwa is used as aligner. input file format is in fastq. this segment needs to be rewritten ###############\n" >&2; set -x
                    qsub_bwar1=$AlignOutputLogs/qsub.bwar1.$SampleName.node$OutputFileSuffix
                    echo "#PBS -V" > $qsub_bwar1
                    echo "#PBS -N ${pipeid}_bwar1_${SampleName}_$OutputFileSuffix" >> $qsub_bwar1
		    echo "#PBS -o $AlignOutputLogs/log.bwar1.$SampleName.node$OutputFileSuffix.ou" >> $qsub_bwar1
		    echo "#PBS -e $AlignOutputLogs/log.bwar1.$SampleName.node$OutputFileSuffix.in" >> $qsub_bwar1
                    echo "#PBS -A $pbsprj" >> $qsub_bwar1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub_bwar1
		    echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_bwar1
                    echo $ccmgres_string >> $qsub_bwar1
                    echo "#PBS -q $pbsqueue" >> $qsub_bwar1
                    echo "#PBS -m ae" >> $qsub_bwar1
                    echo "#PBS -M $email" >> $qsub_bwar1
		    echo "aprun -n 1 -d $thr $scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.R1.sai $AlignOutputDir/$Rone $scriptdir $AlignOutputLogs/log.bwar1.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwar1.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwar1.$SampleName.node$OutputFileSuffix" >> $qsub_bwar1

                    `chmod a+r $qsub_bwar1`
                    jobr1=`qsub $qsub_bwar1`
                    `qhold -h u $jobr1`
                    echo $jobr1 >> $AlignOutputLogs/ALIGNED_$SampleName
                    if [ $paired -eq 1 ]
                    then
                        set +x; echo -e "\n########## bwa aligner. paired-end reads #################\n" >&2; set -x
			qsub_bwar2=$AlignOutputLogs/qsub.bwar2.$SampleName.node$OutputFileSuffix
			echo "#PBS -V" > $qsub_bwar2
			echo "#PBS -N ${pipeid}_bwar2_${SampleName}_$OutputFileSuffix" >> $qsub_bwar2
			echo "#PBS -o $AlignOutputLogs/log.bwar2.$SampleName.node$OutputFileSuffix.ou" >> $qsub_bwar2
			echo "#PBS -e $AlignOutputLogs/log.bwar2.$SampleName.node$OutputFileSuffix.in" >> $qsub_bwar2
			echo "#PBS -A $pbsprj" >> $qsub_bwar2
			echo "#PBS -l walltime=$pbscpu" >> $qsub_bwar2
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_bwar2
                        echo $ccmgres_string >> $qsub_bwar2
			echo "#PBS -q $pbsqueue" >> $qsub_bwar2
			echo "#PBS -m ae" >> $qsub_bwar2
			echo "#PBS -M $email" >> $qsub_bwar2
			echo "$run_string $profiler_string $scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.R2.sai $AlignOutputDir/$Rtwo $scriptdir $AlignOutputLogs/log.bwar2.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwar2.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwar2.$SampleName.node$OutputFileSuffix" >> $qsub_bwar2
			`chmod a+r $qsub_bwar2`
                        jobr2=`qsub $qsub_bwar2`
			`qhold -h u $jobr2`
			echo $jobr2 >> $AlignOutputLogs/ALIGNED_$SampleName

			qsub_bwasampe=$AlignOutputLogs/qsub.bwasampe.$SampleName.node$OutputFileSuffix
			echo "#PBS -V" > $qsub_bwasampe
			echo "#PBS -N ${pipeid}_bwasampe_${SampleName}_$OutputFileSuffix" >> $qsub_bwasampe
			echo "#PBS -o $AlignOutputLogs/log.bwasampe.$SampleName.node$OutputFileSuffix.ou" >> $qsub_bwasampe
			echo "#PBS -e $AlignOutputLogs/log.bwasampe.$SampleName.node$OutputFileSuffix.in" >> $qsub_bwasampe
			echo "#PBS -A $pbsprj" >> $qsub_bwasampe
			echo "#PBS -l walltime=$pbscpu" >> $qsub_bwasampe
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_bwasampe
                        echo $ccmgres_string >> $qsub_bwasampe
			echo "#PBS -q $pbsqueue" >> $qsub_bwasampe
			echo "#PBS -m ae" >> $qsub_bwasampe
			echo "#PBS -M $email" >> $qsub_bwasampe
			echo "#PBS -W depend=afterok:$jobr1:$jobr2" >> $qsub_bwasampe
			echo "$run_string $profiler_string $scriptdir/bwaS2.sh $alignerdir $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.R1.sai $outputsamfileprefix.node$OutputFileSuffix.R2.sai $AlignOutputDir/$Rone $AlignOutputDir/$Rtwo $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $samdir $AlignOutputLogs/log.bwasampe.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwasampe.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwasampe.$SampleName.node$OutputFileSuffix" >> $qsub_bwasampe
			`chmod a+r $qsub_bwasampe`
                        jobwa=`qsub $qsub_bwasampe`
			`qhold -h u $jobwa`
			echo $jobwa >> $AlignOutputLogs/ALIGNED_$SampleName
                    else
                        set +x; echo -e "\n############# bwa aligner. single read #################\n" >&2; set -x
			qsub_bwasamse=$AlignOutputLogs/qsub.bwasamse.$SampleName.node$OutputFileSuffix
			echo "#PBS -V" > $qsub_bwasamse
			echo "#PBS -N ${pipeid}_bwasamse_${SampleName}_$OutputFileSuffix" >> $qsub_bwasamse
			echo "#PBS -o $AlignOutputLogs/log.bwasamse.$SampleName.node$OutputFileSuffix.ou" >> $qsub_bwasamse
			echo "#PBS -e $AlignOutputLogs/log.bwasamse.$SampleName.node$OutputFileSuffix.in" >> $qsub_bwasamse
			echo "#PBS -A $pbsprj" >> $qsub_bwasamse
			echo "#PBS -l walltime=$pbscpu" >> $qsub_bwasamse
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_bwasamse
                        echo $ccmgres_string >> $qsub_bwasamse
			echo "#PBS -q $pbsqueue" >> $qsub_bwasamse
			echo "#PBS -m ae" >> $qsub_bwasamse
			echo "#PBS -M $email" >> $qsub_bwasamse
			echo "#PBS -W depend=afterok:$jobr1" >> $qsub_bwasamse
			echo "aprun -n 1 -d $thr $scriptdir/bwaS3.sh $alignerdir $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.R1.sai $AlignOutputDir/$Rone $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $samdir $AlignOutputLogs/log.bwasamse.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwasamse.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwasamse.$SampleName.node$OutputFileSuffix" >> $qsub_bwasamse
			`chmod a+r $qsub_bwasamse`
                        jobwa=`qsub $qsub_bwasamse`
			`qhold -h u $jobwa`
                        echo $qsub_bwasamse >> $AlignOutputLogs/ALIGNED_$SampleName
                    fi
                elif [ $aligner == "BWA_MEM" ]
                then
                    set +x; echo -e "\n################ bwa mem is used as aligner. input file format is in fastq #################\n" >&2; set -x
                    if [ $paired -eq 1 ]
                    then
                        set +x; echo -e "\n############### bwa mem aligner. paired-end reads ###################\n" >&2; set -x
                        jobfile=$AlignOutputDir/logs/bwamem.$SampleName.node$OutputFileSuffix.jobfile

                        if [ $chunkfastq == "YES" ]
                        then
			   ### MUST FIX QSUB VARIABLE NAME PROPERLY
                           echo "nohup $profiler_string $scriptdir/bwamem_pe.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $AlignOutputDir/$Rone $AlignOutputDir/$Rtwo $scriptdir $samdir $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix > $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in" >> $qsub3
                        elif [ $chunkfastq == "NO" ]
                        then
                           echo "nohup $scriptdir/bwamem_pe_markduplicates.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix $AlignOutputDir/$Rone $AlignOutputDir/$Rtwo $runfile $AlignOutputDir/logs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputDir/logs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $jobfile $RGparms $AlignOutputLogs > $AlignOutputDir/logs/log.bwamem.$SampleName.node$OutputFileSuffix.in" > $jobfile
                           jobfilename=$( basename $jobfile )
                           echo "$AlignOutputDir/logs $jobfilename" >> $AlignOutputLogs/AlignAnisimov.joblist
                        fi
                        `chmod ug=rwx $jobfile`

                     else
                        set +x; echo -e "\n################ bwa mem aligner. single-end reads ################\n" >&2; set -x
                        qsub_bwamem=$AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix
                        echo "#PBS -V" > $qsub_bwamem
                        echo "#PBS -N ${pipeid}_bwamem_${SampleName}_$OutputFileSuffix" >> $qsub_bwamem
                        echo "#PBS -o $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou" >> $qsub_bwamem
                        echo "#PBS -e $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in" >> $qsub_bwamem
                        echo "#PBS -A $pbsprj" >> $qsub_bwamem
                        echo "#PBS -l walltime=$pbscpu" >> $qsub_bwamem
                        echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_bwamem
                        echo "#PBS -q $pbsqueue" >> $qsub_bwamem
                        echo "#PBS -m ae" >> $qsub_bwamem
                        echo "#PBS -M $email" >> $qsub_bwamem

                        if [ $chunkfastq == "YES" ]
                        then
                           echo "aprun -n 1 -d $thr $scriptdir/bwamem_se.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $AlignOutputDir/$Rone $scriptdir $samdir $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix" >> $qsub_bwamem
                        elif [ $chunkfastq == "NO" ]
                        then
                           echo "aprun -n 1 -d $thr $scriptdir/bwamem_se_markduplicates.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir $outputsamfileprefix.node$OutputFileSuffix.bam $AlignOutputDir/$Rone $scriptdir $samdir $samblasterdir $picardir $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix $RGparms" >> $qsub_bwamem
                        fi
   

                       `chmod a+r $qsub_bwamem`
                       jobbwamemse=`qsub $qsub_bwamem`
                       `qhold -h u $jobbwamemse`
                       echo $jobbwamemse >> $AlignOutputLogs/ALIGNED_$SampleName


                    fi
                fi

                if (( $i < 10 ))
                then
                   allfiles=$allfiles" $outputsamfileprefix.node0$i.bam" # this will be used later for merging
                else
                   allfiles=$allfiles" $outputsamfileprefix.node$i.bam" # this will be used later for merging
                fi
		echo `date`
            done # end loop over chunks of the current fastq

######################################################################
############# end loop over chunks of the current fastq ##############
######################################################################
######################################################################

######################################################################
############# WE ARE STILL INSIDE THE LOOP OVER INPUT FASTQ!!! #######
######################################################################


set +x; echo -e "\n\n\n" >&2;
echo -e "##########################################################################################################################################" >&2
echo -e "###########################                                                                            ###################################" >&2
echo -e "###########################   FORM POST-ALIGNMENT QSUBS: MERGINE, SORTING, MARKING DUPLICATES          ###################################" >&2
echo -e "###########################   SKIP THIS BLOCK IF READS WHERE NOT CHUNKED                               ###################################" >&2
echo -e "###########################                                                                            ###################################" >&2
echo -e "##########################################################################################################################################" >&2
echo -e "\n\n\n" >&2; set -x

            cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/pbs.ALIGNED

            if [ $chunkfastq == "YES" ]
            then
   	       #ALIGNED=$( cat $AlignOutputLogs/ALIGNED_* | sed "s/\.[a-z]*//" | tr "\n" ":" )
	       #ALIGNED=$( cat $AlignOutputLogs/ALIGNED_* | sed "s/\..*//" | tr "\n" ":" )


	       listfiles=$( echo $allfiles  | tr " " ":" | sed "s/::/:/g" )
               if [ $sortool == "NOVOSORT" ]
               then
                   set +x; echo -e "\n # merging aligned chunks with novosort \n" >&2; set -x

		   echo "$scriptdir/mergenovo.sh $AlignOutputDir $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $AlignOutputLogs/log.sortmerge_novosort.$SampleName.in $AlignOutputLogs/log.sortmerge_novosort.$SampleName.ou $email $AlignOutputLogs/qsub.sortmerge.novosort.$SampleName" >> $AlignOutputDir/logs/mergenovo.${SampleName}
		   #`chmod a+r $qsub_sortmerge`
		   #mergejob=`qsub $qsub_sortmerge`
		   #`qhold -h u $mergejob`
		   #echo $mergejob  > $AlignOutputLogs/MERGED_$SampleName
               else
                   set +x; echo -e "\n # merging aligned chunks with picard \n" >&2; set -x 
		   qsub_sortmerge=$AlignOutputLogs/qsub.sortmerge.picard.$SampleName
		   echo "#PBS -V" > $qsub_sortmerge
		   echo "#PBS -A $pbsprj" >> $qsub_sortmerge
		   echo "#PBS -N ${pipeid}_sortmerge_picard_$SampleName" >> $qsub_sortmerge
		   echo "#PBS -l walltime=$pbscpu" >> $qsub_sortmerge
		   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_sortmerge
		   echo "#PBS -o $AlignOutputLogs/log.sortmerge.picard.$SampleName.ou" >> $qsub_sortmerge
		   echo "#PBS -e $AlignOutputLogs/log.sortmerge.picard.$SampleName.in" >> $qsub_sortmerge
		   echo "#PBS -q $pbsqueue" >> $qsub_sortmerge
		   echo "#PBS -m ae" >> $qsub_sortmerge
		   echo "#PBS -M $email" >> $qsub_sortmerge
		   echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub_sortmerge
		   echo "aprun -n 1 -d $thr $scriptdir/mergepicard.sh $AlignOutputDir $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $AlignOutputLogs/log.sortmerge.picard.$SampleName.in $AlignOutputLogs/log.sortmerge.picard.$SampleName.ou $email $AlignOutputLogs/qsub.sortmerge.picard.$SampleName" >> $qsub_sortmerge
		   `chmod a+r $qsub_sortmerge`
		   mergejob=`qsub $qsub_sortmerge`
		   `qhold -h u $mergejob`
		   echo $mergejob  > $AlignOutputLogs/MERGED_$SampleName
               fi

            fi
          fi
          (( inputfastqcounter++ )) # was not initialized at beginning of loop, so starts at zero
	done <  $TheInputFile # end loop over input fastq

######################################################################
#############         end loop over iinput fastq        ##############
######################################################################





set +x; echo -e "\n\n\n" >&2
echo "#####################################################################################################################" >&2
echo "#####################################                                        ########################################" >&2
echo "#####################################  SCHEDULE BWA-MEM QSUBS CREATED ABOVE  ########################################" >&2
echo "#####################################                                        ########################################" >&2
echo "#####################################################################################################################" >&2
echo -e "\n\n" >&2; set -x



        if [ $chunkfastq == "NO" -a $aligner == "BWA_MEM" ]
        then
 
           set +x; echo -e "\n ### update autodocumentation script ### \n"; set -x;
           echo -e "# @begin $aligner" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "   # @in datafilestoalign @as fastqc_inputs" >> $outputdir/WorkflowAutodocumentationScript.sh
           AlignedFastqPathTemplate="SampleName/align/SampleName.node00.wdups.sorted.bam"
           echo -e "   # @out input_to_realrecal @as fastq_aligned_sorted_markdup @URI ${AlignedFastqPathTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "# @end $aligner" >> $outputdir/WorkflowAutodocumentationScript.sh


           case $run_method in
           "LAUNCHER")
              # increment so number of nodes = number fo input fastq + 1, even when there is only one input fastq
              # otherwise nowhere for launcher to run, due to the -N 1 option in aprun
              numalignnodes=$(( inputfastqcounter + 1))
   
              set +x; echo -e "\n # run_method is LAUNCHER. scheduling the Anisimov Launcher\n" >&2; set -x
              
              qsubAlignLauncher=$AlignOutputLogs/qsub.align.Anisimov
              echo "#!/bin/bash" > $qsubAlignLauncher
              echo "#PBS -V" >> $qsubAlignLauncher
              echo "#PBS -A $pbsprj" >> $qsubAlignLauncher
              echo "#PBS -N ${pipeid}_align_Anisimov" >> $qsubAlignLauncher
              echo "#PBS -l walltime=$pbscpu" >> $qsubAlignLauncher
              echo "#PBS -l nodes=$numalignnodes:ppn=$thr" >> $qsubAlignLauncher
              echo "#PBS -o $AlignOutputLogs/log.align.Anisimov.ou" >> $qsubAlignLauncher
              echo "#PBS -e $AlignOutputLogs/log.align.Anisimov.in" >> $qsubAlignLauncher
              echo "#PBS -q $pbsqueue" >> $qsubAlignLauncher
              echo "#PBS -m ae" >> $qsubAlignLauncher
              echo "#PBS -M $email" >> $qsubAlignLauncher
              echo "aprun -n $numalignnodes -N 1 -d $thr ~anisimov/scheduler/scheduler.x $AlignOutputLogs/AlignAnisimov.joblist /bin/bash > $AlignOutputLogs/AlignAnisimov.joblist.log" >> $qsubAlignLauncher

              echo "exitcode=\$?" >> $qsubAlignLauncher
              echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubAlignLauncher
              echo "   echo -e \"\n\n AlignAnisimov failed with exit code = \$exitcode \n logfile=$AlignOutputLogs/log.AlignAnisimov.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubAlignLauncher
              echo "   exit 1" >> $qsubAlignLauncher
              echo "fi" >> $qsubAlignLauncher

              AlignAnisimovJoblistId=`qsub $qsubAlignLauncher`
              echo $AlignAnisimovJoblistId >> $TopOutputLogs/pbs.ALIGNED # so that this job could be released in the next section. Should it be held to begin with?
              echo $AlignAnisimovJoblistId > $TopOutputLogs/pbs.MERGED # so that summaryok and start_realrecal_block.sh could depend on this job, in case when there is no merging: a sigle chunk
           ;;           
           "APRUN")
              while read SampleName
              do              
                 set +x; echo -e "\n # run_method is APRUN. scheduling qsubs\n" >&2; set -x
            	 if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
                 then
                     SampleName=$( echo $SampleLine )
	         else
                     SampleName=$( echo -e "$SampleLine" | cut -f 2 )
                 fi
                 AlignOutputDir=$outputdir/${SampleName}/align
                 
                 qsub_bwa=$AlignOutputDir/logs/qsub.bwa.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_bwa


                 ###############################
                 echo -e "\n################# constructing qsub for bwamem\n"
                 echo "#PBS -N ${pipeid}_bwa.${SampleName}" >> $qsub_bwa
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_bwa
                 echo "#PBS -o $AlignOutputDir/logs/log.bwa.${SampleName}.ou" >> $qsub_bwa
                 echo "#PBS -e $AlignOutputDir/logs/log.bwa.${SampleName}.in" >> $qsub_bwa
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_bwa

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
#                 sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/logs/bwamem.${SampleName}.node$OutputFileSuffix.jobfile >> $qsub_bwa
                 echo "aprun -n 1 -N 1 -d $thr /bin/bash $AlignOutputDir/logs/bwamem.${SampleName}.node$OutputFileSuffix.jobfile" >> $qsub_bwa

                 bwa_job=`qsub $qsub_bwa`
                 `qhold -h u $bwa_job`
                 echo $bwa_job >> $TopOutputLogs/pbs.ALIGNED # so that this job could be released in the next section. Should it be held to begin with?
                 echo $bwa_job >> $TopOutputLogs/pbs.MERGED # so that summaryok and start_realrecal_block.sh could depend on this job, in case when there is no merging: a sigle chunk
              done < $TheInputFile # done looping over samples
           ;;
           "QSUB")
              while read SampleName
              do
                 set +x; echo -e "\n # run_method is QSUB. scheduling qsubs\n" >&2; set -x
            	 if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
                 then
                     SampleName=$( echo $SampleLine )
	         else
                     SampleName=$( echo -e "$SampleLine" | cut -f 2 )
                 fi
                 AlignOutputDir=$outputdir/${SampleName}/align
                                  
                 qsub_bwa=$AlignOutputDir/logs/qsub.bwa.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_bwa


                 ###############################
                 echo -e "\n################# constructing qsub for bwamem\n"
                 echo "#PBS -N ${pipeid}_bwa.${SampleName}" >> $qsub_bwa
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_bwa
                 echo "#PBS -o $AlignOutputDir/logs/log.bwa.${SampleName}.ou" >> $qsub_bwa
                 echo "#PBS -e $AlignOutputDir/logs/log.bwa.${SampleName}.in" >> $qsub_bwa
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_bwa

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 echo "/bin/bash $AlignOutputDir/logs/bwamem.${SampleName}.node$OutputFileSuffix.jobfile" >> $qsub_bwa

                 bwa_job=`qsub $qsub_bwa`
                 `qhold -h u $bwa_job`
                 echo $bwa_job >> $TopOutputLogs/pbs.ALIGNED # so that this job could be released in the next section. Should it be held to begin with?
                 echo $bwa_job >> $TopOutputLogs/pbs.MERGED # so that summaryok and start_realrecal_block.sh could depend on this job, in case when there is no merging: a sigle chunk
              done < $TheInputFile # done looping over samples
           ;;
           esac

        fi




set +x; echo -e "\n\n\n" >&2
echo "#####################################################################################################################################" >&2;
echo "#####################################                                                        ########################################" >&2;
echo "#####################################  SCHEDULE NOVOALIGN and MERGENOVO QSUBS CREATED ABOVE  ########################################" >&2;
echo "#####################################                                                        ########################################" >&2;
echo "#####################################################################################################################################" >&2;
echo -e "\n\n" >&2; set -x


        if [ $chunkfastq == "YES" -a $aligner == "NOVOALIGN" -a $sortool == "NOVOSORT" ]
        then

           set +x; echo -e "\n ### update autodocumentation script ### \n"; set -x;
           echo -e "# @begin chunk_fastq" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "   # @in datafilestoalign @as fastqc_inputs" >> $outputdir/WorkflowAutodocumentationScript.sh
           InputFastqPathTemplate="SampleName/align/{left/right}reads_chunk{number}"
           echo -e "   # @out chunks @as chunked_fastq @URI ${InputFastqPathTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "# @end $chunk_fastq" >> $outputdir/WorkflowAutodocumentationScript.sh

           echo -e "# @begin ${aligner}" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "   # @in chunks @as chunked_fastq " >> $outputdir/WorkflowAutodocumentationScript.sh
           AlignedFastqPathTemplate="/SampleName/align/{left/right}reads.node{number}.bam"
           echo -e "   # @out alignment_output @as aligned_fastq @URI ${AlignedFastqPathTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "# @end ${aligner}" >> $outputdir/WorkflowAutodocumentationScript.sh

           echo -e "# @begin merge_${sortool}_${markduplicatestool}" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "   # @in alignment_output @as aligned_fastq" >> $outputdir/WorkflowAutodocumentationScript.sh
           MergedFastqPathTemplate="SampleName/align/SampleName.wdups.sorted.bam"
           echo -e "   # @out input_to_realrecal @as fastq_aligned_sorted_markdup @URI ${MergedFastqPathTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
           echo -e "# @end merge_${sortool}_${markduplicatestool}" >> $outputdir/WorkflowAutodocumentationScript.sh

           case $run_method in
           "APRUN")
              while read SampleName
              do
            	 if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
                 then
                     SampleName=$( echo $SampleLine )
	         else
                     SampleName=$( echo -e "$SampleLine" | cut -f 2 )
                 fi
                 AlignOutputDir=$outputdir/${SampleName}/align

                 # ADD LOOP OVER CHUNKS
                 for i in $(seq 0 $NumChunks)
                 do

                    if (( $i < 10 ))
                    then
                       OutputFileSuffix=0${i}
                    else
                       OutputFileSuffix=${i}
                    fi

                    set +x; echo -e "\n # scheduling qsubs\n" >&2; set -x
                    qsub_novosplit=$AlignOutputDir/logs/qsub.novosplit.${SampleName}.node${OutputFileSuffix}
                    # appending the generic header to the qsub
                    cat $outputdir/qsubGenericHeader > $qsub_novosplit


                    ###############################
                    set +x; echo -e "\n # constructing qsub for novosplit\n" >&2; set -x
                    echo "#PBS -N ${pipeid}_novoalign.${SampleName}.node${OutputFileSuffix}" >> $qsub_novosplit
                    echo "#PBS -l walltime=$pbscpu" >> $qsub_novosplit
                    echo "#PBS -o $AlignOutputDir/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.ou" >> $qsub_novosplit
                    echo "#PBS -e $AlignOutputDir/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.in" >> $qsub_novosplit
                    echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_novosplit

                    # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                    sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/logs/novosplit.${SampleName}.node$OutputFileSuffix >> $qsub_novosplit 

                    novosplit_job=`qsub $qsub_novosplit`
                    `qhold -h u $novosplit_job`
                    echo $novosplit_job >> $TopOutputLogs/pbs.ALIGNED_${SampleName} # so that this job could be released in the next section. Should it be held to begin with?

                 done # done looping over chunks of a sample
                 cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/pbs.ALIGNED
                 alignids=$( cat $TopOutputLogs/ALIGNED_$SampleName | sed "s/\..*//" | tr "\n" ":" )

                 qsub_mergenovo=$AlignOutputDir/logs/qsub.mergenovo.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_mergenovo




                 ###############################
                 set +x; echo -e "\n # constructing qsub for mergenovo\n" >&2; set -x
                 echo "#PBS -N ${pipeid}_mergenovo" >> $qsub_mergenovo
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_mergenovo
                 echo "#PBS -o $AlignOutputLogs/log.mergenovo.ou" >> $qsub_mergenovo
                 echo "#PBS -e $AlignOutputLogs/log.mergenovo.in" >> $qsub_mergenovo
                 # add dependency on novosplit job
                 echo -e "#PBS -W depend=afterok:$alignids" >> $qsub_mergenovo
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_mergenovo

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/logs/mergenovo.${SampleName} >> $qsub_mergenovo 

                 mergenovo_job=`qsub $qsub_mergenovo`
                 `qhold -h u $mergenovo_job`
                 echo $mergenovo_job > $TopOutputLogs/pbs.MERGED # so that this job could be released in the next section, and start_realrecal_block.sh could depend on it

                 # release all held novosplit jobs for this sample
                 `qrls -h u $alignids`

              done < $TheInputFile # done looping over samples
           ;;
           "QSUB")
              while read SampleName
              do
            	 if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
                 then
                     SampleName=$( echo $SampleLine )
	         else
                     SampleName=$( echo -e "$SampleLine" | cut -f 2 )
                 fi
                 AlignOutputDir=$outputdir/${SampleName}/align


                 for i in $(seq 0 $NumChunks)
                 do

                    if (( $i < 10 ))
                    then
                       OutputFileSuffix=0${i}
                    else
                       OutputFileSuffix=${i}
                    fi

                    set +x; echo -e "\n # scheduling qsubs\n" >&2; set -x
                    qsub_novosplit=$AlignOutputDir/logs/qsub.novosplit.${SampleName}.node${OutputFileSuffix}

                    # appending the generic header to the qsub
                    cat $outputdir/qsubGenericHeader > $qsub_novosplit


                    ###############################
                    set +x; echo -e "\n # constructing qsub for novosplit\n" >&2; set -x
                    echo "#PBS -N ${pipeid}_novoalign.${SampleName}.node${OutputFileSuffix}" >> $qsub_novosplit
                    echo "#PBS -l walltime=$pbscpu" >> $qsub_novosplit
                    echo "#PBS -o $AlignOutputDir/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.ou" >> $qsub_novosplit
                    echo "#PBS -e $AlignOutputDir/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.in" >> $qsub_novosplit
                    echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_novosplit
                    cat $AlignOutputDir/logs/novosplit.${SampleName}.node${OutputFileSuffix} >> $qsub_novosplit
                    novosplit_job=`qsub $qsub_novosplit`
                    `qhold -h u $novosplit_job`
                    echo $novosplit_job >> $TopOutputLogs/pbs.ALIGNED # so that this job could be released in the next section. Should it be held to begin with?

                 done # done looping over chunks of a sample
                 cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/pbs.ALIGNED
                 alignids=$( cat $TopOutputLogs/ALIGNED_$SampleName | sed "s/\..*//" | tr "\n" ":" )

                 qsub_mergenovo=$AlignOutputDir/logs/qsub.mergenovo.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_mergenovo


                 ###############################
                 set +x; echo -e "\n # constructing qsub for mergenovo\n" >&2; set -x
                 echo "#PBS -N ${pipeid}_mergenovo" >> $qsub_mergenovo
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_mergenovo
                 echo "#PBS -o $AlignOutputLogs/log.mergenovo.ou" >> $qsub_mergenovo
                 echo "#PBS -e $AlignOutputLogs/log.mergenovo.in" >> $qsub_mergenovo
                 # add dependency on novosplit job
                 echo -e "#PBS -W depend=afterok:$alignids" >> $qsub_mergenovo
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_mergenovo

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/logs/mergenovo.${SampleName} >> $qsub_mergenovo

                 mergenovo_job=`qsub $qsub_mergenovo`
                 `qhold -h u $mergenovo_job`
                 echo $mergenovo_job > $TopOutputLogs/pbs.MERGED # so that this job could be released in the next section, and start_realrecal_block.sh could depend on it

                 # release all held novosplit jobs for this sample
                 `qrls -h u $alignids`


              done < $TheInputFile # done looping over samples
           ;;
           "LAUNCHER")
              # clear out the joblists
              truncate -s 0 $AlignOutputLogs/novosplit.AnisimovJoblist
              truncate -s 0 $AlignOutputLogs/mergenovo.AnisimovJoblist

              while read SampleName
              do
            	 if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
                 then
                     SampleName=$( echo $SampleLine )
	         else
                     SampleName=$( echo -e "$SampleLine" | cut -f 2 )
                 fi
                 AlignOutputDir=$outputdir/${SampleName}/align
                 
                 for i in $(seq 0 $NumChunks)
                 do

                    if (( $i < 10 ))
                    then
                       OutputFileSuffix=0${i}
                    else
                       OutputFileSuffix=${i}
                    fi

                    # creating a qsub out of the job file
                    # need to prepend "nohup" and append log file name, so that logs are properly created when Anisimov launches these jobs
                    novosplit_log=${AlignOutputDir}/logs/log.novosplit.${SampleName}.node$OutputFileSuffix.in
                    awk -v awkvar_novosplitlog=$novosplit_log '{print "nohup "$0" > "awkvar_novosplitlog}' $AlignOutputDir/logs/novosplit.${SampleName}.node${OutputFileSuffix} > $AlignOutputDir/logs/jobfile.novosplit.${SampleName}.node${OutputFileSuffix}
                    echo "$AlignOutputDir/logs/ jobfile.novosplit.${SampleName}.node${OutputFileSuffix}" >> $AlignOutputLogs/novosplit.AnisimovJoblist
                 done # done going through chunks

                 mergenovo_log=${AlignOutputDir}/logs/log.mergenovo.${SampleName}.in
                 awk -v awkvar_mergenovolog=$mergenovo_log '{print "nohup "$0" > "awkvar_mergenovolog}' $AlignOutputDir/logs/mergenovo.${SampleName} > $AlignOutputDir/logs/jobfile.mergenovo.${SampleName}
                 echo "$AlignOutputDir/logs/ jobfile.mergenovo.${SampleName}" >> $AlignOutputLogs/mergenovo.AnisimovJoblist
              done < $TheInputFile # done looping over samples


              set +x; echo -e "\n\n ########### scheduling Anisimov Launcher joblists #############\n\n" >&2; set -x
              qsub_novosplit_anisimov=$AlignOutputLogs/qsub.novosplit.AnisimovLauncher
              qsub_mergenovo_anisimov=$AlignOutputLogs/qsub.mergenovo.AnisimovLauncher

              # appending the generic header to the qsub
              cat $outputdir/qsubGenericHeader > $qsub_novosplit_anisimov
              cat $outputdir/qsubGenericHeader > $qsub_mergenovo_anisimov


              ###############################
              set +x; echo -e "\n ################# constructing qsub for novosplit\n" >&2; set -x
              echo "#PBS -N ${pipeid}_novoalign_Anisimov" >> $qsub_novosplit_anisimov
              echo "#PBS -l walltime=$pbscpu" >> $qsub_novosplit_anisimov
              echo "#PBS -o $AlignOutputLogs/log.novosplit.Anisimov.ou" >> $qsub_novosplit_anisimov
              echo "#PBS -e $AlignOutputLogs/log.novosplit.Anisimov.in" >> $qsub_novosplit_anisimov

              # number of nodes required for alignment is equal to the number of samples times number of chunks into which every sample is broken
              # counter started at zero, so in the end reflects the true number of fsatq samples
              # for some stupid reason NumChunks is actually from 0 to number_of_chunks-1, so we have to increment it now
              (( NumChunks++ ))
              NumberOfNodes=$(( inputfastqcounter * NumChunks )) 
              (( NumberOfNodes++ )) # plus one for the launcher
              echo -e "#PBS -l nodes=$NumberOfNodes:ppn=$thr\n" >> $qsub_novosplit_anisimov
              echo "aprun -n $NumberOfNodes -N 1 -d $thr ~anisimov/scheduler/scheduler.x $AlignOutputLogs/novosplit.AnisimovJoblist /bin/bash > $AlignOutputLogs/novosplit.AnisimovLauncher.log" >> $qsub_novosplit_anisimov

              ### iForge - specific use of launcher - not sure how to handle the switch, so keeping it commented out for now
              ### echo "module load intel/12.0.4" >> $qsub_novosplit_anisimov
              ### echo "module load openmpi-1.4.3-intel-12.0.4" >> $qsub_novosplit_anisimov
              #echo "cat $PBS_NODEFILE | sort -u | awk -v n=1 '{for(i=0;i<n;i++) print \$0}' > ${AlignOutputLogs}/HOSTLIST" >> $qsub_novosplit_anisimov
              ### echo "mpiexec -n $NumberOfNodes  --pernode -machinefile \$PBS_NODEFILE --mca btl tcp,self  ${launcherdir}/scheduler.x $AlignOutputLogs/novosplit.AnisimovJoblist /bin/bash > $AlignOutputLogs/novosplit.AnisimovLauncher.log" >> $qsub_novosplit_anisimov

              novosplit_job=`qsub $qsub_novosplit_anisimov`
              `qhold -h u $novosplit_job`
              echo $novosplit_job >> $TopOutputLogs/pbs.ALIGNED # so that this job could be released in the next section. Should it be held to begin with?


              ###############################
              set +x; echo -e "\n ################# constructing qsub for mergenovo\n" >&2; set -x
              echo "#PBS -N ${pipeid}_mergenovo_Anisimov" >> $qsub_mergenovo_anisimov
              echo "#PBS -l walltime=$pbscpu" >> $qsub_mergenovo_anisimov
              echo "#PBS -o $AlignOutputLogs/log.mergenovo.Anisimov.ou" >> $qsub_mergenovo_anisimov
              echo "#PBS -e $AlignOutputLogs/log.mergenovo.Anisimov.in" >> $qsub_mergenovo_anisimov
              # add dependency on novosplit job
              echo -e "#PBS -W depend=afterok:$novosplit_job" >> $qsub_mergenovo_anisimov

              # number of nodes required for mergenovo is equal to the number of samples plus one for the launcher
              NumberOfNodes=$(( inputfastqcounter + 1 )) # counter started at zero, so in the end reflects the true number of fsatq samples
              echo -e "#PBS -l nodes=$NumberOfNodes:ppn=$thr\n" >> $qsub_mergenovo_anisimov
              echo "aprun -n $NumberOfNodes -N 1 -d $thr ~anisimov/scheduler/scheduler.x $AlignOutputLogs/mergenovo.AnisimovJoblist /bin/bash > $AlignOutputLogs/mergenovo.AnisimovLauncher.log" >> $qsub_mergenovo_anisimov

              mergenovo_job=`qsub $qsub_mergenovo_anisimov`
              `qhold -h u $mergenovo_job`
              echo $mergenovo_job > $TopOutputLogs/pbs.MERGED # so that this job could be released in the next section, and start_realrecal_block.sh could depend on it 

           ;;
           "SERVER")
              # will fill in later
           ;;
           esac
        fi







set +x; echo -e "\n\n\n" >&2
echo "#####################################################################################################################################" >&2
echo "###############################     WRAP UP ALIGNMENT BLOCK                                  ########################################" >&2
echo "###############################     ALL QSUB SCRIPTS BELOW WILL RUN AFTER ALIGNMENT IS DONE  ########################################" >&2
echo "#####################################################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

        
	pbsids=$( cat $TopOutputLogs/pbs.MERGED | sed "s/\..*//" | tr "\n" ":" )
	#extraids=$( cat $TopOutputLogs/pbs.EXTRACTREADS | sed "s/\..*//" | tr "\n" " " )
        mergeids=$( echo $pbsids | tr ":" " " )
        alignids=$( cat $TopOutputLogs/pbs.ALIGNED | sed "s/\..*//" | tr "\n" " " )

        ## generating summary redmine email if analysis ends here
	set +x; echo -e "\n # wrap up and produce summary table if analysis ends here or call realign if analysis continues \n" >&2; set -x;
	if [ $analysis == "ALIGNMENT" -o $analysis == "ALIGN" -o $analysis == "ALIGN_ONLY" ]
	then
	    set +x; echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" >&2; set -x;
            # release all held jobs
            `qrls -h u $alignids`
            `qrls -h u $mergeids`
            #`qrls -h u $extraids`
     
	    lastjobid=""
            cleanjobid=""

            if [ $cleanupflag == "YES" ]
            then 
		set +x; echo -e "\n ###### Removing temporary files  ######\n" >&2; set -x;
		qsub_cleanup=$TopOutputLogs/qsub.cleanup.align
		echo "#PBS -V" > $qsub_cleanup
		echo "#PBS -A $pbsprj" >> $qsub_cleanup
		echo "#PBS -N ${pipeid}_cleanup_aln" >> $qsub_cleanup
		echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
		echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
		echo "#PBS -o $TopOutputLogs/log.cleanup.align.ou" >> $qsub_cleanup
		echo "#PBS -e $TopOutputLogs/log.cleanup.align.in" >> $qsub_cleanup
		echo "#PBS -q $pbsqueue" >> $qsub_cleanup
		echo "#PBS -m ae" >> $qsub_cleanup
		echo "#PBS -M $email" >> $qsub_cleanup
		echo "#PBS -W depend=afterok:$pbsids" >> $qsub_cleanup
		echo "aprun -n 1 -d $thr $scriptdir/cleanup.sh $outputdir $analysis $TopOutputLogs/log.cleanup.align.in $TopOutputLogs/log.cleanup.align.ou $email $TopOutputLogs/qsub.cleanup.align"  >> $qsub_cleanup
		`chmod a+r $qsub_cleanup`
		cleanjobid=`qsub $qsub_cleanup`
		echo $cleanjobid >> $outputdir/logs/pbs.CLEANUP
            fi

            `sleep 30s`
	    set +x; echo -e "\n ###### Generating Summary report   ######\n" >&2; set -x;
	    qsub_summary=$TopOutputLogs/qsub.summary.aln.allok
	    echo "#PBS -V" > $qsub_summary
	    echo "#PBS -A $pbsprj" >> $qsub_summary
	    echo "#PBS -N ${pipeid}_summaryok" >> $qsub_summary
	    echo "#PBS -l walltime=01:00:00" >> $qsub_summary # 1 hour should be more than enough
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
	    echo "#PBS -o $TopOutputLogs/log.summary.aln.ou" >> $qsub_summary
	    echo "#PBS -e $TopOutputLogs/log.summary.aln.in" >> $qsub_summary
	    echo "#PBS -q $pbsqueue" >> $qsub_summary
	    echo "#PBS -m ae" >> $qsub_summary
	    echo "#PBS -M $email" >> $qsub_summary
            if [ `expr ${#cleanjobid}` -gt 0 ]
            then
		echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub_summary
            else
		echo "#PBS -W depend=afterok:$pbsids" >> $qsub_summary
            fi
	    echo "$scriptdir/summary.sh $runfile $email exitok $reportticket"  >> $qsub_summary
	    `chmod a+r $qsub_summary`
	    lastjobid=`qsub $qsub_summary`
	    echo $lastjobid >> $TopOutputLogs/pbs.SUMMARY

	    if [ `expr ${#lastjobid}` -lt 1 ]
	    then

		qsub_summary=$TopOutputLogs/qsub.summary.aln.afterany
		echo "#PBS -V" > $qsub_summary
		echo "#PBS -A $pbsprj" >> $qsub_summary
		echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub_summary
		echo "#PBS -l walltime=01:00:00" >> $qsub_summary # 1 hour should be more than enough
		echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
		echo "#PBS -o $TopOutputLogs/log.summary.aln.afterany.ou" >> $qsub_summary
		echo "#PBS -e $TopOutputLogs/log.summary.aln.afterany.in" >> $qsub_summary
		echo "#PBS -q $pbsqueue" >> $qsub_summary
		echo "#PBS -m ae" >> $qsub_summary
		echo "#PBS -M $email" >> $qsub_summary
		echo "#PBS -W depend=afterany:$pbsids" >> $qsub_summary
		echo "$scriptdir/summary.sh $runfile $email exitnotok $reportticket"  >> $qsub_summary
		`chmod a+r $qsub_summary`
		badjobid=`qsub $qsub_summary`
		echo $badjobid >> $TopOutputLogs/pbs.SUMMARY
	    fi
	fi

	if [ $analysis == "REALIGNMENT" -o $analysis == "REALIGN" -o $analysis == "MULTIPLEXED" ]
	then
            set +x; echo -e "\n ###### analysis continues with realignment   ###### \n" >&2; set -x;
	    qsub_realign=$TopOutputLogs/qsub.start_realrecal_block
	    echo "#PBS -V" > $qsub_realign
	    echo "#PBS -A $pbsprj" >> $qsub_realign
	    echo "#PBS -N ${pipeid}_START_REALRECAL_BLOCK" >> $qsub_realign
	    echo "#PBS -l walltime=01:00:00" >> $qsub_realign # 1 hour should be more than enough
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub_realign
	    echo "#PBS -o $TopOutputLogs/start_realrecal_block.ou" >> $qsub_realign
	    echo "#PBS -e $TopOutputLogs/start_realrecal_block.in" >> $qsub_realign
	    echo "#PBS -q $pbsqueue" >> $qsub_realign
	    echo "#PBS -m ae" >> $qsub_realign
	    echo "#PBS -W depend=afterok:$pbsids" >> $qsub_realign
	    echo "#PBS -M $email" >> $qsub_realign
	    echo "$scriptdir/start_realrecal_block.sh $runfile $TopOutputLogs/start_realrecal_block.in $TopOutputLogs/start_realrecal_block.ou $email $TopOutputLogs/qsub.start_realrecal_block" >> $qsub_realign
	    `chmod a+r $qsub_realign` 
	    `qsub $qsub_realign >> $TopOutputLogs/pbs.REALRECAL`

            # need to release jobs here or realignment will not start
            `qrls -h u $alignids`
            `qrls -h u $mergeids`
            #`qrls -h u $extraids`

	    echo `date`
	fi

	`chmod -R 770 $AlignOutputDir`
	`chmod -R 770 $TopOutputLogs`
	echo `date`

