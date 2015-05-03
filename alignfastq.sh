#!/bin/bash
#
# originally written in collaboration with Mayo Bioinformatics core group
# alignfastq.sh
# align module to be used for input files in fastq format
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
else 
	set -x
	echo `date`
        scriptfile=$0
        runfile=$1
        elog=$2
        olog=$3
        email=$4
        qsubfile=$5
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi






# wrapping commends in echo, so that the output logs would be easier to read: they will have more structure
echo "####################################################################################################"
echo "#####################################                       ########################################"
echo "##################################### PARSING RUN INFO FILE ########################################"
echo "#####################################                       ########################################"
echo "####################################################################################################"




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
        igvdir=$( cat $runfile | grep -w IGVDIR | cut -d '=' -f2 )
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
        sampledir=$( cat $runfile | grep -w SAMPLEDIR | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        markduplicatestool=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        profiling=$( cat $runfile | grep -w PROFILING | cut -d '=' -f2 )
################## MUST PUT IN A CHECK FOR CORRECT SETTING ON THE PROFILING VARIABLE
        profiler=$( cat $runfile | grep -w PROFILER | cut -d '=' -f2 )
        cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )








echo "####################################################################################################"
echo "#####################################                       ########################################"
echo "#####################################  CREATE  DIRECTORIES  ########################################"
echo "#####################################                       ########################################"
echo "####################################################################################################"




        AlignOutputDir=$outputdir/align
        if [ -d $AlignOutputDir ]
        then
           echo "$AlignOutputDir is there; resetting it"
           `rm -r $AlignOutputDir/*`
        else
           mkdir -p $AlignOutputDir
        fi
        `chmod -R 770 $AlignOutputDir/`

        TopOutputLogs=$outputdir/logs
        if [ -d $TopOutputLogs ]
        then
           echo "$TopOutputLogs is there; resetting it"
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




        pipeid=$( cat $TopOutputLogs/CONFIGUREpbs )

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
        igv=$outputdir/$igvdir
        extradir=$outputdir/extractreads







echo "####################################################################################################"
echo "#####################################                       ########################################"
echo "#####################################   PARAMETER   CHECK   ########################################"
echo "#####################################                       ########################################"
echo "####################################################################################################"





############ checking workflow scripts directory
        echo -e "\n\n\n############ checking workflow scripts directory\n"

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $outputdir ]
        then
           mkdir -p $outputdir
        fi



############ checking samples 
        echo -e "\n\n\n############ checking samples\n"

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
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi

        # check that sampledir exists
        if [ ! -d $sampledir ]
        then
           MSG="SAMPLEDIR=$sampledir file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

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

        if [ $chunkfastq != "YES" -a $chunkfastq != "NO" ]
        then
            MSG="CHUNKFASTQ variable must be binary YES/NO; incorrect value encountered"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi





############ checking settings for marking of duplicates 
        echo -e "\n\n\n############ checking settings for marking of duplicates\n"

        if [ $dup != "1" -a $dup != "0" -a $dup != "YES" -a $dup != "NO" ]
        then
           MSG="Invalid value for MARKDUP=$dup"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
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
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
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





############ checking FastQC settings
        echo -e "\n\n\n############ checking FastQC settings\n"

        if [ $fastqcflag != "1" -a $fastqcflag != "0" -a $fastqcflag != "YES" -a $fastqcflag != "NO" ]
        then
           MSG="Invalid value for FASTQCFLAG=$fastqcflag"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit 1;
        else
            if [ $fastqcflag == "1" ]
            then
                $fastqcflag="YES"
            fi
            if [ $fastqcflag == "0" ]
            then
                $fastqcflag="NO"
            fi
        fi





############ checking Cleanup settings
        echo -e "\n\n\n############ checking Cleanup settings\n"

        if [ $cleanupflag != "1" -a $cleanupflag != "0" -a $cleanupflag != "YES" -a $cleanupflag != "NO" ]
        then
           MSG="Invalid value for REMOVETEMPFILES=$cleanupflag"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
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




############ checking computational tools
        echo -e "\n\n\n############ checking computational tools\n"

        if [ $aligner != "NOVOALIGN" -a $aligner != "BWA_ALN" -a $aligner != "BWA_MEM"]
        then
            MSG="ALIGNER=$aligner  is not available at this site"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
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
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           if [ $sortool != "NOVOSORT" -a $sortool != "PICARD" ]
           then
               MSG="Invalid value for SORTOOL=$sortool in configuration file"
               echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
               #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
               exit 1;
           fi
        fi      

        if [ -z $markduplicatestool ]
        then
           MSG="Value for MARKDUPLICATESTOOL must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           if [ $markduplicatestool != "PICARD" -a $markduplicatestool != "SAMBLASTER" ]
            then
               MSG="Invalid value for SORTOOL=$sortool in configuration file"
               echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
               #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
               exit 1;
           fi
        fi

        if [ ! -d $alignerdir ]
        then
           MSG="$alignerdir aligner directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        if [ ! -d $picardir ]
        then
           MSG="PICARDIR=$picardir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        if [ ! -d $samdir ]
        then
           MSG="SAMDIR=$samdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        if [ ! -e $profiler ]
        then
           MSG="PROFILER=$profiler not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi



############ checking presence of references 
        echo -e "\n\n\n############ checking presence of references\n"

        if [ ! -s $refdir/$ref ]
        then
           MSG="$refdir/$ref reference genome not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -s $refdir/$refindexed ] 
        then
           if  [ ! -s $refdir/$refindexed.fa ]
           then
              MSG="$refdir/$refindexed index for reference genome not found"
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
              exit 1;
           fi
        fi



############ checking run method
   if [ $run_method != "LAUNCHER" -a $run_method != "QSUB" -a $run_method != "APRUN" -a $run_method != "SERVER" ]
   then
      MSG="Invalid value for RUNMETHOD=$run_method"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi







############################################################################################################
#####################################                               ########################################
#################################### ALIGNMENT: LOOP OVER SAMPLES  ########################################
#####################################                               ########################################
############################################################################################################

echo -e "\n\n\n#################################### ALIGNMENT: LOOP OVER SAMPLES  ########################################\n\n\n"



        # clear the joblist
        truncate -s 0 $AlignOutputLogs/Anisimov.joblist

        while read SampleName
        do
          echo "processing next sample" 
          # this will evaluate the length of string
          if [ `expr ${#SampleName}` -lt 1 ]
          then
            echo "skipping empty line"
          else

            echo "aligning: $SampleName"

            # form and check the left reads file
	    LeftReadsFastq=$( ls $sampledir/${SampleName}_read1.* )
            if [ ! -e $LeftReadsFastq ]
            then
                    MSG="$LeftReadsFastq left reads file not found"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
                    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit 1;
            fi
	    echo `date`

## I used to check the number of lines in left and right fastq, and compare them, 
## before any further work is scheduled; but, this takes too long.
## wc -l takes 10-15 minutes on a 50X genome file;
## multiply that by two, then by number of genomes - and it takes hours
## just to do that checking; for 45 input genomes the alignfastq.sh
## takes 11:28:45. Very long. After all, if the files are not complete,
## then paired-ended bwa mem and novoalign will abort them; 
## then we can capture that error and compare with results of the FastQC.
## Same goes for when the input file is empty.
##
## Thus I have decided to keep these checks on the number of lines in fastqcommented out for now. 
##
#            NumLinesInLeftFastq=`wc -l $LeftReadsFastq | cut -d ' ' -f 1`
#	    echo `date`
#            if [ $NumLinesInLeftFastq -lt 1 ]
#            then
#                    MSG="$LeftReadsFastq left reads file is empty"
#                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
                    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
#                    exit 1;
#            fi

   
            # form and check the right reads file
            RightReadsFastq=$( ls $sampledir/${SampleName}_read2.* )
            if [ $paired -eq 1 ]
            then
               if [ ! -e $RightReadsFastq ]
               then
                  MSG="$RightReadsFastq right reads file not found"
                  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
                  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  exit 1;
#               else
#                  NumLinesInRightFastq=`wc -l $RightReadsFastq | cut -d ' ' -f 1`
#                  if [ $NumLinesInRightFastq -lt 1 ]
#                  then
#                      MSG="$RightReadsFastq right reads file is empty"
#                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
                      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
#                      exit 1;
#                  elif [ $NumLinesInRightFastq -ne $NumLinesInLeftFastq ]
#                  then 
#                      MSG="number of lines in $RightReadsFastq is not equal to that in $NumLinesInLeftFastq"
#                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
                   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
#                      exit 1;
#                  fi
               fi
            fi


            # form the read group values for the aligner 
            sID=$SampleName
            sPU=$SampleName
            sSM=$SampleName
            sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
            sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
            sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
            RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )





###########################   CHECK FASTQ QUALITY ###########################################

            if [ $fastqcflag == "YES" ]
            then
		echo "calculating quality values for fastq file"
                
                # check that fastqc tool is there
		if [ ! -d $fastqcdir ]
		then
		    MSG="FASTQCDIR=$fastqcdir directory not found"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
		fi

                # create necessary directories
		if [ ! -d $outputdir/fastqc ]
                then
                    mkdir $outputdir/fastqc
                    `chmod -R 770 $outputdir/fastqc`
		fi
                if [ ! -d $TopOutputLogs/fastqc ]
                then
                    mkdir $TopOutputLogs/fastqc
                    `chmod -R 770 $TopOutputLogs/fastqc`
                fi


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
		echo "aprun -n 1 -d $thr $scriptdir/fastq.sh $fastqcdir $outputdir/fastqc $fastqcparms $LeftReadsFastq $TopOutputLogs/fastqc/log.fastqc_R1_${SampleName}.in $TopOutputLogs/log.fastqc_R1_${SampleName}.ou $email $TopOutputLogs/fastqc/qsub.fastqc_R1_$SampleName" >> $qsub_fastqcR1
		`chmod a+r $qsub_fastqcR1`
                `qsub $qsub_fastqcR1 >> $TopOutputLogs/FASTQCpbs`

		if [ $paired -eq 1 ]
		then
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
		    echo "aprun -n 1 -d $thr $scriptdir/fastq.sh $fastqcdir $outputdir/fastqc $fastqcparms $RightReadsFastq $TopOutputLogs/fastqc/log.fastqc_R2_${SampleName}.in $TopOutputLogs/fastqc/log.fastqc_R2_$SampleName.ou $email $TopOutputLogs/qsub.fastqc_R2_$SampleName" >> $qsub_fastqcR2
		    `chmod a+r $qsub_fastqcR2`
                    `qsub $qsub_fastqcR2 >> $TopOutputLogs/FASTQCpbs`
		fi
            else
		echo "quality information for fastq files will NOT be calculated."
            fi
            
            ## done with generating quality info for each read file





###########################   CHUNK INPUT FASTQ  ###########################################


            # results of alignment for each sample will be placed into its own subfolder 
            if [ ! -d $AlignOutputDir/$SampleName ]
            then
		mkdir $AlignOutputDir/$SampleName
		outputsamfileprefix=$AlignOutputDir/$SampleName/$SampleName
	    else
		outputsamfileprefix=$AlignOutputDir/$SampleName/$SampleName
	    fi
	    `chmod -R 770 $AlignOutputDir/$SampleName/`
            # when running multiple samples via Anisimov, there will be lots of qsubs, logs and jobfiles
            # best keep them in the sample subfolder
            if [ ! -d $AlignOutputDir/$SampleName/logs ]
            then
		mkdir $AlignOutputDir/$SampleName/logs
            fi



            sortedplain=$outputsamfileprefix.wrg.sorted.bam
            outsortnodup=$outputsamfileprefix.nodups.sorted.bam
            outsortwdup=$outputsamfileprefix.wdups.sorted.bam




            cd $AlignOutputDir/$SampleName


            # create new names for chunks of fastq
            # doing this outside the chunkfastq conditional, because
            # non-chunking is handled by creating a symbolic link from chunk 0 to original fastq
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
	   	echo "mod is 0; no reads for last chunk file, one idle node"
	   	let NumChunks-=1
               fi

               echo "splitting read file1=$LeftReadsFastq"
               `split -l $NumLinesPerChunk -a 2 -d $LeftReadsFastq $LeftReadsChunkNamePrefix`
               exitcode=$?
               if [ $exitcode -ne 0 ]
               then
                      MSG="splitting of read file $LeftReadsFastq failed. exitcode=$exitcode"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                      exit $exitcode;
               fi
               if [ $paired -eq 1 ]
               then
                   echo "splitting read file2=$RightReadsFastq"
 	     	   `split -l $NumLinesPerChunk -a 2 -d $RightReadsFastq $RightReadsChunkNamePrefix`
	   	   exitcode=$?
	   	   if [ $exitcode -ne 0 ]
                   then
                       MSG="splitting of read file $RightReadsFastq failed.  exitcode=$exitcode"
                       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                       exit $exitcode;
	   	   fi
               fi
            else
               NumChunks=0
               echo "copying original fastq into a single chunk takes too long for whole genome"
               echo "instead, we will create symbolic link to the original file"
               echo "adding double-0 to the chunk name, because expecting to chunk files into tens of pieces"
               echo "need consistency in naming: two digits for chunk number"
               #   `cp $LeftReadsFastq ${LeftReadsChunkNamePrefix}0`
               ln -s $LeftReadsFastq leftreads_chunk00

               if [ $paired -eq 1 ]
               then
               #   `cp $RightReadsFastq ${RightReadsChunkNamePrefix}0`
                  ln -s $RightReadsFastq rightreads_chunk00
               fi
            fi

            ## done chunking input fastq




###########################   FORM ALIGNMENT QSUBS  ###########################################

            # if want to profile, then must set up  the environment
            # for now, let us do without this; too much else going on: Aug 22, 2014
            # in other words, do not set the option in the runfile
            #if [ $profiling == 'memprof' ]
            #then
            #   ccmgres_string="#PBS -l gres=ccm"
            #   run_string="module add ccm; ccmrun "
            #   profiler_string="$profiler "
            #else
            #   ccmgres_string=""
            #   run_string="aprun -n 1 -d $thr "
            #   profiler_string=""
            #fi
            

            allfiles=""
            # begin loop over chunks of the current input fastq
            for i in $(seq 0 $NumChunks)
            do
                echo "step 1: aligning chunk $i... "
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
                   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
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
                        #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
			exit 1;
                    fi
                fi
                if [ $aligner == "NOVOALIGN"  ]
		then
                    echo -e "\n###############   novoalign is used as aligner. input file in fastq format ################\n"
                    #qsub_novosplit=$AlignOutputLogs/qsub_novosplit.novoaln.$SampleName.node$i
                    #echo "#PBS -V" > $qsub_novosplit
                    #echo "#PBS -A $pbsprj" >> $qsub_novosplit
                    #echo "#PBS -N ${pipeid}_novoaln_${SampleName}_$i" >> $qsub_novosplit
		    #echo "#PBS -l walltime=$pbscpu" >> $qsub_novosplit
		    #echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_novosplit
		    #echo $ccmgres_string >> $qsub_novosplit
		    #echo "#PBS -o $AlignOutputLogs/log.novoaln.$SampleName.node$i.ou" >> $qsub_novosplit
		    #echo "#PBS -e $AlignOutputLogs/log.novoaln.$SampleName.node$i.in" >> $qsub_novosplit
                    #echo "#PBS -q $pbsqueue" >> $qsub_novosplit
                    #echo "#PBS -m ae" >> $qsub_novosplit
                    #echo "#PBS -M $email" >> $qsub_novosplit
                    if [ $paired -eq 1 ]
                    then
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $scriptdir $runfile $paired $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.novoaln.$SampleName.node$OutputFileSuffix" >> $AlignOutputDir/${SampleName}/logs/novosplit.${SampleName}.node$OutputFileSuffix
                    else
			echo "$scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $scriptdir $runfile $paired $AlignOutputDir/$SampleName/$Rone $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.novoaln.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.novoaln.$SampleName.node$OutputFileSuffix" >> $AlignOutputDir/${SampleName}/logs/novosplit.${SampleName}.node$OutputFileSuffix
                    fi
                    #`chmod a+r $qsub`_novosplit
                    #jobnovo=`qsub $qsub`_novosplit
                    #`qhold -h u $jobnovo`
		    #echo $jobnovo >> $AlignOutputLogs/ALIGNED_$SampleName
		elif [ $aligner == "BWA_ALN" ] 
                then
                    echo "bwa is used as aligner. input file format is in fastq"
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
		    echo "aprun -n 1 -d $thr $scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.R1.sai $AlignOutputDir/$SampleName/$Rone $scriptdir $AlignOutputLogs/log.bwar1.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwar1.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwar1.$SampleName.node$OutputFileSuffix" >> $qsub_bwar1

                    `chmod a+r $qsub_bwar1`
                    jobr1=`qsub $qsub_bwar1`
                    `qhold -h u $jobr1`
                    echo $jobr1 >> $AlignOutputLogs/ALIGNED_$SampleName
                    if [ $paired -eq 1 ]
                    then
                        echo "bwa aligner. paired-end reads"
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
			echo "$run_string $profiler_string $scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.R2.sai $AlignOutputDir/$SampleName/$Rtwo $scriptdir $AlignOutputLogs/log.bwar2.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwar2.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwar2.$SampleName.node$OutputFileSuffix" >> $qsub_bwar2
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
			echo "$run_string $profiler_string $scriptdir/bwaS2.sh $alignerdir $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.R1.sai $outputsamfileprefix.node$OutputFileSuffix.R2.sai $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $samdir $AlignOutputLogs/log.bwasampe.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwasampe.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwasampe.$SampleName.node$OutputFileSuffix" >> $qsub_bwasampe
			`chmod a+r $qsub_bwasampe`
                        jobwa=`qsub $qsub_bwasampe`
			`qhold -h u $jobwa`
			echo $jobwa >> $AlignOutputLogs/ALIGNED_$SampleName
                    else
                        echo "bwa aligner. single read"
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
			echo "aprun -n 1 -d $thr $scriptdir/bwaS3.sh $alignerdir $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.R1.sai $AlignOutputDir/$SampleName/$Rone $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $samdir $AlignOutputLogs/log.bwasamse.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwasamse.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwasamse.$SampleName.node$OutputFileSuffix" >> $qsub_bwasamse
			`chmod a+r $qsub_bwasamse`
                        jobwa=`qsub $qsub_bwasamse`
			`qhold -h u $jobwa`
                        echo $qsub_bwasamse >> $AlignOutputLogs/ALIGNED_$SampleName
                    fi
                elif [ $aligner == "BWA_MEM" ]
                then
                    echo "bwa mem is used as aligner. input file format is in fastq"
                    if [ $paired -eq 1 ]
                    then
                        echo "bwa mem aligner. paired-end reads"
                        jobfile=$AlignOutputDir/$SampleName/logs/bwamem.$SampleName.node$OutputFileSuffix.jobfile

                        if [ $chunkfastq == "YES" ]
                        then
			   ### MUST FIX QSUB VARIABLE NAME PROPERLY
                           echo "$run_string $profiler_string $scriptdir/bwamem_pe.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $scriptdir $samdir $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix" >> $qsub3
                        elif [ $chunkfastq == "NO" ]
                        then
                           echo "nohup $profiler_string $scriptdir/bwamem_pe_markduplicates.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $runfile $AlignOutputDir/$SampleName/logs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputDir/$SampleName/logs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $jobfile $RGparms $AlignOutputLogs > $AlignOutputDir/$SampleName/logs/log.bwamem.$SampleName.node$OutputFileSuffix.in" > $jobfile
                           jobfilename=$( basename $jobfile )
                           echo "$AlignOutputDir/$SampleName/logs $jobfilename" >> $AlignOutputLogs/AlignAnisimov.joblist
                        fi


                     else
                        echo "bwa mem aligner. single-end reads"
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
                           echo "aprun -n 1 -d $thr $scriptdir/bwamem_se.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.sam $outputsamfileprefix.node$OutputFileSuffix.bam $AlignOutputDir/$SampleName/$Rone $scriptdir $samdir $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix" >> $qsub_bwamem
                        elif [ $chunkfastq == "NO" ]
                        then
                           echo "aprun -n 1 -d $thr $scriptdir/bwamem_se_markduplicates.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$OutputFileSuffix.bam $AlignOutputDir/$SampleName/$Rone $scriptdir $samdir $samblasterdir $picardir $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.in $AlignOutputLogs/log.bwamem.$SampleName.node$OutputFileSuffix.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$OutputFileSuffix $RGparms" >> $qsub_bwamem
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
            done
            # end loop over chunks of the current input fastq





##########################################################################################################################################
###########################                                                                    ###########################################
echo -e "\n\n\n ###########   FORM POST-ALIGNMENT QSUBS: MERGINE, SORTING, MARKING DUPLICATES  ################################### \n\n\n"
###########################                                                                    ###########################################
##########################################################################################################################################

            cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/ALIGNEDpbs

            if [ $chunkfastq == "YES" ]
            then
   	       #ALIGNED=$( cat $AlignOutputLogs/ALIGNED_* | sed "s/\.[a-z]*//" | tr "\n" ":" )
	       #ALIGNED=$( cat $AlignOutputLogs/ALIGNED_* | sed "s/\..*//" | tr "\n" ":" )

	       listfiles=$( echo $allfiles  | tr " " ":" | sed "s/::/:/g" )
               if [ $sortool == "NOVOSORT" ]
               then
                   echo -e "\n merging aligned chunks with novosort\n"
		   #qsub_sortmerge=$AlignOutputLogs/qsub.sortmerge.novosort.$SampleName
		   #echo "#PBS -V" > $qsub_sortmerge
		   #echo "#PBS -A $pbsprj" >> $qsub_sortmerge
		   #echo "#PBS -N ${pipeid}_sortmerge_novosort_$SampleName" >> $qsub_sortmerge
		   #echo "#PBS -l walltime=$pbscpu" >> $qsub_sortmerge
		   #echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_sortmerge
		   #echo "#PBS -o $AlignOutputLogs/log.sortmerge_novosort.$SampleName.ou" >> $qsub_sortmerge
		   #echo "#PBS -e $AlignOutputLogs/log.sortmerge_novosort.$SampleName.in" >> $qsub_sortmerge
		   #echo "#PBS -q $pbsqueue" >> $qsub_sortmerge
		   #echo "#PBS -m ae" >> $qsub_sortmerge
		   #echo "#PBS -M $email" >> $qsub_sortmerge
		   #echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub_sortmerge
		   echo "$scriptdir/mergenovo.sh $AlignOutputDir/$SampleName $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $AlignOutputLogs/log.sortmerge_novosort.$SampleName.in $AlignOutputLogs/log.sortmerge_novosort.$SampleName.ou $email $AlignOutputLogs/qsub.sortmerge.novosort.$SampleName" >> $AlignOutputDir/${SampleName}/logs/mergenovo.${SampleName}
		   #`chmod a+r $qsub_sortmerge`
		   #mergejob=`qsub $qsub_sortmerge`
		   #`qhold -h u $mergejob`
		   #echo $mergejob  > $AlignOutputLogs/MERGED_$SampleName
               else
                   echo "merging aligned chunks with picard"
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
		   echo "aprun -n 1 -d $thr $scriptdir/mergepicard.sh $AlignOutputDir/$SampleName $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $AlignOutputLogs/log.sortmerge.picard.$SampleName.in $AlignOutputLogs/log.sortmerge.picard.$SampleName.ou $email $AlignOutputLogs/qsub.sortmerge.picard.$SampleName" >> $qsub_sortmerge
		   `chmod a+r $qsub_sortmerge`
		   mergejob=`qsub $qsub_sortmerge`
		   `qhold -h u $mergejob`
		   echo $mergejob  > $AlignOutputLogs/MERGED_$SampleName
               fi

               ## COMMENTING THIS OUT AS A MAYO IDIOSYNCRACY
	       #echo `date`
	       #echo -e "\n extract reads specified in CHRINDEX param\n"
	       #qsub_extractreads=$AlignOutputLogs/qsub.extractreadsbam.$SampleName
   	       #echo "#PBS -V" > $qsub_extractreads
	       #echo "#PBS -A $pbsprj" >> $qsub_extractreads
	       #echo "#PBS -N ${pipeid}_extrbam_$SampleName" >> $qsub_extractreads
	       #echo "#PBS -l walltime=$pbscpu" >> $qsub_extractreads
	       #echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_extractreads
	       #echo "#PBS -o $AlignOutputLogs/log.extractreadsbam.$SampleName.ou" >> $qsub_extractreads
	       #echo "#PBS -e $AlignOutputLogs/log.extractreadsbam.$SampleName.in" >> $qsub_extractreads
	       #echo "#PBS -q $pbsqueue" >> $qsub_extractreads
	       #echo "#PBS -m ae" >> $qsub_extractreads
	       #echo "#PBS -M $email" >> $qsub_extractreads
	       #echo "#PBS -W depend=afterok:$mergejob" >> $qsub_extractreads
	       #echo "aprun -n 1 -d $thr $scriptdir/extract_reads_bam.sh $AlignOutputDir/$SampleName $outsortwdup $runfile $AlignOutputLogs/log.extractreadsbam.$SampleName.in $AlignOutputLogs/log.extractreadsbam.$SampleName.ou $email  $AlignOutputLogs/qsub.extractreadsbam.$SampleName $igv $extradir" >> $qsub_extractreads
	       #`chmod a+r $qsub_extractreads`
	       #extrajob=`qsub $qsub_extractreads`
               #`qhold -h u $extrajob`
               #echo $extrajob >> $TopOutputLogs/EXTRACTREADSpbs

	       #cat $AlignOutputLogs/MERGED_$SampleName >> $TopOutputLogs/MERGEDpbs
            fi
          fi
          (( inputfastqcounter++ )) # was not initialized at beginning of loop, so starts at zero
	done <  $outputdir/SAMPLENAMES.list
        # end loop over input fastq






#####################################################################################################################################
#####################################                                                        ########################################
echo -e "\n\n\n######################  SCHEDULE NOVOALIGN and MERGENOVO QSUBS CREATED ABOVE  #################################\n\n\n"
#####################################                                                        ########################################
#####################################################################################################################################


        if [ $chunkfastq == "YES" -a $aligner == "NOVOALIGN" -a $sortool == "NOVOSORT" ]
        then

           case $run_method in
           "APRUN")
              while read SampleName
              do
                 # ADD LOOP OVER CHUNKS
                 for i in $(seq 0 $NumChunks)
                 do

                    if (( $i < 10 ))
                    then
                       OutputFileSuffix=0${i}
                    else
                       OutputFileSuffix=${i}
                    fi

                    echo -e "\n\nscheduling qsubs\n\n"
                    qsub_novosplit=$AlignOutputDir/${SampleName}/logs/qsub.novosplit.${SampleName}.node${OutputFileSuffix}
                    # appending the generic header to the qsub
                    cat $outputdir/qsubGenericHeader > $qsub_novosplit


                    ###############################
                    echo -e "\n################# constructing qsub for novosplit\n"
                    echo "#PBS -N ${pipeid}_novoalign.${SampleName}.node${OutputFileSuffix}" >> $qsub_novosplit
                    echo "#PBS -l walltime=$pbscpu" >> $qsub_novosplit
                    echo "#PBS -o $AlignOutputDir/${SampleName}/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.ou" >> $qsub_novosplit
                    echo "#PBS -e $AlignOutputDir/${SampleName}/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.in" >> $qsub_novosplit
                    echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_novosplit

                    # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                    sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/${SampleName}/logs/novosplit.${SampleName}.node$OutputFileSuffix >> $qsub_novosplit 

                    novosplit_job=`qsub $qsub_novosplit`
                    `qhold -h u $novosplit_job`
                    echo $novosplit_job >> $TopOutputLogs/ALIGNEDpbs_${SampleName} # so that this job could be released in the next section. Should it be held to begin with?

                 done # done looping over chunks of a sample
                 cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/ALIGNEDpbs
                 alignids=$( cat $TopOutputLogs/ALIGNED_$SampleName | sed "s/\..*//" | tr "\n" ":" )

                 qsub_mergenovo=$AlignOutputDir/${SampleName}/logs/qsub.mergenovo.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_mergenovo




                 ###############################
                 echo -e "\n ################# constructing qsub for mergenovo\n"
                 echo "#PBS -N ${pipeid}_mergenovo" >> $qsub_mergenovo
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_mergenovo
                 echo "#PBS -o $AlignOutputLogs/log.mergenovo.ou" >> $qsub_mergenovo
                 echo "#PBS -e $AlignOutputLogs/log.mergenovo.in" >> $qsub_mergenovo
                 # add dependency on novosplit job
                 echo -e "#PBS -W depend=afterok:$alignids" >> $qsub_mergenovo
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_mergenovo

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/${SampleName}/logs/mergenovo.${SampleName} >> $qsub_mergenovo 

                 mergenovo_job=`qsub $qsub_mergenovo`
                 `qhold -h u $mergenovo_job`
                 echo $mergenovo_job > $TopOutputLogs/MERGEDpbs # so that this job could be released in the next section, and start_realrecal_block.sh could depend on it

                 # release all held novosplit jobs for this sample
                 `qrls -h u $alignids`

              done < $outputdir/SAMPLENAMES.list # done looping over samples
           ;;
           "QSUB")
              while read SampleName
              do
                 for i in $(seq 0 $NumChunks)
                 do

                    if (( $i < 10 ))
                    then
                       OutputFileSuffix=0${i}
                    else
                       OutputFileSuffix=${i}
                    fi

                    echo -e "\n\nscheduling qsubs\n\n"
                    qsub_novosplit=$AlignOutputDir/${SampleName}/logs/qsub.novosplit.${SampleName}.node${OutputFileSuffix}

                    # appending the generic header to the qsub
                    cat $outputdir/qsubGenericHeader > $qsub_novosplit


                    ###############################
                    echo -e "\n################# constructing qsub for novosplit\n"
                    echo "#PBS -N ${pipeid}_novoalign.${SampleName}.node${OutputFileSuffix}" >> $qsub_novosplit
                    echo "#PBS -l walltime=$pbscpu" >> $qsub_novosplit
                    echo "#PBS -o $AlignOutputDir/${SampleName}/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.ou" >> $qsub_novosplit
                    echo "#PBS -e $AlignOutputDir/${SampleName}/logs/log.novosplit.${SampleName}.node${OutputFileSuffix}.in" >> $qsub_novosplit
                    echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_novosplit
                    cat $AlignOutputDir/${SampleName}/logs/novosplit.${SampleName}.node${OutputFileSuffix} >> $qsub_novosplit
                    novosplit_job=`qsub $qsub_novosplit`
                    `qhold -h u $novosplit_job`
                    echo $novosplit_job >> $TopOutputLogs/ALIGNEDpbs # so that this job could be released in the next section. Should it be held to begin with?

                 done # done looping over chunks of a sample
                 cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/ALIGNEDpbs
                 alignids=$( cat $TopOutputLogs/ALIGNED_$SampleName | sed "s/\..*//" | tr "\n" ":" )

                 qsub_mergenovo=$AlignOutputDir/${SampleName}/logs/qsub.mergenovo.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_mergenovo




                 ###############################
                 echo -e "\n ################# constructing qsub for mergenovo\n"
                 echo "#PBS -N ${pipeid}_mergenovo" >> $qsub_mergenovo
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_mergenovo
                 echo "#PBS -o $AlignOutputLogs/log.mergenovo.ou" >> $qsub_mergenovo
                 echo "#PBS -e $AlignOutputLogs/log.mergenovo.in" >> $qsub_mergenovo
                 # add dependency on novosplit job
                 echo -e "#PBS -W depend=afterok:$alignids" >> $qsub_mergenovo
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_mergenovo

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/${SampleName}/logs/mergenovo.${SampleName} >> $qsub_mergenovo

                 mergenovo_job=`qsub $qsub_mergenovo`
                 `qhold -h u $mergenovo_job`
                 echo $mergenovo_job > $TopOutputLogs/MERGEDpbs # so that this job could be released in the next section, and start_realrecal_block.sh could depend on it

                 # release all held novosplit jobs for this sample
                 `qrls -h u $alignids`


              done < $outputdir/SAMPLENAMES.list # done looping over samples
           ;;
           "LAUNCHER")
              # clear out the joblists
              truncate -s 0 $AlignOutputLogs/novosplit.AnisimovJoblist
              truncate -s 0 $AlignOutputLogs/mergenovo.AnisimovJoblist

              while read SampleName
              do
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
                    novosplit_log=${AlignOutputDir}/${SampleName}/logs/log.novosplit.${SampleName}.node$OutputFileSuffix.in
                    awk -v awkvar_novosplitlog=$novosplit_log '{print "nohup "$0" > "awkvar_novosplitlog}' $AlignOutputDir/${SampleName}/logs/novosplit.${SampleName}.node${OutputFileSuffix} > $AlignOutputDir/${SampleName}/logs/jobfile.novosplit.${SampleName}.node${OutputFileSuffix}
                    echo "$AlignOutputDir/${SampleName}/logs/ jobfile.novosplit.${SampleName}.node${OutputFileSuffix}" >> $AlignOutputLogs/novosplit.AnisimovJoblist
                 done # done going through chunks

                 mergenovo_log=${AlignOutputDir}/${SampleName}/logs/log.mergenovo.${SampleName}.in
                 awk -v awkvar_mergenovolog=$mergenovo_log '{print "nohup "$0" > "awkvar_mergenovolog}' $AlignOutputDir/${SampleName}/logs/mergenovo.${SampleName} > $AlignOutputDir/${SampleName}/logs/jobfile.mergenovo.${SampleName}
                 echo "$AlignOutputDir/${SampleName}/logs/ jobfile.mergenovo.${SampleName}" >> $AlignOutputLogs/mergenovo.AnisimovJoblist
              done < $outputdir/SAMPLENAMES.list # done looping over samples


              echo -e "\n\nscheduling Anisimov Launcher joblists\n\n"
              qsub_novosplit_anisimov=$AlignOutputLogs/qsub.novosplit.AnisimovLauncher
              qsub_mergenovo_anisimov=$AlignOutputLogs/qsub.mergenovo.AnisimovLauncher

              # appending the generic header to the qsub
              cat $outputdir/qsubGenericHeader > $qsub_novosplit_anisimov
              cat $outputdir/qsubGenericHeader > $qsub_mergenovo_anisimov


              ###############################
              echo -e "\n################# constructing qsub for novosplit\n"
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
              echo $novosplit_job >> $TopOutputLogs/ALIGNEDpbs # so that this job could be released in the next section. Should it be held to begin with?


              ###############################
              echo -e "\n ################# constructing qsub for mergenovo\n"
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
              echo $mergenovo_job > $TopOutputLogs/MERGEDpbs # so that this job could be released in the next section, and start_realrecal_block.sh could depend on it 

           ;;
           "SERVER")
              # will fill in later
           ;;
           esac
        else 
           MSG="something went wrong with scheduling alignment: chunking/aligner/sort_tool options not correctly set"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
        fi








#####################################################################################################################
#####################################                                        ########################################
#####################################  SCHEDULE BWA-MEM QSUBS CREATED ABOVE  ########################################
#####################################                                        ########################################
#####################################################################################################################


        if [ $chunkfastq == "NO" ]
        then

           case $run_method in
           "APRUN")
              while read SampleName
              do
                 echo -e "\n\nscheduling qsubs\n\n"
                 qsub_bwa=$AlignOutputDir/${SampleName}/logs/qsub.bwa.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_bwa


                 ###############################
                 echo -e "\n################# constructing qsub for novosplit\n"
                 echo "#PBS -N ${pipeid}_bwa.${SampleName}" >> $qsub_bwa
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_bwa
                 echo "#PBS -o $AlignOutputDir/${SampleName}/logs/log.bwa.${SampleName}.ou" >> $qsub_bwa
                 echo "#PBS -e $AlignOutputDir/${SampleName}/logs/log.bwa.${SampleName}.in" >> $qsub_bwa
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_bwa

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 sed "1!b;s/^/aprun -n 1 -N 1 -d $thr /" $AlignOutputDir/${SampleName}/logs/bwamem.${SampleName}.node$OutputFileSuffix.jobfile >> $qsub_bwa

                 bwa_job=`qsub $qsub_bwa`
                 `qhold -h u $bwa_job`
                 echo $bwa_job >> $TopOutputLogs/ALIGNEDpbs # so that this job could be released in the next section. Should it be held to begin with?
              done
           ;;
           "QSUB")
              while read SampleName
              do
                 echo -e "\n\nscheduling qsubs\n\n"
                 qsub_bwa=$AlignOutputDir/${SampleName}/logs/qsub.bwa.${SampleName}
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_bwa


                 ###############################
                 echo -e "\n################# constructing qsub for novosplit\n"
                 echo "#PBS -N ${pipeid}_bwa.${SampleName}" >> $qsub_bwa
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_bwa
                 echo "#PBS -o $AlignOutputDir/${SampleName}/logs/log.bwa.${SampleName}.ou" >> $qsub_bwa
                 echo "#PBS -e $AlignOutputDir/${SampleName}/logs/log.bwa.${SampleName}.in" >> $qsub_bwa
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_bwa

                 # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
                 echo $AlignOutputDir/${SampleName}/logs/bwamem.${SampleName}.node$OutputFileSuffix.jobfile >> $qsub_bwa

                 bwa_job=`qsub $qsub_bwa`
                 `qhold -h u $bwa_job`
                 echo $bwa_job >> $TopOutputLogs/ALIGNEDpbs # so that this job could be released in the next section. Should it be held to begin with?
              done
           ;;
           "LAUNCHER")

              # increment so number of nodes = number fo input fastq + 1, even when there is only one input fastq
              # otherwise nowhere for launcher to run, due to the -N 1 option in aprun
              numalignnodes=$(( inputfastqcounter + 1))
   
              # scheduling the Anisimov Launcher
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
   
              # if not planning to profile in the launcher (string is empty)
              # then schedule aprun as usual
   #           if [ -z $profiler_string ]
   #           then
                 echo "aprun -n $numalignnodes -N 1 -d $thr ~anisimov/scheduler/scheduler.x $AlignOutputLogs/AlignAnisimov.joblist /bin/bash > $AlignOutputLogs/AlignAnisimov.joblist.log" >> $qsubAlignLauncher
              # otherwise use ccm
   #           else 
   #              echo $ccmgres_string >> $qsubAlignLauncher
   #              echo "$run_string ~anisimov/scheduler/scheduler.x $AlignOutputLogs/AlignAnisimov.joblist /bin/bash > $AlignOutputLogs/AlignAnisimov.joblist.log" >> $qsubAlignLauncher
   #           fi

              AlignAnisimovJoblistId=`qsub $qsubAlignLauncher`

              echo $AlignAnisimovJoblistId >> $TopOutputLogs/ALIGNEDpbs # so that this job could be released in the next section. Should it be held to begin with?
              echo $AlignAnisimovJoblistId > $TopOutputLogs/MERGEDpbs # so that summaryok and start_realrecal_block.sh could depend on this job, in case when there is no merging: a sigle chunk
           ;;
           esac

        fi




########################################################################################################
#####################################                           ########################################
#####################################  WRAP UP ALIGNMENT BLOCK  ########################################
#####################################                           ########################################
########################################################################################################

        
	pbsids=$( cat $TopOutputLogs/MERGEDpbs | sed "s/\..*//" | tr "\n" ":" )
	#extraids=$( cat $TopOutputLogs/EXTRACTREADSpbs | sed "s/\..*//" | tr "\n" " " )
        mergeids=$( echo $pbsids | tr ":" " " )
        alignids=$( cat $TopOutputLogs/ALIGNEDpbs | sed "s/\..*//" | tr "\n" " " )

        ## generating summary redmine email if analysis ends here
	echo "wrap up and produce summary table if analysis ends here or call realign if analysis continues"
	if [ $analysis == "ALIGNMENT" -o $analysis == "ALIGN" -o $analysis == "ALIGN_ONLY" ]
	then
            # release all held jobs
            `qrls -h u $alignids`
            `qrls -h u $mergeids`
            `qrls -h u $extraids`
     
	    lastjobid=""
            cleanjobid=""

            if [ $cleanupflag == "YES" ]
            then 
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
		echo $cleanjobid >> $outputdir/logs/CLEANUPpbs
            fi

            `sleep 30s`
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
	    echo "$scriptdir/summary.sh $outputdir $email exitok $reportticket"  >> $qsub_summary
	    `chmod a+r $qsub_summary`
	    lastjobid=`qsub $qsub_summary`
	    echo $lastjobid >> $TopOutputLogs/SUMMARYpbs

	    if [ `expr ${#lastjobid}` -lt 1 ]
	    then
		echo "at least one job aborted"
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
		echo "$scriptdir/summary.sh $outputdir $email exitnotok $reportticket"  >> $qsub_summary
		`chmod a+r $qsub_summary`
		badjobid=`qsub $qsub_summary`
		echo $badjobid >> $TopOutputLogs/SUMMARYpbs
	    fi
	fi

	if [ $analysis == "REALIGNMENT" -o $analysis == "REALIGN" ]
	then
            echo " analysis continues with realignment"
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
	    `qsub $qsub_realign >> $TopOutputLogs/REALRECALpbs`

            # need to release jobs here or realignment will not start
            `qrls -h u $alignids`
            `qrls -h u $mergeids`
            `qrls -h u $extraids`
	    echo `date`
	fi

	`chmod -R 770 $AlignOutputDir`
	`chmod -R 770 $TopOutputLogs`
	echo `date`
fi
