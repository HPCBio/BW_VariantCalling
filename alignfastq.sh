#!/bin/sh
#
# originally written in collaboration with Mayo Bioinformatics core group
# alignfastq.sh
# align module to be used for input files in fastq format
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=lmainzer@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
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







####################################################################################################
#####################################                       ########################################
##################################### PARSING RUN INFO FILE ########################################
#####################################                       ########################################
####################################################################################################



echo -e "\n\n\n##################################### PARSING RUN INFO FILE ########################################\n\n\n"


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
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        sampledir=$( cat $runfile | grep -w SAMPLEDIR | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        markduplicatestool=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        profiling=$( cat $runfile | grep -w PROFILING | cut -d '=' -f2 )
        profiler=$( cat $runfile | grep -w PROFILER | cut -d '=' -f2 )
        cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )








####################################################################################################
#####################################                       ########################################
#####################################  CREATE  DIRECTORIES  ########################################
#####################################                       ########################################
####################################################################################################

echo -e "\n\n\n#####################################  CREATE  DIRECTORIES  ########################################\n\n\n"



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




        pipeid=$( cat $TopOutputLogs/MAINpbs )

        chunks=`expr $nodes "-" 1`
        if [ $NumChunks -lt 1 ]
        then
            chunks=$nodes
        fi
        nthreads=`expr $thr "-" 1`
        if [ $nthreads -lt 1 ]
        then
            nthreads=$thr
        fi
        igv=$outputdir/$igvdir
        extradir=$outputdir/extractreads







####################################################################################################
#####################################                       ########################################
#####################################   PARAMETER   CHECK   ########################################
#####################################                       ########################################
####################################################################################################





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

        if [ $type == "GENOME" -o $type == "WHOLE_GENOME" -o $type == "WHOLEGENOME" -o $type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $type == "EXOME" -o $type == "WHOLE_EXOME" -o $type == "WHOLEEXOME" -o $type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for TYPE=$type in configuration file."
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

        # construct a list of samplenames, check that files actually exist
        numsamples=0
        truncate -s 0 $TopOutputLogs/SAMPLENAMES.tmp.list
        for fastqfile in $sampledir/*
        do
            # strip path, which read (left/right), and extension from input files
            # and put that info into the samplenames file
            samplename=$( basename $fastqfile | sed 's/_read.\?\.[^.]*$//' )
            echo -e "$samplename" >> $TopOutputLogs/SAMPLENAMES.tmp.list
            let numsamples+=1
        done
        # paired-ended fastq will produce duplicate lines in the samplenames file, so remove the duplicates
        uniq  $TopOutputLogs/SAMPLENAMES.tmp.list >  $TopOutputLogs/SAMPLENAMES.list
        rm  $TopOutputLogs/SAMPLENAMES.tmp.list

        if [ $numsamples -lt 1 ]
            then
              MSG="No samples found in SAMPLEDIR=$sampledir."
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
              #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
              exit 1;
            fi
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

        if [ $aligner != "NOVOALIGN" -a $aligner != "BWA" -a $aligner != "BWA_MEM"]
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
            alignparms=$( cat $runfile | grep -w NOVOPARAMS | cut -d '=' -f2 | tr " " "_" )_-c_${thr}
        fi
        if [ $aligner == "BWA" ]
        then
            alignerdir=$( cat $runfile | grep -w BWADIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w BWAPARAMS | cut -d '=' -f2 | tr " " "_" )_-t_${thr}
        fi
        if [ $aligner == "BWA_MEM" ]
        then
            alignerdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
            refindexed=$( cat $runfile | grep -w BWAMEMINDEX | cut -d '=' -f2 )
            alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 | tr " " "_" )_-t_${thr}
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
          if [ `expr ${#SampleName}` -lt 7 ]
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
            if [ ! -d $AlignOutputDir/$SampleName/AnisimovLogs ]
            then
		mkdir $AlignOutputDir/$SampleName/AnisimovLogs
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
               if [ $NumChunks -lt 1 ]
               then
                   NumChunks=1
               fi

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
               `split -l $NumLinesPerChunk -a 1 -d $LeftReadsFastq $LeftReadsChunkNamePrefix`
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
 	     	   `split -l $NumLinesPerChunk -a 1 -d $RightReadsFastq $RightReadsChunkNamePrefix`
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
               #   `cp $LeftReadsFastq ${LeftReadsChunkNamePrefix}0`
               ln -s $LeftReadsFastq leftreads_chunk0

               if [ $paired -eq 1 ]
               then
               #   `cp $RightReadsFastq ${RightReadsChunkNamePrefix}0`
                  ln -s $RightReadsFastq rightreads_chunk0
               fi
            fi

            ## done chunking input fastq




###########################   FORM ALIGNMENT QSUBS  ###########################################

            # if want to profile, then must set up  the environment
            if [ $profiling == 'memprof' ]
            then
               ccmgres_string="#PBS -l gres=ccm"
               run_string="module add ccm; ccmrun "
               profiler_string="$profiler "
            else
               ccmgres_string=""
               run_string="aprun -n 1 -d $thr "
               profiler_string=""
            fi
            

            allfiles=""
            # begin loop over chunks of the current input fastq
            for i in $(seq 0 $NumChunks)
            do
                echo "step 1: aligning chunk $i... "
		echo `date`
                Rone=${LeftReadsChunkNamePrefix}$i
                if [ ! -s $Rone ]
                then
                   MSG="chunk $i of read file $LeftReadsFastq file not found"
                   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                   exit 1;
                fi
		if [ $paired -eq 1 ]
		then
                    Rtwo=${RightReadsChunkNamePrefix}$i
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
                    echo "novoalign is used as aligner. input file in fastq format"
                    qsub=$AlignOutputLogs/qsub.novoaln.$SampleName.node$i
                    echo "#PBS -V" > $qsub
                    echo "#PBS -A $pbsprj" >> $qsub
                    echo "#PBS -N ${pipeid}_novoaln_${SampleName}_$i" >> $qsub
		    echo "#PBS -l walltime=$pbscpu" >> $qsub
		    echo "#PBS -l nodes=1:ppn=$thr" >> $qsub
		    echo $ccmgres_string >> $qsub
		    echo "#PBS -o $AlignOutputLogs/log.novoaln.$SampleName.node$i.ou" >> $qsub
		    echo "#PBS -e $AlignOutputLogs/log.novoaln.$SampleName.node$i.in" >> $qsub
                    echo "#PBS -q $pbsqueue" >> $qsub
                    echo "#PBS -m ae" >> $qsub
                    echo "#PBS -M $email" >> $qsub
                    if [ $paired -eq 1 ]
                    then
			echo "$run_string $profiler_string $scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.sam $outputsamfileprefix.node$i.bam $scriptdir $samdir $paired $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $AlignOutputLogs/log.novoaln.$SampleName.node$i.in $AlignOutputLogs/log.novoaln.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.novoaln.$SampleName.node$i" >> $qsub
                    else
			echo "$run_string $profiler_string $scriptdir/novosplit.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.sam $outputsamfileprefix.node$i.bam $scriptdir $samdir $paired $AlignOutputDir/$SampleName/$Rone $AlignOutputLogs/log.novoaln.$SampleName.node$i.in $AlignOutputLogs/log.novoaln.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.novoaln.$SampleName.node$i" >> $qsub
                    fi
                    `chmod a+r $qsub`
                    jobnovo=`qsub $qsub`
                    `qhold -h u $jobnovo`
		    echo $jobnovo >> $AlignOutputLogs/ALIGNED_$SampleName
		elif [ $aligner == "BWA" ] 
                then
                    echo "bwa is used as aligner. input file format is in fastq"
                    qsub1=$AlignOutputLogs/qsub.bwar1.$SampleName.node$i
                    echo "#PBS -V" > $qsub1
                    echo "#PBS -N ${pipeid}_bwar1_${SampleName}_$i" >> $qsub1
		    echo "#PBS -o $AlignOutputLogs/log.bwar1.$SampleName.node$i.ou" >> $qsub1
		    echo "#PBS -e $AlignOutputLogs/log.bwar1.$SampleName.node$i.in" >> $qsub1
                    echo "#PBS -A $pbsprj" >> $qsub1
		    echo "#PBS -l walltime=$pbscpu" >> $qsub1
		    echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
                    echo $ccmgres_string >> $qsub1
                    echo "#PBS -q $pbsqueue" >> $qsub1
                    echo "#PBS -m ae" >> $qsub1
                    echo "#PBS -M $email" >> $qsub1
		    echo "aprun -n 1 -d $thr $scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.R1.sai $AlignOutputDir/$SampleName/$Rone $scriptdir $AlignOutputLogs/log.bwar1.$SampleName.node$i.in $AlignOutputLogs/log.bwar1.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwar1.$SampleName.node$i" >> $qsub1

                    `chmod a+r $qsub1`
                    jobr1=`qsub $qsub1`
                    `qhold -h u $jobr1`
                    echo $jobr1 >> $AlignOutputLogs/ALIGNED_$SampleName
                    if [ $paired -eq 1 ]
                    then
                        echo "bwa aligner. paired-end reads"
			qsub2=$AlignOutputLogs/qsub.bwar2.$SampleName.node$i
			echo "#PBS -V" > $qsub2
			echo "#PBS -N ${pipeid}_bwar2_${SampleName}_$i" >> $qsub2
			echo "#PBS -o $AlignOutputLogs/log.bwar2.$SampleName.node$i.ou" >> $qsub2
			echo "#PBS -e $AlignOutputLogs/log.bwar2.$SampleName.node$i.in" >> $qsub2
			echo "#PBS -A $pbsprj" >> $qsub2
			echo "#PBS -l walltime=$pbscpu" >> $qsub2
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub2
                        echo $ccmgres_string >> $qsub2
			echo "#PBS -q $pbsqueue" >> $qsub2
			echo "#PBS -m ae" >> $qsub2
			echo "#PBS -M $email" >> $qsub2
			echo "$run_string $profiler_string $scriptdir/bwaS1.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.R2.sai $AlignOutputDir/$SampleName/$Rtwo $scriptdir $AlignOutputLogs/log.bwar2.$SampleName.node$i.in $AlignOutputLogs/log.bwar2.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwar2.$SampleName.node$i" >> $qsub2
			`chmod a+r $qsub2`
                        jobr2=`qsub $qsub2`
			`qhold -h u $jobr2`
			echo $jobr2 >> $AlignOutputLogs/ALIGNED_$SampleName

			qsub3=$AlignOutputLogs/qsub.bwasampe.$SampleName.node$i
			echo "#PBS -V" > $qsub3
			echo "#PBS -N ${pipeid}_bwasampe_${SampleName}_$i" >> $qsub3
			echo "#PBS -o $AlignOutputLogs/log.bwasampe.$SampleName.node$i.ou" >> $qsub3
			echo "#PBS -e $AlignOutputLogs/log.bwasampe.$SampleName.node$i.in" >> $qsub3
			echo "#PBS -A $pbsprj" >> $qsub3
			echo "#PBS -l walltime=$pbscpu" >> $qsub3
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub3
                        echo $ccmgres_string >> $qsub3
			echo "#PBS -q $pbsqueue" >> $qsub3
			echo "#PBS -m ae" >> $qsub3
			echo "#PBS -M $email" >> $qsub3
			echo "#PBS -W depend=afterok:$jobr1:$jobr2" >> $qsub3
			echo "$run_string $profiler_string $scriptdir/bwaS2.sh $alignerdir $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.R1.sai $outputsamfileprefix.node$i.R2.sai $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $outputsamfileprefix.node$i.sam $outputsamfileprefix.node$i.bam $samdir $AlignOutputLogs/log.bwasampe.$SampleName.node$i.in $AlignOutputLogs/log.bwasampe.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwasampe.$SampleName.node$i" >> $qsub3
			`chmod a+r $qsub3`
                        jobwa=`qsub $qsub3`
			`qhold -h u $jobwa`
			echo $jobwa >> $AlignOutputLogs/ALIGNED_$SampleName
                    else
                        echo "bwa aligner. single read"
			qsub3=$AlignOutputLogs/qsub.bwasamse.$SampleName.node$i
			echo "#PBS -V" > $qsub3
			echo "#PBS -N ${pipeid}_bwasamse_${SampleName}_$i" >> $qsub3
			echo "#PBS -o $AlignOutputLogs/log.bwasamse.$SampleName.node$i.ou" >> $qsub3
			echo "#PBS -e $AlignOutputLogs/log.bwasamse.$SampleName.node$i.in" >> $qsub3
			echo "#PBS -A $pbsprj" >> $qsub3
			echo "#PBS -l walltime=$pbscpu" >> $qsub3
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub3
                        echo $ccmgres_string >> $qsub3
			echo "#PBS -q $pbsqueue" >> $qsub3
			echo "#PBS -m ae" >> $qsub3
			echo "#PBS -M $email" >> $qsub3
			echo "#PBS -W depend=afterok:$jobr1" >> $qsub3
			echo "aprun -n 1 -d $thr $scriptdir/bwaS3.sh $alignerdir $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.R1.sai $AlignOutputDir/$SampleName/$Rone $outputsamfileprefix.node$i.sam $outputsamfileprefix.node$i.bam $samdir $AlignOutputLogs/log.bwasamse.$SampleName.node$i.in $AlignOutputLogs/log.bwasamse.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwasamse.$SampleName.node$i" >> $qsub3
			`chmod a+r $qsub3`
                        jobwa=`qsub $qsub3`
			`qhold -h u $jobwa`
                        echo $qsub3 >> $AlignOutputLogs/ALIGNED_$SampleName
                    fi
                elif [ $aligner == "BWA_MEM" ]
                then
                    echo "bwa mem is used as aligner. input file format is in fastq"
                    if [ $paired -eq 1 ]
                    then
                        echo "bwa mem aligner. paired-end reads"
                        jobfile=$AlignOutputDir/$SampleName/AnisimovLogs/bwamem.$SampleName.node$i.jobfile

                        if [ $chunkfastq == "YES" ]
                        then
                           echo "$run_string $profiler_string $scriptdir/bwamem_pe.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.sam $outputsamfileprefix.node$i.bam $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $scriptdir $samdir $AlignOutputLogs/log.bwamem.$SampleName.node$i.in $AlignOutputLogs/log.bwamem.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$i" >> $qsub4
                        elif [ $chunkfastq == "NO" ]
                        then
                           echo "nohup $profiler_string $scriptdir/bwamem_pe_markduplicates.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i $AlignOutputDir/$SampleName/$Rone $AlignOutputDir/$SampleName/$Rtwo $runfile $AlignOutputDir/$SampleName/AnisimovLogs/log.bwamem.$SampleName.node$i.in $AlignOutputDir/$SampleName/AnisimovLogs/log.bwamem.$SampleName.node$i.ou $email $jobfile $RGparms $AlignOutputLogs > $AlignOutputDir/$SampleName/AnisimovLogs/log.bwamem.$SampleName.node$i.in" > $jobfile
                           jobfilename=$( basename $jobfile )
                           echo "$AlignOutputDir/$SampleName/AnisimovLogs $jobfilename" >> $AlignOutputLogs/AlignAnisimov.joblist
                        fi


                     else
                        echo "bwa mem aligner. single-end reads"
                        qsub5=$AlignOutputLogs/qsub.bwamem.$SampleName.node$i
                        echo "#PBS -V" > $qsub5
                        echo "#PBS -N ${pipeid}_bwamem_${SampleName}_$i" >> $qsub5
                        echo "#PBS -o $AlignOutputLogs/log.bwamem.$SampleName.node$i.ou" >> $qsub5
                        echo "#PBS -e $AlignOutputLogs/log.bwamem.$SampleName.node$i.in" >> $qsub5
                        echo "#PBS -A $pbsprj" >> $qsub5
                        echo "#PBS -l walltime=$pbscpu" >> $qsub5
                        echo "#PBS -l nodes=1:ppn=$thr" >> $qsub5
                        echo "#PBS -q $pbsqueue" >> $qsub5
                        echo "#PBS -m ae" >> $qsub5
                        echo "#PBS -M $email" >> $qsub5

                        if [ $chunkfastq == "YES" ]
                        then
                           echo "aprun -n 1 -d $thr $scriptdir/bwamem_se.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.sam $outputsamfileprefix.node$i.bam $AlignOutputDir/$SampleName/$Rone $scriptdir $samdir $AlignOutputLogs/log.bwamem.$SampleName.node$i.in $AlignOutputLogs/log.bwamem.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$i" >> $qsub5
                        elif [ $chunkfastq == "NO" ]
                        then
                           echo "aprun -n 1 -d $thr $scriptdir/bwamem_se_markduplicates.sh $alignerdir $alignparms $refdir/$refindexed $AlignOutputDir/$SampleName $outputsamfileprefix.node$i.bam $AlignOutputDir/$SampleName/$Rone $scriptdir $samdir $samblasterdir $picardir $AlignOutputLogs/log.bwamem.$SampleName.node$i.in $AlignOutputLogs/log.bwamem.$SampleName.node$i.ou $email $AlignOutputLogs/qsub.bwamem.$SampleName.node$i $RGparms" >> $qsub5
                        fi
   

                       `chmod a+r $qsub5`
                       jobbwamemse=`qsub $qsub4`
                       `qhold -h u $jobbwamemse`
                       echo $jobbwamemse >> $AlignOutputLogs/ALIGNED_$SampleName
                       echo $jobbwamemse >> $TopOutputLogs/ALIGN_NCSA_jobids


                    fi
                fi

                allfiles=$allfiles" $outputsamfileprefix.node$i.bam"

                (( inputfastqcounter++ ))
		echo `date`
            done
            # end loop over chunks of the current input fastq






###########################   FORM POST-ALIGNMENT QSUBS: MERGINE, SORTING, MARKING DUPLICATES  ###########################################

            cat $AlignOutputLogs/ALIGNED_$SampleName >> $TopOutputLogs/ALIGNEDpbs

            if [ $chunkfastq == "YES" ]
            then
   	       #ALIGNED=$( cat $AlignOutputLogs/ALIGNED_* | sed "s/\.[a-z]*//" | tr "\n" ":" )
	       ALIGNED=$( cat $AlignOutputLogs/ALIGNED_* | sed "s/\..*//" | tr "\n" ":" )

	       listfiles=$( echo $allfiles  | tr " " ":" | sed "s/::/:/g" )
               if [ $sortool == "NOVOSORT" ]
               then
                   echo "merging aligned chunks with novosort"
		   qsub1=$AlignOutputLogs/qsub.sortmerge.novosort.$SampleName
		   echo "#PBS -V" > $qsub1
		   echo "#PBS -A $pbsprj" >> $qsub1
		   echo "#PBS -N ${pipeid}_sortmerge_novosort_$SampleName" >> $qsub1
		   echo "#PBS -l walltime=$pbscpu" >> $qsub1
		   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
		   echo "#PBS -o $AlignOutputLogs/log.sortmerge_novosort.$SampleName.ou" >> $qsub1
		   echo "#PBS -e $AlignOutputLogs/log.sortmerge_novosort.$SampleName.in" >> $qsub1
		   echo "#PBS -q $pbsqueue" >> $qsub1
		   echo "#PBS -m ae" >> $qsub1
		   echo "#PBS -M $email" >> $qsub1
		   echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
		   echo "aprun -n 1 -d $thr $scriptdir/mergenovo.sh $AlignOutputDir/$SampleName $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $AlignOutputLogs/log.sortmerge_novosort.$SampleName.in $AlignOutputLogs/log.sortmerge_novosort.$SampleName.ou $email $AlignOutputLogs/qsub.sortmerge.novosort.$SampleName" >> $qsub1
		   `chmod a+r $qsub1`
		   mergejob=`qsub $qsub1`
		   `qhold -h u $mergejob`
		   echo $mergejob  >> $AlignOutputLogs/MERGED_$SampleName
               else
                   echo "merging aligned chunks with picard"
		   qsub1=$AlignOutputLogs/qsub.sortmerge.picard.$SampleName
		   echo "#PBS -V" > $qsub1
		   echo "#PBS -A $pbsprj" >> $qsub1
		   echo "#PBS -N ${pipeid}_sortmerge_picard_$SampleName" >> $qsub1
		   echo "#PBS -l walltime=$pbscpu" >> $qsub1
		   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
		   echo "#PBS -o $AlignOutputLogs/log.sortmerge.picard.$SampleName.ou" >> $qsub1
		   echo "#PBS -e $AlignOutputLogs/log.sortmerge.picard.$SampleName.in" >> $qsub1
		   echo "#PBS -q $pbsqueue" >> $qsub1
		   echo "#PBS -m ae" >> $qsub1
		   echo "#PBS -M $email" >> $qsub1
		   echo "#PBS -W depend=afterok:$ALIGNED" >> $qsub1
		   echo "aprun -n 1 -d $thr $scriptdir/mergepicard.sh $AlignOutputDir/$SampleName $listfiles $outsortwdup $outsortnodup $sortedplain $dupparms $RGparms $runfile $AlignOutputLogs/log.sortmerge.picard.$SampleName.in $AlignOutputLogs/log.sortmerge.picard.$SampleName.ou $email $AlignOutputLogs/qsub.sortmerge.picard.$SampleName" >> $qsub1
		   `chmod a+r $qsub1`
		   mergejob=`qsub $qsub1`
		   `qhold -h u $mergejob`
		   echo $mergejob  >> $AlignOutputLogs/MERGED_$SampleName
               fi

	       echo `date`
	       echo "extract reads specified in CHRINDEX param"
	       qsub5=$AlignOutputLogs/qsub.extractreadsbam.$SampleName
   	       echo "#PBS -V" > $qsub5
	       echo "#PBS -A $pbsprj" >> $qsub5
	       echo "#PBS -N ${pipeid}_extrbam_$SampleName" >> $qsub5
	       echo "#PBS -l walltime=$pbscpu" >> $qsub5
	       echo "#PBS -l nodes=1:ppn=$thr" >> $qsub5
	       echo "#PBS -o $AlignOutputLogs/log.extractreadsbam.$SampleName.ou" >> $qsub5
	       echo "#PBS -e $AlignOutputLogs/log.extractreadsbam.$SampleName.in" >> $qsub5
	       echo "#PBS -q $pbsqueue" >> $qsub5
	       echo "#PBS -m ae" >> $qsub5
	       echo "#PBS -M $email" >> $qsub5
	       echo "#PBS -W depend=afterok:$mergejob" >> $qsub5
	       echo "aprun -n 1 -d $thr $scriptdir/extract_reads_bam.sh $AlignOutputDir/$SampleName $outsortwdup $runfile $AlignOutputLogs/log.extractreadsbam.$SampleName.in $AlignOutputLogs/log.extractreadsbam.$SampleName.ou $email  $AlignOutputLogs/qsub.extractreadsbam.$SampleName $igv $extradir" >> $qsub5
	       `chmod a+r $qsub5`
	       extrajob=`qsub $qsub5`
               `qhold -h u $extrajob`
               echo $extrajob >> $TopOutputLogs/EXTRACTREADSpbs

	       cat $AlignOutputLogs/MERGED_$SampleName >> $TopOutputLogs/MERGEDpbs
            fi
          fi
	done <  $TopOutputLogs/SAMPLENAMES.list
        # end loop over input fastq






#############################################################################################################
#####################################                                ########################################
#####################################  SCHEDULE QSUBS CREATED ABOVE  ########################################
#####################################                                ########################################
#############################################################################################################


        if [ $chunkfastq == "NO" ]
        then

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
           if [ -z $profiler_string ]
           then
              echo "aprun -n $numalignnodes -N 1 -d $thr ~anisimov/scheduler/scheduler.x $AlignOutputLogs/AlignAnisimov.joblist /bin/bash > $AlignOutputLogs/AlignAnisimov.joblist.log" >> $qsubAlignLauncher
           # otherwise use ccm
           else 
              echo $ccmgres_string >> $qsubAlignLauncher
              echo "$run_string ~anisimov/scheduler/scheduler.x $AlignOutputLogs/AlignAnisimov.joblist /bin/bash > $AlignOutputLogs/AlignAnisimov.joblist.log" >> $qsubAlignLauncher
           fi

           AlignAnisimovJoblistId=`qsub $qsubAlignLauncher`

           echo $AlignAnisimovJoblistId >> $TopOutputLogs/ALIGNEDpbs # so that this job could be released in the next section. Should it be held to begin with?
           echo $AlignAnisimovJoblistId >> $TopOutputLogs/MERGEDpbs # so that summaryok could depend on this job, in case when there is no merging: a sigle chunk
           echo $AlignAnisimovJoblistId >> $TopOutputLogs/ALIGN_NCSA_jobids # so that sort jobs in realignment/recalibration block could find the job ids to depend upon 

        fi




########################################################################################################
#####################################                           ########################################
#####################################  WRAP UP ALIGNMENT BLOCK  ########################################
#####################################                           ########################################
########################################################################################################

        
	pbsids=$( cat $TopOutputLogs/MERGEDpbs | sed "s/\..*//" | tr "\n" ":" )
	extraids=$( cat $TopOutputLogs/EXTRACTREADSpbs | sed "s/\..*//" | tr "\n" " " )
        mergeids=$( echo $pbsids | tr ":" " " )
        alignids=$( cat $TopOutputLogs/ALIGNEDpbs | sed "s/\..*//" | tr "\n" " " )
	echo $pbsids >> $TopOutputLogs/ALIGN_NCSA_jobids

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
		qsub6=$TopOutputLogs/qsub.cleanup.align
		echo "#PBS -V" > $qsub6
		echo "#PBS -A $pbsprj" >> $qsub6
		echo "#PBS -N ${pipeid}_cleanup_aln" >> $qsub6
		echo "#PBS -l walltime=$pbscpu" >> $qsub6
		echo "#PBS -l nodes=1:ppn=1" >> $qsub6
		echo "#PBS -o $TopOutputLogs/log.cleanup.align.ou" >> $qsub6
		echo "#PBS -e $TopOutputLogs/log.cleanup.align.in" >> $qsub6
		echo "#PBS -q $pbsqueue" >> $qsub6
		echo "#PBS -m ae" >> $qsub6
		echo "#PBS -M $email" >> $qsub6
		echo "#PBS -W depend=afterok:$pbsids" >> $qsub6
		echo "aprun -n 1 -d $thr $scriptdir/cleanup.sh $outputdir $analysis $TopOutputLogs/log.cleanup.align.in $TopOutputLogs/log.cleanup.align.ou $email $TopOutputLogs/qsub.cleanup.align"  >> $qsub6
		`chmod a+r $qsub6`
		cleanjobid=`qsub $qsub6`
		echo $cleanjobid >> $outputdir/logs/CLEANUPpbs
            fi

            `slepp 30s`
	    qsub4=$TopOutputLogs/qsub.summary.aln.allok
	    echo "#PBS -V" > $qsub4
	    echo "#PBS -A $pbsprj" >> $qsub4
	    echo "#PBS -N ${pipeid}_summaryok" >> $qsub4
	    echo "#PBS -l walltime=$pbscpu" >> $qsub4
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub4
	    echo "#PBS -o $TopOutputLogs/log.summary.aln.ou" >> $qsub4
	    echo "#PBS -e $TopOutputLogs/log.summary.aln.in" >> $qsub4
	    echo "#PBS -q $pbsqueue" >> $qsub4
	    echo "#PBS -m ae" >> $qsub4
	    echo "#PBS -M $email" >> $qsub4
            if [ `expr ${#cleanjobid}` -gt 0 ]
            then
		echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub4
            else
		echo "#PBS -W depend=afterok:$pbsids" >> $qsub4
            fi
	    echo "aprun -n 1 -d 1 $scriptdir/summary.sh $outputdir $email exitok"  >> $qsub4
	    `chmod a+r $qsub4`
	    lastjobid=`qsub $qsub4`
	    echo $lastjobid >> $TopOutputLogs/SUMMARYpbs

	    if [ `expr ${#lastjobid}` -lt 1 ]
	    then
		echo "at least one job aborted"
		qsub5=$TopOutputLogs/qsub.summary.aln.afterany
		echo "#PBS -V" > $qsub5
		echo "#PBS -A $pbsprj" >> $qsub5
		echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub5
		echo "#PBS -l walltime=$pbscpu" >> $qsub5
		echo "#PBS -l nodes=1:ppn=1" >> $qsub5
		echo "#PBS -o $TopOutputLogs/log.summary.aln.afterany.ou" >> $qsub5
		echo "#PBS -e $TopOutputLogs/log.summary.aln.afterany.in" >> $qsub5
		echo "#PBS -q $pbsqueue" >> $qsub5
		echo "#PBS -m ae" >> $qsub5
		echo "#PBS -M $email" >> $qsub5
		echo "#PBS -W depend=afterany:$pbsids" >> $qsub5
		echo "aprun -n 1 -d $thr $scriptdir/summary.sh $outputdir $email exitnotok"  >> $qsub5
		`chmod a+r $qsub5`
		badjobid=`qsub $qsub5`
		echo $badjobid >> $TopOutputLogs/SUMMARYpbs
	    fi
	fi

	if [ $analysis == "REALIGNMENT" -o $analysis == "REALIGN" ]
	then
            echo " analysis continues with realignment"
	    qsub2=$TopOutputLogs/qsub.main.realn
	    echo "#PBS -V" > $qsub2
	    echo "#PBS -A $pbsprj" >> $qsub2
	    echo "#PBS -N ${pipeid}_MAINrealn" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	    echo "#PBS -o $TopOutputLogs/MAINrealn.ou" >> $qsub2
	    echo "#PBS -e $TopOutputLogs/MAINrealn.in" >> $qsub2
	    echo "#PBS -q $pbsqueue" >> $qsub2
	    echo "#PBS -m ae" >> $qsub2
	    echo "#PBS -M $email" >> $qsub2
            #echo "#PBS -W depend=afterany:$pbsids" >> $qsub2
	    echo "$scriptdir/realign.sh $runfile $TopOutputLogs/MAINrealn.in $TopOutputLogs/MAINrealn.ou $email $TopOutputLogs/qsub.main.realn" >> $qsub2
	    `chmod a+r $qsub2` 
	    `qsub $qsub2 >> $TopOutputLogs/MAINREALNpbs`

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
