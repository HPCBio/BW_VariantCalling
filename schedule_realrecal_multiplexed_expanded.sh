#!/bin/bash
# written in collaboration with Mayo Bioinformatics core group
# and H3A Africa
#
# script to produce improved bams - one per sample 
# from  multiplexed aligned bams 
#redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 7 ]
then
   MSG="parameter mismatch."
   echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
   exit 1;
fi

   set -x
   echo `date`
   scriptfile=$0
   outputdir=$1
   runfile=$2
   flag=$3
   elog=$4
   olog=$5
   email=$6
   qsubfile=$7
   LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

   if [ ! -s $runfile ]
   then
      MSG="$runfile configuration file not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
   fi

# wrapping commends in echo, so that the output logs would be easier to read: they will have more structure

   echo -e "####################################################################################################
   #####################################                       ######################################## 
   ##################################### PARSING RUN INFO FILE ########################################
   #####################################  SANITY CHECK         ######################################## 
   ####################################################################################################" 

   reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
   pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
   thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
   input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   run_method=$( cat $runfile | grep -w RUNMETHOD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
   scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
   refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
   ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
   picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
   samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
   gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
   tabixdir=$( cat $runfile | grep -w TABIXDIR | cut -d '=' -f2 )
   vcftoolsdir=$( cat $runfile | grep -w VCFTOOLSDIR | cut -d '=' -f2 )
   dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
   #kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
   targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
   indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
   realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
   multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
   chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
   indices=$( echo $chrindex | sed 's/:/ /g' )
   sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
   sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
   javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
   skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
   cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   
   echo "#################################################################################"
   echo -e cheching type of analysis
   echo "#################################################################################"

   if [ $analysis != "MULTIPLEXED" ]
   then
      MSG="ANALYSIS=$analysis Program=$scriptfile Invalid pipeline program for this type of analysis. This program is for the MULTIPLEXED case only"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi

   echo "#################################################################################"
   echo -e type of sample
   echo "#################################################################################"

   if [ $multisample != "YES" -a $multisample != "1" ]
   then
      MSG="MULTISAMPLE=$multisample Invalid value for this type of analysis. This program is for the MULTIPLEXED case only"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi

   echo "#################################################################################"
   echo -e checking input type
   echo "#################################################################################"

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
         echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         exit 1;
      fi
   fi


   echo "#################################################################################"
   echo -e checking cleanup options
   echo "#################################################################################"

   if [ $cleanupflag != "1" -a $cleanupflag != "0" -a $cleanupflag != "YES" -a $cleanupflag != "NO" ]
   then
      MSG="Invalid value for REMOVETEMPFILES=$cleanupflag"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
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

   echo "#################################################################################"
   echo -e skip/include variant calling module
   echo "#################################################################################"

   if [ $skipvcall != "1" -a $skipvcall != "0" -a $skipvcall != "YES" -a $skipvcall != "NO" ]
   then
      MSG="Invalid value for SKIPVCALL=$skipvcall"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   else
      if [ $skipvcall == "1" ]
      then
         $skipvcall="YES"
      fi
      if [ $skipvcall == "0" ]
      then
         $skipvcall="NO"
      fi
   fi

   echo "#################################################################################"
   echo -e java for gatk
   echo "#################################################################################"

   if [ -z $javamodule ]
   then
      MSG="Value for JAVAMODULE must be specified in configuration file"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

   echo "#################################################################################"
   echo -e directories for output and for tools
   echo "#################################################################################"

   if [ ! -d $outputdir ]
   then
      MSG="$outputdir ROOT directory for this run of the pipeline not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

   if [ ! -d $picardir ]
   then
      MSG="$picardir picard directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi
   if [ ! -d $samdir ]
   then
      MSG="$samdir samtools directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi
   if [ ! -d $gatk ]
   then
      MSG="$gatk GATK directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

   if [ ! -d $refdir ]
   then
      MSG="$refdir reference genome directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi
   if [ ! -s $refdir/$ref ]
   then
      MSG="$ref reference genome not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi
  if [ ! -s $refdir/$dbSNP ]
   then
      MSG="$refdir/$dbSNP dbSNP for reference genome not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi
   if [ ! -d $refdir/$indeldir ]
   then
      MSG="$indeldir indel directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi
   if [ ! -s $sampleinfo ]
   then
      MSG="$sampleinfo SAMPLEINFORMATION file not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

#   if [ -s $refdir/$dbSNP ]
#   then
#      realparms="-known:$refdir/$dbSNP"
#      recalparms="--knownSites:$refdir/$dbSNP"
#   else
#      MSG="dbSNP not found"
#      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
#      exit 1;
#   fi
#   if [ -s $refdir/$kgenome ]
#   then
#      realparms=$realparms":-known:$refdir/$kgenome"
#      recalparms=$recalparms":--knownSites:$refdir/$kgenome"
#   fi


   echo "#################################################################################"
   echo -e checking scheduler
   echo "#################################################################################"

   if [ $run_method != "LAUNCHER" -a $run_method != "QSUB" -a $run_method != "APRUN" -a $run_method != "SERVER" ]
   then
      MSG="Invalid value for RUNMETHOD=$run_method"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi

   echo "#################################################################################"
   echo -e checking files produced from info sheet
   echo "#################################################################################"

   if [ ! -s $outputdir/SAMPLENAMES.list ]
   then
      MSG="$outputdir/SAMPLENAMES.list SAMPLENAMES.list file not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      exit 1;
   fi

   if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
   then
      MSG="$outputdir/SAMPLENAMES_multiplexed.list SAMPLENAMES_multiplexed.list file not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      exit 1;
   fi

   if [ ! -s $outputdir/SAMPLEGROUPS.list ]
   then
      MSG="$outputdir/SAMPLEGROUPS.list SAMPLEGROUPS.list file not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      exit 1;
   fi

   echo "####################################################################################################"
   echo "####################################################################################################"
   echo "###########################      params ok. Now          creating log folders        ###############"
   echo "####################################################################################################"
   echo "####################################################################################################"


   TopOutputLogs=$outputdir/logs
   RealignOutputLogs=$outputdir/logs/realign
   VcallOutputLogs=$outputdir/logs/variant
   if [ ! -d $TopOutputLogs ]
   then
      MSG="$TopOutputLogs realign directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      exit 1;
   fi
   if [ ! -d $RealignOutputLogs ]
   then
      mkdir -p $RealignOutputLogs
   fi
   if [ ! -d $VcallOutputLogs ]
   then
      mkdir -p $VcallOutputLogs
   fi
   pipeid=$( cat $TopOutputLogs/CONFIGUREpbs )

   if [ $schedule == "LAUNCHER" ]
   then
      truncate -s 0 $RealignOutputLogs/realrecal.AnisimovJoblist
   fi

   echo "####################################################################################################"
   echo "####################################################################################################"
   echo "###########################   generating regions, intervals, known/knownSites        ###############"
   echo "####################################################################################################"
   echo "####################################################################################################"

   echo `date`
   i=1
   echo -e "region variable, to be used by vcallgatk-unifiedGenotyper"
   region=$refdir/$dbSNP

   echo -e "checking that there are indels to use for the realign-recalibration analysis "
   echo -e "realparms will have indels for  gatk-IndelRealigner"
   echo -e "recalparms will have indels for  gatk-Recalibration"
   for chr in $indices
   do
      #echo "this chunk with region setup does not apply to H3A dataset"
      #if [ -d $targetkit ]
      #then
      #   if [ `cat $targetkit/${chr}.bed | wc -l` -gt 0 ]
      #   then
      #      region[$i]="-L:$targetkit/${chr}.bed"
      #   else
      #      region[$i]="-L:$chr"
      #   fi
      #else
      #   region[$i]="-L:$chr"
      #fi
      cd $refdir/$indeldir
      indels=`find $PWD -type f -name "${chr}.*.vcf"`
      realparms[$i]=$( echo $indels | sed "s/\/projects/:known:\/projects/g" | sed "s/ //g" |tr "\n" ":" )
      recalparms[$i]=$( echo $indels | sed "s/\/projects/:knownSites:\/projects/g" | sed "s/ //g" | tr "\n" ":" )
      (( i++ ))
   done
   echo `date`


   echo "####################################################################################################"
   echo "####################################################################################################"
   echo "#####################################  CREATE  DIRECTORIES   #######################################"
   echo "#################################  4 RESULTS BY SAMPLE/LANE      ####################################"
   echo "####################################################################################################"
   echo "####################################################################################################"



   echo -e "the new directory structure needs to be created here:"
   echo -e "FOR EACH SAMPLE"
   echo -e "outputdir/sample/realign   outputdir/sample/variant   outputdir/sample/logs"
   echo -e "FOR EACH LANE"
   echo -e "outputdir/sample/lane/align outputdir/sample/lane/realign  outputdir/sample/lane/logs"
   echo -e "soft links are created to alignment results obtained in previous steps\n\n\n"

   while read SampleLine
   do
      if [ `expr ${#SampleLine}` -gt 1 ]
      then
          ## parsing non-empty line
	  sample=$( echo "$SampleLine" | cut -f 1 )
	  lane=$( echo "$SampleLine" | cut -f 3 )
          echo -e "################################################################################"
          echo -e "################################################################################"
	  echo -e "########## first, let's checking that alignment info exists  ###################"
          echo -e "################################################################################"

	  alignedfile=`find $outputdir/align/$lane/ -name "*.wdups.sorted.bam"`
	  alignedfilehdr=`find $outputdir/align/$lane/ -name "*.wdups.sorted.bam.header"`
	  alignedfilebai=`find $outputdir/align/$lane/ -name "*.wdups.sorted.bam.bai"`
	  if [ -s $alignedfile -a -s $alignedfilehdr -a -s $alignedfilebai ]
	  then
              echo -e "alignment files for this lane $lane were found"
          else
              MSG="No aligned bam file(s) found at $outputdir/align/${lane}"
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
              exit 1;              
	  fi
          echo -e "################################################################################"
          echo -e "################################################################################"
          echo -e "########## now we can create the folders for the rest of the analysis  #########"
          echo -e "################################################################################"

	  if [ ! -d $outputdir/${sample} ]
	  then
              echo -e "first time we see this sample $sample"
              echo -e "creating output folders for it"
              mkdir -p $outputdir/${sample}
              mkdir -p $outputdir/${sample}/realign/logs
              mkdir -p $outputdir/${sample}/variant/logs
              mkdir -p $outputdir/${sample}/logs
              mkdir -p $outputdir/${sample}/$lane
              ln -s    $outputdir/align/$lane $outputdir/${sample}/$lane/align
              mkdir -p $outputdir/${sample}/${lane}/realign/logs
              mkdir -p $outputdir/${sample}/${lane}/logs
	  elif [ ! -d $outputdir/${sample}/${lane} ]
	  then
              echo -e "NOT the first time we see this sample $sample"
              echo -e "populating $sample with folders for lane $lane"
              mkdir -p $outputdir/${sample}/$lane
              ln -s $outputdir/align/$lane $outputdir/${sample}/$lane/align
              mkdir -p $outputdir/${sample}/${lane}/realign/logs
              mkdir -p $outputdir/${sample}/${lane}/logs
	  else
             MSG="lane $lane for sample $samplealready exists. Error. Lane name has to be unique"
	     echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
             exit 1;
          fi
      fi
   done  <  $outputdir/SAMPLENAMES_multiplexed.list
   # end loop over sampleLanes 



   echo -e "#############################################################################################################
   ###################################                                #############################################
   ###################################   nested loop1      starts    here    ###################################### 
   ########################         call to realig-recalibrate per lane per chr      ##############################
   ########################   outer loop by Lane; inner loop by chromosome           ##############################
   ########################                                                      ##################################
   ################################################################################################################"


lanecounter=1 
while read SampleLine
do
    if [ `expr ${#SampleLine}` -gt 1 ]
    then
        echo -e "################################################################################"
        echo -e "################################################################################"
	echo -e "######  processing next non-empty line in SAMPLENAMES_multiplexed.list #########"
        echo -e "################################################################################"
        echo -e "################################################################################"
        # now parsing the line just being read
	sample=$( echo "$SampleLine" | cut -f 1 )
	lane=$( echo "$SampleLine" | cut -f 3 )

        echo -e "################################################################################"
	echo -e "######    realigning recalibrating per lane: $lane sample: $sample   ###########"
        echo -e "######    gathering input files and other parameters                 ###########"
        echo -e "################################################################################"

	echo `date`
        RealignOutputDir=$outputdir/${sample}/${lane}/realign
        AlignOutputDir=$outputdir/${sample}/${lane}/align
        RealignLog=$outputdir/${sample}/realign/logs
        cd $AlignOutputDir
        input_bam=`find ./ -name "*.wdups.sorted.bam"`
	input_bam=$( echo $input_bam | sed "s/\.\///g" | tr "\n" " " | sed "s/ //g" )
        inputfile=$AlignOutputDir/$input_bam
        RGline=${input_bam}.RGline

        if [ `expr ${#input_bam}` -lt 1 -o ! -s $input_bam ]
        then
            MSG="No bam file(s) found to perform realign-recalibrate at $outputdir/${lane}/align"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            exit 1;
        fi
        echo -e "################################################################################"      
        echo -e  "now putting together the RG line for RGparms"
        echo -e "################################################################################"

        if [ `expr ${#RGline}` -lt 1 -o ! -s $RGline ]
        then
           echo -e "RGparms line needs to be recreated from scratch" 
           RGparms=$( grep "^@RG" *wdups*header | sed 's/@RG//' | tr ":" "=" | tr " " ":" | tr "\t" ":" | sed "s/ //g" )
        else 
           RGparms=$( cat $RGline )
        fi

        echo -e "################################################################################"      
        echo -e  "now about to enter the loop on chromosomes "
        echo -e "################################################################################"

	chromosomecounter=1
	for chr in $indices
	do
            echo -e "################################################################################"
	    echo "generating real-recal calls for lane=$lane chr=${chr} ..."
            echo -e "################################################################################"

            echo "$scriptdir/realrecalLane.sh $lane $RealignOutputDir ${chr}.realrecal.${lane}.output.bam $chr $inputfile $RGparms ${realparms[$chr]} ${recalparms[$chr]} $runfile $RealignLog/log.realrecalLane.$lane.$chr.in $RealignLog/log.realrecalLane.$lane.$chr.ou $email $RealignLog/realrecalLane.${lane}.${chr}" > $RealignLog/realrecalLane.${lane}.${chr}
           

         # end loop over chromosomes
	 echo `date`
	 ((  chromosomecounter++ )) 
	done # done going through chromosomes 
    fi   # skipping empty lines

    # end loop over sampleLanes
    (( lanecounter++ ))
done <  $outputdir/SAMPLENAMES_multiplexed.list



echo -e "\n\n\n #############################################################################################################
                ###################################                                ##########################################
                ###################################       loop2      starts    here    ###################################### 
                ##################################  MERGE1: COMBINE LANES INTO SAMPLES     ##################################
                ########################       MERGE2: COMBINE CHR INTO SAMPLE and VERIFYBAM    #############################
                #############################################################################################################"

groupcounter=1 
while read SampleLine
do
    if [ `expr ${#SampleLine}` -gt 1 ]
    then
	echo "processing next sample, non-empty line"
	echo `date`

	sample=$( echo "$SampleLine" | cut -f 1 )
	lanes=$( echo "$SampleLine" | cut -f 2 )
        lanes=$( echo $lanes | tr " " ":" )

        echo -e "#############################################################################################################"
	echo -e "  setting up params for the commands"
        echo -e "#############################################################################################################"

        RealignLog=$outputdir/${sample}/realign/logs
        RealignOutputDir=$outputdir/${sample}/realign
        Merged1OutPrefix=$outputdir/${sample}/realign/${sample}.mergedLanes
        Merged2OutPrefix=$outputdir/${sample}/realign/${sample}.mergedVerified
        sLB="multiplexed"

        echo -e "#############################################################################################################"
        echo -e "#######################   MERGE1 assemble the mergeLanes call for sample=$sample ############################\n"
        echo -e "#############################################################################################################"

        echo "$scriptdir/mergeLanesPerSample.sh $runfile $sample $sLB $lanes ${Merged1OutPrefix}.bam $outputdir $RealignLog/log.mergeLanes.$sample.in $RealignLog/log.mergeLanes.$sample.ou $email $RealignLog/mergeLanes.${sample}" > $RealignLog/mergeLanes.${sample}


        echo -e "#############################################################################################################"
        echo -e "#######################   MERGE2 assemble the mergeAndVerify call for sample=$sample    #####################\n"
        echo -e "#############################################################################################################"

        echo "$scriptdir/verifySample.sh $runfile $sample $RealignOutputDir ${Merged2OutPrefix}.bam $outputdir $RealignLog/log.mergeAndVerify.$sample.in $RealignLog/log.mergeAndVerify.$sample.ou $email $RealignLog/mergeAndVerify.${sample}" > $RealignLog/mergeAndVerify.${sample}

    fi   # skipping empty lines

    # end loop over sampleLanes
    (( groupcounter++ ))
done <  $outputdir/SAMPLEGROUPS.list


echo -e "\n\n\n #############################################################################################################
                ###################################                                ##########################################
                ###################################       loop3      starts    here    ###################################### 
                ##################################       second REALIGN  -- BY SAMPLE              ##########################
                ##################################       followed by vcall -- BY SAMPLE            ##########################
                ##################################  outer loop by Sample inner loop by chr         ##########################
                ########################                                                      ###############################
                #############################################################################################################"

samplecounter=1 
while read SampleLine
do
    if [ `expr ${#SampleLine}` -gt 1 ]
    then
        echo -e "################################################################################"
	echo "processing next sample, non-empty line"
        echo -e "################################################################################"


	sample=$( echo "$SampleLine" | cut -f 1 )
        RealignLog=$outputdir/${sample}/realign/logs	
        RealignOutputDir=$outputdir/${sample}/realign
        sampleOutPrefix=${sample}.mergedLanes
        bamfile=$RealignOutputDir/$sampleOutPrefix.bam
        VcallOutputDir=$outputdir/$sample/variant

        echo "sample=$sample about to enter inner loop by chr..."
	chromosomecounter=1
	for chr in $indices
	do
            outfile=${chr}.$sampleOutPrefix.realigned.calmd.bam
            echo "#############################################################################################################################"
	    echo "##############   generating second realign call for sample=$sample chr=${chr} with bam=$infile ##############################"
            echo "#############################################################################################################################"

            echo "$scriptdir/realignSample.sh $RealignOutputDir $outfile $chr $bamfile  ${realparms[$chr]} $sample $runfile $RealignLog/log.realignSample.$sample.$chr.in $RealignLog/log.realignSample.$sample.$chr.ou $email $RealignLog/realignSample.${sample}.${chr}" > $RealignLog/realignSample.${sample}.${chr}
           
           if [ $skipvcall == "NO" ]
           then
               echo "#############################################################################################################################"
               echo "##############  generating  vcallgatk call for sample=$sample chr=${chr} with bam=$infile ###################################"
               echo "#############################################################################################################################"

	       echo "$scriptdir/vcallgatk.sh $VcallOutputDir $RealignOutputDir  $outfile $chr $region $runfile $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.ou $email $VcallOutputDir/logs/vcallgatk.${sample}.${chr}" >> $VcallOutputDir/logs/vcallgatk.${sample}.${chr}
	   fi

         # end loop over chromosomes
	 echo `date`
	 ((  chromosomecounter++ )) 
	done # done going through chromosomes 
    fi   # skipping empty lines

    # end loop over sampleLanes
    (( samplecounter++ ))
done <  $outputdir/SAMPLEGROUPS.list


echo -e "\n\n   #############################################################################################################
                ###################################                             #############################################
                ###################################   now schedule these jobs   #############################################
                ###################################   just considering the launcher option for now ##########################
                #############################################################################################################"

   case $run_method in
   "LAUNCHER")
      echo -e "#################################################################################################"
      echo -e "LAUNCHER. This case is now been tested. Proceed with caution"
      echo -e "gathering jobs generated in loop1"

      echo -e "moved the following code from the loop to outside the loop "
      echo -e "no need to repeat this calculation "
      echo -e "#################################################################################################"

 
      echo -e "#################################################################################################
               #########     colecting jobs from loop 1 realign-recalibrate by lane then by chr ################
               #################################################################################################"

      RealignOutputLogs=$outputdir/logs/realign
      VcallOutputLogs=$outputdir/logs/variant


      while read SampleLine
      do
	  if [ `expr ${#SampleLine}` -gt 1 ]
	  then
	      echo "#################################################################################################"
	      echo "#############   processing next non-empty line in SAMPLENAMES_multiplexed.list.##################"
	      echo "#############   realign-recalibrate by lane then by chromosome                 ##################"
	      echo "#################################################################################################"
              # now parsing the line just being read
	      sample=$( echo "$SampleLine" | cut -f 1 )
	      lane=$( echo "$SampleLine" | cut -f 3 )

              RealignOutputDir=$outputdir/$sample/$lane/realign
              truncate -s 0 $RealignOutputLogs/realrecalLane.${lane}.AnisimovJoblist
     
	      echo "#################################################################################################"
              echo "##########   collecting all jobs for sample=$sample lane=$lane on all chromosomes ###############"
	      echo "#############   realign-recalibrate by lane then by chromosome                 ##################"
	      echo "#################################################################################################"

	      for chr in $indices
	      do
                   # creating a qsub out of the job file
                   # need to prepend "nohup" and append log file name, so that logs are properly created when Anisimov launches these jobs 

		  realrecal_log=$outputdir/$sample/realign/logs/log.realrecalLane.$lane.$chr.in
		  awk -v awkvar_realrecallog=$realrecal_log '{print "nohup "$0" > "awkvar_realrecallog}' $outputdir/$sample/realign/logs/realrecalLane.${lane}.${chr} > $outputdir/$sample/realign/logs/jobfile.realrecalLane.${lane}.${chr}
		  echo "$outputdir/$sample/realign/logs/ jobfile.realrecalLane.${lane}.${chr}" >> $RealignOutputLogs/realrecalLane.${lane}.AnisimovJoblist

              done
              # end loop over chr

	      echo "#################################################################################################"
	      echo "#############   realign-recalibrate by lane then by chromosome                 ##################"
              echo -e "############ generating qsub for the corresponding  Anisimov Launcher joblists ###############"
	      echo "#################################################################################################"

              qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecalLane.${lane}.AnisimovLauncher

              cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov

              echo "#PBS -N ${pipeid}_realrecalxlane_${lane}" >> $qsub_realrecal_anisimov
              echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal_anisimov
              echo "#PBS -o $RealignOutputLogs/log.realrecalLane.${lane}.ou" >> $qsub_realrecal_anisimov
              echo -e "#PBS -e $RealignOutputLogs/log.realrecalLane.${lane}.in\n" >> $qsub_realrecal_anisimov

              echo -e "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_realrecal_anisimov
              echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realrecalLane.${lane}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realrecalLane.${lane}.AnisimovLauncher.log" >> $qsub_realrecal_anisimov

          fi
          #end loop over lanes
      done <  $outputdir/SAMPLENAMES_multiplexed.list

      echo -e "#################################################################################################
      ######### colecting jobs from loop2 and loop3 mergeLanes mergeAndVerify RealignSample Vcallgatk      #####
      ##########################################################################################################"

      RealignOutputLogs=$outputdir/logs/realign
      VcallOutputLogs=$outputdir/logs/variant

      while read SampleLine
      do
	  if [ `expr ${#SampleLine}` -gt 1 ]
	  then
	      echo "#################################################################################################"
	      echo "#############     processing next sample, non-empty line     ####################################"
	      echo "#################################################################################################"
              echo "#############     putting together the AnisimovJobLists      ####################################"
	      echo "#################################################################################################"
              echo "#############     constructing the qsub for mergeLanes       ####################################"
	      echo "#################################################################################################"
	      echo "#################################################################################################"

	      sample=$( echo "$SampleLine" | cut -f 1 )

              RealignOutputDir=$outputdir/$sample/realign
              truncate -s 0 $RealignOutputLogs/mergeLanes.${sample}.AnisimovJoblist
              truncate -s 0 $RealignOutputLogs/mergeAndVerify.${sample}.AnisimovJoblist
              truncate -s 0 $RealignOutputLogs/realignSample.${sample}.AnisimovJoblist
     
	      if [ $skipvcall == "NO" ]
	      then
		  VcallOutputDir=$outputdir/$sample/variant
		  truncate -s 0 $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist
              fi

              echo -e "#########################################################################
              ############### mergelanes per sample   ##########################################
              ##################################################################################"


	      mergeLanes_log=$outputdir/$sample/realign/logs/log.mergelLanes.$sample.in
	      awk -v awkvar_mergeLanes_log=$mergeLanes_log '{print "nohup "$0" > "awkvar_mergeLanes_log}' $outputdir/$sample/realign/logs/mergeLanes.${sample} > $outputdir/$sample/realign/logs/jobfile.mergeLanes.${sample}
	      echo "$outputdir/$sample/realign/logs/ jobfile.mergeLanes.${sample}" >> $RealignOutputLogs/mergeLanes.${sample}.AnisimovJoblist
              qsub_mergeLanes_anisimov=$RealignOutputLogs/qsub.mergeLanes.${sample}.AnisimovLauncher
              cat $outputdir/qsubGenericHeader > $qsub_mergeLanes_anisimov


	      echo "#PBS -N ${pipeid}_mergeLanes_${sample}" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -l walltime=$pbscpu" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -o $RealignOutputLogs/log.mergeLanes.${sample}.ou" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -e $RealignOutputLogs/log.mergeLanes.${sample}.in" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_mergeLanes_anisimov
	      echo "aprun -n 1 -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/mergeLanes.${sample}.AnisimovJoblist /bin/bash > $RealignOutputLogs/mergeLanes.${sample}.AnisimovLauncher.log" >> $qsub_mergeLanes_anisimov



              echo -e "#########################################################################
              ############### mergeAnVerify per sample   ##########################################
              ##################################################################################"


	      mergeAndVerify_log=$outputdir/$sample/realign/logs/log.mergeAndVerify.$sample.in
	      awk -v awkvar_mergeAndVerify_log=$mergeAndVerify_log '{print "nohup "$0" > "awkvar_mergeAndVerify_log}' $outputdir/$sample/realign/logs/mergeAndVerify.${sample} > $outputdir/$sample/realign/logs/jobfile.mergeAndVerify.${sample}
	      echo "$outputdir/$sample/realign/logs/ jobfile.mergeAndVerify.${sample}" >> $RealignOutputLogs/mergeAndVerify.${sample}.AnisimovJoblist
              qsub_mergeAndVerify_anisimov=$RealignOutputLogs/qsub.mergeAndVerify.${sample}.AnisimovLauncher
              cat $outputdir/qsubGenericHeader > $qsub_mergeAndVerify_anisimov

	      echo "#PBS -N ${pipeid}_mergeAndVerify_${sample}" >> $qsub_mergeAndVerify_anisimov
	      echo "#PBS -l walltime=$pbscpu" >> $qsub_mergeAndVerify_anisimov
	      echo -e "#PBS -o $RealignOutputLogs/log.mergeAndVerify.${sample}.ou" >> $qsub_mergeAndVerify_anisimov
	      echo -e "#PBS -e $RealignOutputLogs/log.mergeAndVerify.${sample}.in" >> $qsub_mergeAndVerify_anisimov
	      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_mergeAndVerify_anisimov
	      echo "aprun -n 1 -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/mergeAndVerify.${sample}.AnisimovJoblist /bin/bash > $RealignOutputLogs/mergeAndVerify.${sample}.AnisimovLauncher.log" >> $qsub_mergeAndVerify_anisimov


              echo -e "###########################################################################
              ############### next loop on chromosomes will be used for  #########################
              ############### all realignSample  jobs on chromosomes     #########################
              ############### and all vcallgatk   jobs on chromosomes    #########################
              ####################################################################################"

	      for chr in $indices
	      do
		  echo -e "###################################################################### 
	      #############   sample=$sample  chr=$chr   ########################################
              ###################################################################################\n"
  
		  realignSample_log=$outputdir/$sample/realign/logs/log.realignSample.$sample.$chr.in
		  awk -v awkvar_realignSample=$realignSample_log '{print "nohup "$0" > "awkvar_realignSample}' $outputdir/$sample/realign/logs/realignSample.${sample}.${chr} > $outputdir/$sample/realign/logs/jobfile.realignSample.${sample}.${chr}
		  echo "$outputdir/$sample/realign/logs/ jobfile.realignSample.${sample}.${chr}" >> $RealignOutputLogs/realignSample.${sample}.AnisimovJoblist

                  if [ $skipvcall == "NO" ]
		  then
		      vcall_log=$outputdir/$sample/variant/logs/log.vcallgatk.$sample.$chr.in
		      awk -v awkvar_vcallog=$vcall_log '{print "nohup "$0" > "awkvar_vcallog}' $outputdir/$sample/variant/logs/vcallgatk.${sample}.${chr} > $outputdir/$sample/variant/logs/jobfile.vcallgatk.${sample}.${chr}
		      echo "$outputdir/$sample/variant/logs/ jobfile.vcallgatk.${sample}.${chr}" >> $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist
                  fi
	      done

              echo -e "########################################################################
	      #############   constructing the qsub call for realignSample ####################
              #################################################################################\n"  

	      qsub_realignSample_anisimov=$RealignOutputLogs/qsub.realignSample.${sample}.AnisimovLauncher
              cat $outputdir/qsubGenericHeader > $qsub_realignSample_anisimov

	      echo "#PBS -N ${pipeid}_realignSample_${sample}" >> $qsub_realignSample_anisimov
	      echo "#PBS -l walltime=$pbscpu" >> $qsub_realignSample_anisimov
	      echo -e "#PBS -o $RealignOutputLogs/log.realignSample.${sample}.ou" >> $qsub_realignSample_anisimov
	      echo -e "#PBS -e $RealignOutputLogs/log.realignSample.${sample}.in\n" >> $qsub_realignSample_anisimov
	      echo "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_realignSample_anisimov
	      echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realignSample.${sample}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realignSample.${sample}.AnisimovLauncher.log" >> $qsub_realignSample_anisimov

	      if [ $skipvcall == "NO" ]
	      then

		  echo -e "########################################################################
	      #############   constructing the qsub call for vcallgatk     ####################
              #################################################################################\n"  

                  qsub_vcallgatk_anisimov=$VcallOutputLogs/qsub.vcallgatk.${sample}.AnisimovLauncher
		  cat $outputdir/qsubGenericHeader > $qsub_vcallgatk_anisimov
		  echo "#PBS -N ${pipeid}_vcallgatk_${sample}" >> $qsub_vcallgatk_anisimov
		  echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk_anisimov
		  echo "#PBS -o $VcallOutputLogs/log.vcallgatk.${sample}.ou" >> $qsub_vcallgatk_anisimov
		  echo -e "#PBS -e $VcallOutputLogs/log.vcallgatk.${sample}.in\n" >> $qsub_vcallgatk_anisimov
		  echo -e "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_vcallgatk_anisimov
		  echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist /bin/bash > $VcallOutputLogs/vcallgatk.${sample}.AnisimovLauncher.log" >> $qsub_vcallgatk_anisimov
	      fi
	  fi
      done <  $outputdir/SAMPLEGROUPS.list



echo -e "#############################################################################################
#############################################################################################
#############################################################################################
####            arranging execution order with job dependencies #############################
#############################################################################################
#############################################################################################"

      RealignOutputLogs=$outputdir/logs/realign
      VcallOutputLogs=$outputdir/logs/variant
      
      truncate -s 0 $RealignOutputLogs/REALXSAMPLEXCHRpbs
      truncate -s 0 $RealignOutputLogs/REALRECALXLANEXCHRpbs
      truncate -s 0 $RealignOutputLogs/MERGEXSAMPLEpbs
      truncate -s 0 $RealignOutputLogs/MERGEANDVERIFYXSAMPLEpbs

      if [ $skipvcall == "NO" ]
      then
          truncate -s 0 $VcallOutputLogs/VCALLGATKpbs
      fi
      
      cd $RealignOutputLogs 
      # first we schedule the realrecalxlane
      # FIX: gather all jobids for a single sample and then create the PBS depend line for sample
      
      while read SampleLine
      do
	  if [ `expr ${#SampleLine}` -gt 1 ]
	  then
	      echo "processing next non-empty line in SAMPLENAMES_multiplexed.list."
              # now parsing the line just being read
	      sample=$( echo "$SampleLine" | cut -f 1 )
	      lanes=$( echo "$SampleLine" | cut -f 2 )
	      truncate -s 0 $RealignOutputLogs/${sample}_REALRECALpbs
	      for lane in $lanes
	      do
              	realrecalXlane_job=`qsub $RealignOutputLogs/qsub.realrecalLane.${lane}.AnisimovLauncher`
              	`qhold -h u $realrecalXlane_job` 
              	echo $realrecalXlane_job >> $RealignOutputLogs/REALRECALXLANEXCHRpbs 
              	echo $realrecalXlane_job >> $RealignOutputLogs/${sample}_REALRECALpbs              	
              
              done
              realrecalxsamplexlane=$( cat $RealignOutputLogs/${sample}_REALRECALpbs | sed "s/\..*//" | tr "\n" ":" )
              sed -i "2i #PBS -W depend=afterok:$realrecalxsamplexlane" $RealignOutputLogs/qsub.mergeLanes.${sample}.AnisimovLauncher
     
              mergeSample_job=`qsub $RealignOutputLogs/qsub.mergeLanes.${sample}.AnisimovLauncher`
              `qhold -h u $mergeSample_job` 
              sed -i "2i #PBS -W depend=afterok:$mergeSample_job" $RealignOutputLogs/qsub.realignSample.${sample}.AnisimovLauncher    

              RealignSample_job=`qsub $RealignOutputLogs/qsub.realignSample.${sample}.AnisimovLauncher`
              `qhold -h u $RealignSample_job`
              sed -i "2i #PBS -W depend=afterok:$RealignSample_job" $RealignOutputLogs/qsub.mergeAndVerify.${sample}.AnisimovLauncher

              mergeAndVerify_job=`qsub $RealignOutputLogs/qsub.mergeAndVerify.${sample}.AnisimovLauncher`
              `qhold -h u $mergeAndVerify`
              sed -i "2i #PBS -W depend=afterok:$mergeAndVerify_job" $VcallOutputLogs/qsub.vcallgatk.${sample}.AnisimovLauncher

              VcallSample_job=`qsub $VcallOutputLogs/qsub.vcallgatk.${sample}.AnisimovLauncher`  
              `qhold -h u $VcallSample_job`

              echo $mergeSample_job     >> $RealignOutputLogs/MERGEXSAMPLEpbs
              echo $RealignSample_job   >> $RealignOutputLogs/REALXSAMPLEXCHRpbs  
              echo $VcallSample_job     >> $VcallOutputLogs/VCALLGATKpbs               
              echo $mergeAndVerify_job  >> $RealignOutputLogs/MERGEANDVERIFYXSAMPLEpbs
          fi
      done <  $outputdir/SAMPLEGROUPS.list
   ;;
   esac

echo -e "#############################################################################################
#############################################################################################
#############################################################################################
###########################  wrap up and produce summary table   ############################
#############################################################################################
#############################################################################################"

   if [ $skipvcall == "NO" ]
   then
      summarydependids=$( cat $VcallOutputLogs/VCALLGATKpbs | sed "s/\..*//" | tr "\n" ":" )
   else
      summarydependids=$( cat $RealignOutputLogs/REALXSAMPLEXCHRpbs | sed "s/\..*//" | tr "\n" ":" )
   fi

   lastjobid=""
   cleanjobid=""

   if [ $cleanupflag == "YES" ]
   then
       qsub_cleanup=$TopOutputLogs/qsub.cleanup
       echo "#PBS -V" > $qsub_cleanup
       echo "#PBS -A $pbsprj" >> $qsub_cleanup
       echo "#PBS -N ${pipeid}_cleanup" >> $qsub_cleanup
       echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
       echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
       echo "#PBS -o $TopOutputLogs/log.realign.ou" >> $qsub_cleanup
       echo "#PBS -e $TopOutputLogs/log.realign.in" >> $qsub_cleanup
       echo "#PBS -q $pbsqueue" >> $qsub_cleanup
       echo "#PBS -m a" >> $qsub_cleanup
       echo "#PBS -M $email" >> $qsub_cleanup
       echo "#PBS -W depend=afterok:$summarydependids" >> $qsub_cleanup
       echo "aprun -n 1 -d $thr $scriptdir/cleanup.sh $outputdir $analysis $TopOutputLogs/log.cleanup.in $TopOutputLogs/log.cleanup.ou $email $TopOutputLogs/qsub.cleanup" >> $qsub_cleanup
       `chmod a+r $qsub_cleanup`
       cleanjobid=`qsub $qsub_cleanup`
       echo $cleanjobid >> $outputdir/logs/CLEANUPpbs
   fi

   `sleep 30s`
   qsub_summary=$TopOutputLogs/qsub.summary.allok
   echo "#PBS -V" > $qsub_summary
   echo "#PBS -A $pbsprj" >> $qsub_summary
   echo "#PBS -N ${pipeid}_summaryok" >> $qsub_summary
   echo "#PBS -l walltime=01:00:00" >> $qsub_summary # 1 hour should be more than enough
   echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
   echo "#PBS -o $TopOutputLogs/log.summary.ou" >> $qsub_summary
   echo "#PBS -e $TopOutputLogs/log.summary.in" >> $qsub_summary
   echo "#PBS -q $pbsqueue" >> $qsub_summary
   echo "#PBS -m a" >> $qsub_summary
   echo "#PBS -M $email" >> $qsub_summary
   if [ `expr ${#cleanjobid}` -gt 0 ]
   then
       echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub_summary
   else
       echo "#PBS -W depend=afterok:$summarydependids" >> $qsub_summary
   fi
   echo "$scriptdir/summary.sh $outputdir $email exitok $reportticket"  >> $qsub_summary
   `chmod a+r $qsub_summary`
   lastjobid=`qsub $qsub_summary`
   echo $lastjobid >> $TopOutputLogs/SUMMARYpbs

   if [ `expr ${#lastjobid}` -lt 1 ]
   then
       echo "at least one job aborted"
       qsub_summary=$TopOutputLogs/qsub.summary.afterany
       echo "#PBS -V" > $qsub_summary
       echo "#PBS -A $pbsprj" >> $qsub_summary
       echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub_summary
       echo "#PBS -l walltime=01:00:00" >> $qsub_summary # 1 hour should be more than enough
       echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
       echo "#PBS -o $TopOutputLogs/log.summary.afterany.ou" >> $qsub_summary
       echo "#PBS -e $TopOutputLogs/log.summary.afterany.in" >> $qsub_summary
       echo "#PBS -q $pbsqueue" >> $qsub_summary
       echo "#PBS -m a" >> $qsub_summary
       echo "#PBS -M $email" >> $qsub_summary
       echo "#PBS -W depend=afterany:$summarydependids" >> $qsub_summary
       echo "$scriptdir/summary.sh $outputdir $email exitnotok $reportticket"  >> $qsub_summary
       `chmod a+r $qsub_summary`
       badjobid=`qsub $qsub_summary`
       echo $badjobid >> $TopOutputLogs/SUMMARYpbs
   fi



echo -e "##########################################################################################
#############################################################################################
############################## release all jobs now #########################################
#############################################################################################"

     realsampleids=$( cat $RealignOutputLogs/REALXSAMPLEXCHRpbs | sed "s/\..*//" | tr "\n" " " )
     realrecalids=$( cat $RealignOutputLogs/REALRECALXLANEXCHRpbs | sed "s/\..*//" | tr "\n" " " )
     mergeids=$( cat $RealignOutputLogs/MERGEXSAMPLEpbs | sed "s/\..*//" | tr "\n" " " )
     mergeVerifyids=$( cat $RealignOutputLogs/MERGEANDVERIFYXSAMPLEpbs | sed "s/\..*//" | tr "\n" " " )


     if [ $skipvcall == "NO" ]
     then
         vcallids=$( cat $VcallOutputLogs/VCALLGATKpbs | sed "s/\..*//" | tr "\n" " " )
     fi

    `qrls -h u $realrecalids`   
    `qrls -h u $mergeids`
    `qrls -h u $realsampleids`
    `qrls -h u $mergeVerifyids`

     if [ $skipvcall == "NO" ]
     then
        `qrls -h u $vcallids`
     fi

