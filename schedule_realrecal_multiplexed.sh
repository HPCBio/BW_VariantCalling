#!/bin/bash
# written in collaboration with Mayo Bioinformatics core group
# and H3A Africa
#
# script to produce improved bams - one per sample 
# from  multiplexed aligned bams 
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 7 ]
then
   MSG="parameter mismatch."
   echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
   exit 1;
fi
   umask 0027
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

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

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
   #targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
   indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
   snpdir=$( cat $runfile | grep -w SNPDIR | cut -d '=' -f2 )   
   realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
   multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
   chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
   indices=$( echo $chrindex | sed 's/:/ /g' )
   sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
   sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
   javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
   skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
   cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

   set +x; echo -e "\n\n#############      cheching type of analysis ###############\n\n" >&2; set -x;
   
   if [ $analysis != "MULTIPLEXED" ]
   then
      MSG="ANALYSIS=$analysis Program=$scriptfile Invalid pipeline program for this type of analysis. This program is for the MULTIPLEXED case only"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi

   set +x; echo -e "\n\n#############      cheching type of sample  ###############\n\n" >&2; set -x;

   if [ $multisample != "YES" -a $multisample != "1" ]
   then
      MSG="MULTISAMPLE=$multisample Invalid value for this type of analysis. This program is for the MULTIPLEXED case only"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi

   set +x; echo -e "\n\n#############      cheching input type      ###############\n\n" >&2; set -x;

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

   set +x; echo -e "\n\n#############      cheching cleanup options   ###############\n\n" >&2; set -x;

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

   set +x; echo -e "\n\n#############skip/include variant calling module ###############\n\n" >&2; set -x;

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

   set +x; echo -e "\n\n#############version of java for gatk ###############\n\n" >&2; set -x;

   if [ -z $javamodule ]
   then
      MSG="Value for JAVAMODULE must be specified in configuration file"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

   set +x; echo -e "\n\n###########directories for output and for tools ###############\n\n" >&2; set -x;


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

   if [ ! -d $refdir/$snpdir ]
   then
      MSG="$snpdir SNP directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi


   set +x; echo -e "\n\n###########    checking scheduler     ###############\n\n" >&2; set -x;

   if [ $run_method != "LAUNCHER" -a $run_method != "QSUB" -a $run_method != "APRUN" -a $run_method != "SERVER" ]
   then
      MSG="Invalid value for RUNMETHOD=$run_method"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   fi

   set +x; echo -e "\n\n########### checking files produced from info sheet ###############\n\n" >&2; set -x;

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

   set +x; echo -e "\n\n" >&2;   
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo "###########################      params ok. Now          creating log folders        ###############" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;  

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


   set +x; echo -e "\n\n" >&2; 
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "###########################   generating regions, intervals, known/knownSites        ############" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x; 


   echo -e "region is the array with snps per chr which will be used by vcallgatk-unifiedGenotyper/haplotypeCaller"
   echo -e "realparms is the array with indels per chr which will be used for  gatk-IndelRealigner"
   echo -e "recalparms is the array with indels per chr which will be used for  gatk-Recalibration"

   echo `date`
   i=1

   for chr in $indices
   do
       cd $refdir/$snpdir
       region[$i]=$( find $PWD -type f -name "${chr}.*.vcf.gz" | sed "s/^/:knownSites:/g" | tr "\n" ":" )

       cd $refdir/$indeldir
       realparms[$i]=$( find $PWD -type f -name "${chr}.*.vcf" | sed "s/^/:known:/g" | tr "\n" ":" )
       recalparms[$i]=$( find $PWD -type f -name "${chr}.*.vcf" | sed "s/^/:knownSites:/g" | tr "\n" ":" )
       (( i++ ))
   done
   echo `date`

   set +x; echo -e "\n\n" >&2;
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo "#####################################  CREATE  DIRECTORIES                    ######################" >&2
   echo "#####################################  FOR RESULTS BY SAMPLE AND BY LANE      ######################" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;

   echo -e "the new directory structure where outputdir: $outputdir"
   echo -e "FOR EACH SAMPLE"
   echo -e "outputdir/sample/realign   outputdir/sample/variant   outputdir/sample/logs"
   echo -e "FOR EACH LANE"
   echo -e "outputdir/sample/lane/align outputdir/sample/lane/realign  outputdir/sample/lane/logs"
   echo -e "soft links are created to alignment results obtained in previous steps"

   while read SampleLine
   do
      if [ `expr ${#SampleLine}` -gt 1 ]
      then
      
          set +x; echo -e "\n\n###### parsing non-empty line, we expect to find five fields###############\n\n" >&2; set -x;
          
	  sample=$( echo "$SampleLine" | cut -f 1 )
	  lane=$( echo "$SampleLine" | cut -f 2 )
	  
          set +x; echo -e "\n\n###### first, let's checking that alignment info exists     ###############\n\n" >&2; set -x;

	  alignedfile=`find $outputdir/$lane/align -name "*.wdups.sorted.bam"`
	  alignedfilehdr=`find $outputdir/$lane/align -name "*.wdups.sorted.bam.header"`
	  alignedfilebai=`find $outputdir/$lane/align -name "*.wdups.sorted.bam.bai"`

	  if [ -s $alignedfile -a -s $alignedfilehdr -a -s $alignedfilebai ]
	  then
              echo -e "alignment files for this lane $lane were found"
          else
              MSG="No aligned bam file(s) found at $outputdir/${lane}/align"
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
              exit 1;              
	  fi
	  
          set +x; echo -e "\n\n###### now we can create the folders for the rest of the analysis  ##############\n\n" >&2; set -x;

	  if [ ! -d $outputdir/${sample} ]
	  then
              echo -e "first time we see this sample $sample"
              echo -e "creating output folders for it"

              mkdir -p $outputdir/${sample}/realign/logs
              mkdir -p $outputdir/${sample}/variant/logs
              mkdir -p $outputdir/${sample}/logs
              mkdir -p $outputdir/${sample}/$lane
              ln -s    $outputdir/$lane/align $outputdir/${sample}/$lane/align
              mkdir -p $outputdir/${sample}/${lane}/realign/logs
              mkdir -p $outputdir/${sample}/${lane}/logs
	  elif [ ! -d $outputdir/${sample}/${lane} ]
	  then
              echo -e "NOT the first time we see this sample $sample"
              echo -e "populating $sample with folders for lane $lane"
              mkdir -p $outputdir/${sample}/$lane
              ln -s $outputdir/$lane/align $outputdir/${sample}/$lane/align
              mkdir -p $outputdir/${sample}/${lane}/realign/logs
              mkdir -p $outputdir/${sample}/${lane}/logs
	  else
             MSG="lane $lane for sample $sample already exists. Error. Lane name has to be unique"
	     echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
             exit 1;
          fi
      fi
   done  <  $outputdir/SAMPLENAMES_multiplexed.list
   # end loop over sampleLanes 



   set +x; echo -e "\n\n" >&2; 
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo "###################################  NESTED LOOP1      starts    here     ##########################" >&2
   echo "########################         call to realign-recalibrate per lane per chr      #################" >&2
   echo "########################   outer loop by Lane; inner loop by chromosome            #################" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;

lanecounter=1 
while read SampleLine
do
    if [ `expr ${#SampleLine}` -gt 1 ]
    then

        set +x; echo -e "\n\n###### processing next non-empty line in SAMPLENAMES_multiplexed.list ##############\n\n" >&2; set -x;


        # now parsing the line just being read
	sample=$( echo "$SampleLine" | cut -f 1 )
	lane=$( echo "$SampleLine" | cut -f 2 )

        set +x; echo -e "\n\n###### realigning recalibrating per lane: $lane sample: $sample      ##############\n\n" >&2; set -x;

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
        
        set +x; echo -e "\n\n###### now putting together the RG line for RGparms                ##############\n\n" >&2; set -x;
        
        if [ `expr ${#RGline}` -lt 1 -o ! -s $RGline ]
        then
           echo -e "RGparms line needs to be recreated from scratch" 
           RGparms=$( grep "^@RG" *wdups*header | sed 's/@RG//' | tr ":" "=" | tr " " ":" | tr "\t" ":" | sed "s/ //g" )
        else 
           RGparms=$( cat $RGline )
        fi

        set +x; echo -e "\n\n###### now entering the loop by chr                                ##############\n\n" >&2; set -x;

       
	chromosomecounter=1
	for chr in $indices
	do

            set +x; echo -e "\n\n###### generating real-recal calls for lane=$lane chr=${chr}   ##############\n\n" >&2; set -x;

            echo "$scriptdir/realrecalLane.sh $lane $RealignOutputDir ${chr}.realrecal.${lane}.output.bam $chr $inputfile $RGparms ${region[$chromosomecounter]} ${realparms[$chromosomecounter]} ${recalparms[$chromosomecounter]} $runfile $RealignLog/log.realrecalLane.$lane.$chr.in $RealignLog/log.realrecalLane.$lane.$chr.ou $email $RealignLog/realrecalLane.${lane}.${chr}" > $RealignLog/realrecalLane.${lane}.${chr}
           

         # end loop over chromosomes
	 echo `date`
	 ((  chromosomecounter++ )) 
	done # done going through chromosomes 
    fi   # skipping empty lines

    # end loop over sampleLanes
    (( lanecounter++ ))
done <  $outputdir/SAMPLENAMES_multiplexed.list


   set +x; echo -e "\n\n" >&2;
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo "###################################   NESTED   LOOP2  starts   here    #############################" >&2
   echo "##################################       second REALIGN  -- BY SAMPLE              #################" >&2
   echo "##################################       followed by vcall -- BY SAMPLE            #################" >&2
   echo "##################################  outer loop by Sample inner loop by chr        ##################" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;
        
        
samplecounter=1 
while read SampleLine
do
    if [ `expr ${#SampleLine}` -gt 1 ]
    then
    
        set +x; echo -e "\n\n###### processing next sample, non-empty line                     ##############\n\n" >&2; set -x;
    
	sample=$( echo "$SampleLine" | cut -f 1 )
        RealignLog=$outputdir/${sample}/realign/logs	
        RealignOutputDir=$outputdir/${sample}/realign
        bamfile=$outputdir/${sample}/realign/${sample}.mergedLanes.bam

        set +x; echo -e "\n\n###### sample=$sample about to enter inner loop by                ##############\n\n" >&2; set -x;


	chromosomecounter=1
	for chr in $indices
	do
            outfile=$outputdir/${sample}/realign/${chr}.${sample}.realignedSample.calmd.bam

            set +x; echo -e "\n\n###### generating realignment calls for sample=$sample chr=${chr} with bam=$bamfile ########\n\n" >&2; set -x;

            echo "$scriptdir/realignSample.sh $RealignOutputDir $outfile $chr $bamfile  ${realparms[$chromosomecounter]} $sample $runfile $RealignLog/log.realignSample.$sample.$chr.in $RealignLog/log.realignSample.$sample.$chr.ou $email $RealignLog/realignSample.${sample}.${chr}" > $RealignLog/realignSample.${sample}.${chr}
           
           if [ $skipvcall == "NO" ]
           then
               set +x; echo -e "\n\n###### generating variant calls for sample=$sample chr=${chr} with bam=$bamfile ########\n\n" >&2; set -x;

               VcallOutputDir=$outputdir/$sample/variant
	       echo "$scriptdir/vcallgatk.sh $VcallOutputDir $RealignOutputDir  $outfile $chr ${region[$chromosomecounter]} $runfile $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.ou $email $VcallOutputDir/logs/vcallgatk.${sample}.${chr}" >> $VcallOutputDir/logs/vcallgatk.${sample}.${chr}
	   fi

         set +x; echo -e "\n\n###### bottom of the loop over chrs                             ##############\n\n" >&2; set -x;
	 echo `date`
	 ((  chromosomecounter++ )) 
	done # done going through chromosomes 
	
    fi   # skipping empty lines

    set +x; echo -e "\n\n###### bottom of the loop over samples                              ##############\n\n" >&2; set -x;

    (( samplecounter++ ))
done <  $outputdir/SAMPLEGROUPS.list

set +x; echo -e "\n\n" >&2; 
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "###################################   QSUB and AnisimovJobLists for  these jobs         ############" >&2
echo "###################################   just considering the launcher option for now #################" >&2
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


   case $run_method in
   "LAUNCHER")
   
   
        set +x; echo -e "\n\n" >&2; 
        echo "####################################################################################################" >&2
	echo -e "LAUNCHER" >&2
        echo "####################################################################################################" >&2

	echo "####################################################################################################" >&2
        echo "####################################################################################################" >&2
        echo " realrecal per lane: next block collects jobs in *AnisimovJobList and generates the qsub script " >&2
        echo "####################################################################################################" >&2
        echo -e "\n\n" >&2; set -x;	

      RealignOutputLogs=$outputdir/logs/realign
      VcallOutputLogs=$outputdir/logs/variant

      echo -e "\n\n ###### loop1 by sample  ###### \n\n"

      while read SampleLine
      do
	  if [ `expr ${#SampleLine}` -gt 1 ]
	  then

              set +x; echo -e "\n\n###### processing next non-empty line in SAMPLENAMES_multiplexed.list ##############\n\n" >&2; set -x;

              # now parsing the line just being read
	      sample=$( echo "$SampleLine" | cut -f 1 )
	      lane=$( echo "$SampleLine" | cut -f 2 )

              RealignOutputDir=$outputdir/$sample/$lane/realign
              truncate -s 0 $RealignOutputLogs/realrecalLane.${lane}.AnisimovJoblist


              set +x; echo -e "\n\n" >&2; 
              echo "####################################################################################################" >&2
              echo -e "collecting all jobs for sample=$sample lane=$lane" >&2
              echo -e "first stop: realign-recalibrate for lane=$lane on all chromosomes" >&2
              echo "####################################################################################################" >&2
              echo -e "\n\n" >&2; set -x;
        
	      for chr in $indices
	      do

                  set +x; echo -e "\n\n###### collecting info for chr $chr ##############\n\n" >&2; set -x;
                  # creating a qsub out of the job file
                  # need to prepend "nohup" and append log file name, so that logs are properly created when Anisimov launches these jobs 

		  realrecal_log=$outputdir/$sample/realign/logs/log.realrecalLane.$lane.$chr.in
                  truncate -s 0 $outputdir/$sample/realign/logs/log.realrecalLane.$lane.$chr.in

		  awk -v awkvar_realrecallog=$realrecal_log '{print "nohup "$0" > "awkvar_realrecallog}' $outputdir/$sample/realign/logs/realrecalLane.${lane}.${chr} > $outputdir/$sample/realign/logs/jobfile.realrecalLane.${lane}.${chr}
		  echo "$outputdir/$sample/realign/logs/ jobfile.realrecalLane.${lane}.${chr}" >> $RealignOutputLogs/realrecalLane.${lane}.AnisimovJoblist

              done

              set +x; echo -e "\n\n" >&2; 
              echo "####################################################################################################" >&2              
              echo "########### end loop over chr" >&2  
              echo "####################################################################################################" >&2 
              echo -e "######## QSUB for this Anisimov Launcher joblists"  >&2 
              echo "####################################################################################################" >&2 
              echo -e "\n\n" >&2; set -x;
              
              
              qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecalLane.${lane}.AnisimovLauncher

              set +x; echo -e "\n\n###### # appending the generic header to the qsub ##############\n\n" >&2; set -x;
                  

              cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov

              echo "#PBS -N ${pipeid}_realrecalxlane_${lane}" >> $qsub_realrecal_anisimov
              echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal_anisimov
              echo "#PBS -o $RealignOutputLogs/log.realrecalLane.${lane}.ou" >> $qsub_realrecal_anisimov
              echo "#PBS -e $RealignOutputLogs/log.realrecalLane.${lane}.in" >> $qsub_realrecal_anisimov

              echo "#PBS -l nodes=$chromosomecounter:ppn=$thr" >> $qsub_realrecal_anisimov
              # the actual command
              echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realrecalLane.${lane}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realrecalLane.${lane}.AnisimovLauncher.log" >> $qsub_realrecal_anisimov

          fi

          set +x; echo -e "\n\n###### bottom of the loop over samples                              ##############\n\n" >&2; set -x;


      done <  $outputdir/SAMPLENAMES_multiplexed.list

      set +x; echo -e "\n\n" >&2; 
      echo "####################################################################################################" >&2
      echo "####################################################################################################" >&2
      echo "  next block collects jobs in *AnisimovJobList and generates the qsub script " >&2
      echo "  for mergeLanesPerSample realignSample verifySample vcallgatk " >&2
      echo "####################################################################################################" >&2
      echo "####################################################################################################" >&2
      echo -e "\n\n" >&2; set -x;

      RealignOutputLogs=$outputdir/logs/realign
      VcallOutputLogs=$outputdir/logs/variant

      while read SampleLine
      do
	  if [ `expr ${#SampleLine}` -gt 1 ]
	  then

              set +x; echo -e "\n\n###### processing next non-empty line                               ##############\n\n" >&2; set -x;


	      sample=$( echo "$SampleLine" | cut -f 1 )
	      lanes=$( echo "$SampleLine" | cut -f 2 )
              lanes=$( echo $lanes | tr " " ":" )

              RealignLog=$outputdir/${sample}/realign/logs
              RealignOutputDir=$outputdir/${sample}/realign
              sampleOutPrefix=$outputdir/${sample}/realign/${sample}
              sLB="multiplexed"

              truncate -s 0 $RealignOutputLogs/mergeLanes.${sample}.AnisimovJoblist
              truncate -s 0 $RealignOutputLogs/realignSample.${sample}.AnisimovJoblist
              truncate -s 0 $RealignOutputLogs/verifySample.${sample}.AnisimovJoblist
     
	      if [ $skipvcall == "NO" ]
	      then
		  VcallOutputDir=$outputdir/$sample/variant
		  truncate -s 0 $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist
              fi

              set +x; echo -e "\n\n######  constructing the qsub for mergeLanes                        ##############\n\n" >&2; set -x;

              qsub_mergeLanes_anisimov=$RealignOutputLogs/qsub.mergeLanes.${sample}.AnisimovLauncher
              cat $outputdir/qsubGenericHeader > $qsub_mergeLanes_anisimov

	      echo "#PBS -N ${pipeid}_mergeLanes_${sample}" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -l walltime=$pbscpu" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -o $RealignLog/log.mergeLanes.${sample}.ou" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -e $RealignLog/log.mergeLanes.${sample}.in" >> $qsub_mergeLanes_anisimov
	      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_mergeLanes_anisimov
              echo "aprun -n 1 -N 1 -d $thr $scriptdir/mergeLanesPerSample.sh $runfile $sample $sLB $lanes ${sampleOutPrefix}.mergedLanes.bam $outputdir $RealignLog/log.mergeLanes.$sample.in $RealignLog/log.mergeLanes.$sample.ou $email $qsub_mergeLanes_anisimov" >> $qsub_mergeLanes_anisimov

              set +x; echo -e "\n\n######  constructing the qsub for verifysample                       ##############\n\n" >&2; set -x;


              qsub_verifySample_anisimov=$RealignOutputLogs/qsub.verifySample.${sample}.AnisimovLauncher
              cat $outputdir/qsubGenericHeader > $qsub_verifySample_anisimov

	      echo "#PBS -N ${pipeid}_verifySample_${sample}" >> $qsub_verifySample_anisimov
	      echo "#PBS -l walltime=$pbscpu" >> $qsub_verifySample_anisimov
	      echo "#PBS -o $RealignLog/log.verifySample.${sample}.ou" >> $qsub_verifySample_anisimov
	      echo "#PBS -e $RealignLog/log.verifySample.${sample}.in" >> $qsub_verifySample_anisimov
	      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_verifySample_anisimov
              echo "aprun -n 1 -N 1 -d $thr $scriptdir/verifySample.sh $runfile $sample $outputdir/${sample}/realign  ${sample}.verified $outputdir $RealignLog/log.verifySample.$sample.in $RealignLog/log.verifySample.$sample.ou $email $qsub_verifySample_anisimov" >> $qsub_verifySample_anisimov


              set +x; echo -e "\n\n######  collecting jobs in AnsimovJobLists for  realignSample and vcallgatk ##############\n\n" >&2; set -x;

	      for chr in $indices
	      do
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

              set +x; echo -e "\n\n######  constructing the qsub job  for realignSample                   ##############\n\n" >&2; set -x;


	      qsub_realignSample_anisimov=$RealignOutputLogs/qsub.realignSample.${sample}.AnisimovLauncher
              cat $outputdir/qsubGenericHeader > $qsub_realignSample_anisimov

	      echo "#PBS -N ${pipeid}_realignSample_${sample}" >> $qsub_realignSample_anisimov
	      echo "#PBS -l walltime=$pbscpu" >> $qsub_realignSample_anisimov
	      echo "#PBS -o $RealignOutputLogs/log.realignSample.${sample}.ou" >> $qsub_realignSample_anisimov
	      echo "#PBS -e $RealignOutputLogs/log.realignSample.${sample}.in" >> $qsub_realignSample_anisimov
	      echo "#PBS -l nodes=$chromosomecounter:ppn=$thr" >> $qsub_realignSample_anisimov
	      echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realignSample.${sample}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realignSample.${sample}.AnisimovLauncher.log" >> $qsub_realignSample_anisimov


	      if [ $skipvcall == "NO" ]
	      then	      
              set +x; echo -e "\n\n######  constructing the qsub job  for varcall with   ${sample}         ##############\n\n" >&2; set -x;

                  qsub_vcallgatk_anisimov=$VcallOutputLogs/qsub.vcallgatk.${sample}.AnisimovLauncher
		  cat $outputdir/qsubGenericHeader > $qsub_vcallgatk_anisimov

		  echo "#PBS -N ${pipeid}_vcallgatk_${sample}" >> $qsub_vcallgatk_anisimov
		  echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk_anisimov
		  echo "#PBS -o $VcallOutputLogs/log.vcallgatk.${sample}.ou" >> $qsub_vcallgatk_anisimov
		  echo "#PBS -e $VcallOutputLogs/log.vcallgatk.${sample}.in" >> $qsub_vcallgatk_anisimov
		  echo "#PBS -l nodes=$chromosomecounter:ppn=$thr" >> $qsub_vcallgatk_anisimov
                  echo "ulimit -s unlimited " >> $qsub_vcallgatk_anisimov
                  echo "export APRUN_XFER_LIMITS=1 " >> $qsub_vcallgatk_anisimov
		  echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist /bin/bash > $VcallOutputLogs/vcallgatk.${sample}.AnisimovLauncher.log" >> $qsub_vcallgatk_anisimov
	      fi
	  fi
      done <  $outputdir/SAMPLEGROUPS.list

      set +x; echo -e "\n\n" >&2; 
      echo "####################################################################################################" >&2
      echo "####################################################################################################" >&2
      echo "####              arranging execution order with job dependencies            #######################" >&2
      echo "####################################################################################################" >&2      
      echo "####################################################################################################" >&2
      echo -e "\n\n" >&2; set -x;
        
      RealignOutputLogs=$outputdir/logs/realign
      VcallOutputLogs=$outputdir/logs/variant
      
      truncate -s 0 $RealignOutputLogs/REALXSAMPLEXCHRpbs
      truncate -s 0 $RealignOutputLogs/REALRECALXLANEXCHRpbs
      truncate -s 0 $RealignOutputLogs/MERGEXSAMPLEpbs
      truncate -s 0 $RealignOutputLogs/VERIFYXSAMPLEpbs

      if [ $skipvcall == "NO" ]
      then
          truncate -s 0 $VcallOutputLogs/VCALLGATKpbs
      fi
      
      cd $RealignOutputLogs 
      # first we schedule the realrecalxlane
      # FIX: gather all jobids for a single sample and then create the PBS depend line for sample

      set +x; echo -e "\n\n######  next loop to schedule the realrecalxlane         ##############\n\n" >&2; set -x;
      
      while read SampleLine
      do
	  if [ `expr ${#SampleLine}` -gt 1 ]
	  then
	  
              set +x; echo -e "\n\n###### processing next non-empty line in SAMPLEGROUPS.list        ##############\n\n" >&2; set -x;	  

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
              sed -i "2i #PBS -W depend=afterok:$RealignSample_job" $RealignOutputLogs/qsub.verifySample.${sample}.AnisimovLauncher    

              VerifySample_job=`qsub $RealignOutputLogs/qsub.verifySample.${sample}.AnisimovLauncher`
              `qhold -h u $VerifySample_job`

              sed -i "2i #PBS -W depend=afterok:$VerifySample_job" $VcallOutputLogs/qsub.vcallgatk.${sample}.AnisimovLauncher             
              VcallSample_job=`qsub $VcallOutputLogs/qsub.vcallgatk.${sample}.AnisimovLauncher`  
              `qhold -h u $VcallSample_job`

              echo $mergeSample_job   >> $RealignOutputLogs/MERGEXSAMPLEpbs
              echo $RealignSample_job >> $RealignOutputLogs/REALXSAMPLEXCHRpbs
              echo $VerifySample_job  >> $RealignOutputLogs/VERIFYXSAMPLEpbs
              echo $VcallSample_job   >> $VcallOutputLogs/VCALLGATKpbs               
          fi
      done <  $outputdir/SAMPLEGROUPS.list
   ;;
   esac

   set +x; echo -e "\n\n" >&2; 
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo "###########################    wrap up and produce summary table        ############################" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;
        

   if [ $skipvcall == "NO" ]
   then
      summarydependids=$( cat $VcallOutputLogs/VCALLGATKpbs | sed "s/\..*//" | tr "\n" ":" )
   else
      summarydependids=$( cat $RealignOutputLogs/VERIFYXSAMPLEpbs | sed "s/\..*//" | tr "\n" ":" )
   fi

   lastjobid=""
   cleanjobid=""

   if [ $cleanupflag == "YES" ]
   then
       qsub_cleanup=$TopOutputLogs/qsub.cleanup
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
       #`chmod a+r $qsub_cleanup`
       cleanjobid=`qsub $qsub_cleanup`
       echo $cleanjobid >> $outputdir/logs/CLEANUPpbs
   fi

   `sleep 30s`
   qsub_summary=$TopOutputLogs/qsub.summary.afterany
   echo "#PBS -A $pbsprj" >> $qsub_summary
   echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub_summary
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
       echo "#PBS -W depend=afterany:$summarydependids" >> $qsub_summary
   fi
   echo "$scriptdir/summary.sh $runfile $email exitafterany $reportticket"  >> $qsub_summary
   #`chmod a+r $qsub_summary`
   lastjobid=`qsub $qsub_summary`
   echo $lastjobid >> $TopOutputLogs/SUMMARYpbs


   echo "rewrite summary report if all jobs finished ok"
   qsub_summaryany=$TopOutputLogs/qsub.summary.afterok
   echo "#PBS -A $pbsprj" >> $qsub_summaryany
   echo "#PBS -N ${pipeid}_summary_afterok" >> $qsub_summaryany
   echo "#PBS -l walltime=01:00:00" >> $qsub_summaryany # 1 hour should be more than enough
   echo "#PBS -l nodes=1:ppn=1" >> $qsub_summaryany
   echo "#PBS -o $TopOutputLogs/log.summary.afterok.ou" >> $qsub_summaryany
   echo "#PBS -e $TopOutputLogs/log.summary.afterok.in" >> $qsub_summaryany
   echo "#PBS -q $pbsqueue" >> $qsub_summaryany
   echo "#PBS -m a" >> $qsub_summaryany
   echo "#PBS -M $email" >> $qsub_summaryany
   if [ `expr ${#cleanjobid}` -gt 0 ]
   then
       echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub_summaryany
   else
       echo "#PBS -W depend=afterok:$summarydependids" >> $qsub_summaryany
   fi
   echo "$scriptdir/summary.sh $runfile $email exitok $reportticket"  >> $qsub_summaryany
   #`chmod a+r $qsub_summaryany`
   badjobid=`qsub $qsub_summaryany`
   echo $badjobid >> $TopOutputLogs/SUMMARYpbs



   set +x; echo -e "\n\n" >&2; 
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo "###########################    release  all   job                       ############################" >&2
   echo "####################################################################################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;
        

     realsampleids=$( cat $RealignOutputLogs/REALXSAMPLEXCHRpbs | sed "s/\..*//" | tr "\n" " " )
     verifysampleids=$( cat $RealignOutputLogs/VERIFYXSAMPLEpbs | sed "s/\..*//" | tr "\n" " " )
     realrecalids=$( cat $RealignOutputLogs/REALRECALXLANEXCHRpbs | sed "s/\..*//" | tr "\n" " " )
     mergeids=$( cat $RealignOutputLogs/MERGEXSAMPLEpbs | sed "s/\..*//" | tr "\n" " " )
     if [ $skipvcall == "NO" ]
     then
         vcallids=$( cat $VcallOutputLogs/VCALLGATKpbs | sed "s/\..*//" | tr "\n" " " )
     fi

     `qrls -h u $realsampleids`    
     `qrls -h u $verifysampleids`
     `qrls -h u $realrecalids`   
     `qrls -h u $mergeids`
     if [ $skipvcall == "NO" ]
     then
         `qrls -h u $vcallids`
     fi

