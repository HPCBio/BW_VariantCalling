#!/bin/bash
# written in collaboration with Mayo Bioinformatics core group
#
#  script to realign and recalibrate the aligned file(s)
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 7 ]
then
   MSG="parameter mismatch."
   echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
   exit 1;
else
   set -x
   echo `date`
######################## realigndir and listfiles removed
######################## RealignOutputLogs replaced with outputdir: top level output folder
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
                ####################################################################################################"
                #####################################                       ######################################## 
echo -e "\n\n\n ##################################### PARSING RUN INFO FILE ######################################## \n\n\n"
                #####################################                       ######################################## 
                #################################################################################################### 





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
   realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
   multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
   samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 )
   chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
   #indices=$( echo $chrindex | sed 's/^/chr/' | sed 's/:/ chr/g' )
   indices=$( echo $chrindex | sed 's/:/ /g' )
   sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
   sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
   sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
   javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
   skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
   cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

   echo "#################################################################################"
   echo -e cheching type of analysis
   echo "#################################################################################"

   if [ $analysis == "MULTIPLEXED" ]
   then
      MSG="ANALYSIS=$analysis Program=$scriptfile Invalid pipeline program for this type of analysis. This program is for the NON-MULTIPLEXED case only"
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

  # if [ -s $refdir/$dbSNP ]
  # then
  #    realparms="-known:$refdir/$dbSNP"
  #    recalparms="--knownSites:$refdir/$dbSNP"
  # else
  #    MSG="dbSNP not found"
  #    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
  #    exit 1;
  # fi
   
   if [ ! -d $refdir/$indeldir ]
   then
      MSG="$indeldir indel directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi


   echo "####################################################################################################"
   echo "####################################################################################################"
   echo "###########################      parameters ok. Now      creating log folders        ###############"
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

   echo -e "\n####################################################################################################"
   echo -e "\n####################################################################################################"
   echo -e "\n###########################   generating regions, intervals, known/knownSites        ###############"
   echo -e "\n####################################################################################################"
   echo -e "\n####################################################################################################"
   echo -e "\n ###### region is the array with snps per chr which will be used by vcallgatk-unifiedGenotyper/haplotypeCaller"
   echo -e "\n ###### realparms is the array with indels per chr which will be used for  gatk-IndelRealigner"
   echo -e "\n ###### recalparms is the array with indels per chr which will be used for  gatk-Recalibration"
   echo -e "\n####################################################################################################\n\n"

   echo `date`
   i=1

   for chr in $indices
   do
       cd $refdir/vcf_per_chr
       snps=`find $PWD -type f -name "${chr}.*.vcf.gz"`
       region[$i]=$( echo $snps | sed "s/\/projects/:knownSites:\/projects/g" | sed "s/ //g" | tr "\n" ":" )

       cd $refdir/$indeldir
       indels=`find $PWD -type f -name "${chr}.*.vcf"`
       realparms[$i]=$( echo $indels | sed "s/\/projects/:known:\/projects/g" | sed "s/ //g" |tr "\n" ":" )
       recalparms[$i]=$( echo $indels | sed "s/\/projects/:knownSites:\/projects/g" | sed "s/ //g" | tr "\n" ":" )
       (( i++ ))
   done
   echo `date`

   echo -e "\n####################################################################################################"
   echo -e "\n####################### done  generating regions, intervals, known/knownSites        ###############"
   echo -e "\n####################################################################################################"


                ####################################################################################################
                #####################################                       ########################################
echo -e "\n\n\n #####################################  CREATE OUTPUT  DIRECTORIES  ################################# \n\n\n"
                #####################################                       ########################################
                ####################################################################################################

   if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
   then 
        TheInputFile=$outputdir/SAMPLENAMES_multiplexed.list
   else
        TheInputFile=$outputdir/SAMPLENAMES.list
   fi


   while read SampleLine
   do
      if [ `expr ${#SampleLine}` -gt 1 ]
      then
          ## parsing non-empty line
          
          if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
          then
              echo -e "this line has five fields, we just need the first field to create folders"
	      sample=$( echo "$SampleLine" | cut -f 1 )
	  else
              echo -e "this line has one field, we need to parse the sample name out of the filename"
	      sample=$( echo "$SampleLine" )
	  fi
	  
          echo -e "######################################################################"
	  echo -e "########## first, let's checking that alignment info exists  #########"

	  alignedfile=`find $outputdir/$sample/align -name "*.wdups.sorted.bam"`
	  alignedfilehdr=`find $outputdir/$sample/align -name "*.wdups.sorted.bam.header"`
	  alignedfilebai=`find $outputdir/$sample/align -name "*.wdups.sorted.bam.bai"`

	  if [ -s $alignedfile -a -s $alignedfilehdr -a -s $alignedfilebai ]
	  then
              echo -e "alignment files for this sample $sample were found at $outputdir/${sample}/align"
          else
              MSG="No aligned bam file(s) found at $outputdir/${sample}/align"
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
              exit 1;              
	  fi

          echo -e "################################################################################"
          echo -e "########## now we can create the folders for the rest of the analysis  #########"
	  if [ -d $outputdir/${sample} ]
	  then
              echo -e "creating output folders for sample=$sample"

              mkdir -p $outputdir/${sample}/realign/logs
              mkdir -p $outputdir/${sample}/variant/logs
          fi

          # resetting logs
	  rm $outputdir/${sample}/realign/logs/*
          rm $outputdir/${sample}/variant/logs/*

      fi  # end processing non-empty lines
   done  <  $TheInputFile
   # end loop over samples




                #############################################################################################################
                ###################################                                ##########################################
echo -e "\n\n\n ###################################       main loops starts here   ###################################### 
                ########################   outer loops by sample; inner loops by chromosome/region   #################### \n\n\n"
                ########################                                                      ###############################
                #############################################################################################################



      echo -e "\n####################################################################################################"
      echo -e "\n##########                      loop by sample starts here                              ###########"
      echo -e "\n####################################################################################################"

      samplecounter=1 
      while read SampleLine
      do
          if [ `expr ${#SampleLine}` -gt 1 ]
          then
              ## processing non-empty line
	      echo `date`
	      echo -e " ######     processing this line: $SampleLine #####\n"

              if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
              then
                   echo -e "\n ###### this line has five fields, we just need the first field to create folders"
	           sample=$( echo "$SampleLine" | cut -f 1 )
	      else
                   echo -e "\n ###### this line has one field, we will use filename as samplename"
	           sample=$( echo "$SampleLine" )
	      fi

              echo -e "\n####################################################################################################"	       
              echo -e "\n######################################  get aligned bam files  #####################################"
              echo -e "\n####################################################################################################"
              echo -e "\n## we already checked that these file exist when we were creating folders in the previous block  ###"
              echo -e "\n####################################################################################################"
              
              cd $outputdir/${sample}/align
	      aligned_bam=`find ./ -name "*.wdups.sorted.bam"`
	      
              # now check that there is only one bam file
              aligned_bam=$outputdir/${sample}/align/${aligned_bam}
              aligned_bam=( $aligned_bam ) # recast variable as an array and count the number of members
              if [ ${#aligned_bam[@]} -ne 1 ]
              then
                 MSG="more than one bam file found in $outputdir/${sample}/align"
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
                 exit 1;
              fi

              aligned_bam="${aligned_bam[*]}" # recast variable back as string	      
              header=${aligned_bam}.header
 
              echo -e "\n####################################################################################################"
              echo -e "\n##########                   now putting together the RG line for RGparms                ###########"
	      echo -e "\n####################################################################################################"
 
              #truncate -s 0 ${aligned_bam}.RGline
	      #`grep "RG.*ID" $header | cut -f 2` >> ${aligned_bam}.RGline
              RGline=${aligned_bam}.RGline
              #RGline=$( sed -i 's/ID/RGID/' "$RGline" )
              #RGline=$( sed -i 's/ / RG/g' "$RGline" )
              #RGparms=$( cat "$RGline" )
              
              if [ `expr ${#RGline}` -lt 1 -o ! -s $RGline ]
              then
                 echo -e "RGparms line needs to be recreated from scratch" 
                 RGparms=$( grep "^@RG" *wdups*header | sed 's/@RG//' | tr ":" "=" | tr " " ":" | tr "\t" ":" | sed "s/ //g" )
              else 
                 RGparms=$( cat $RGline )
              fi

 
              echo -e "\n####################################################################################################"
              echo -e "\n##########                   now defining paths to output folders                        ###########"
	      echo -e "\n####################################################################################################"
              
              RealignOutputDir=$outputdir/${sample}/realign
              VcallOutputDir=$outputdir/${sample}/variant
              RealignLog=$outputdir/${sample}/realign/logs
              truncate -s 0 $RealignOutputLogs/verifySample.${sample}.AnisimovJoblist
              truncate -s 0 $RealignOutputLogs/realrecal.${sample}.AnisimovJoblist

              if [ $skipvcall == "NO" ]
              then
                  truncate -s 0 $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist
              fi
              
              echo -e "\n####################################################################################################"
              echo -e "\n##########                   loop by chromosome starts here                              ###########"
	      echo -e "\n####################################################################################################"

	      chromosomecounter=1
	      for chr in $indices
	      do
		  echo -e "\n #########################################################################################################\n"
		  echo -e "\n ######    generating real-recal, vcallgatk  calls for chr=${chr}                                   ######\n"
		  echo -e "\n #########################################################################################################\n"      
		  echo `date`
                  realrecaloutputfile=${chr}.realrecal.${sample}.calmd.bam
                  echo -e "\n#######################   assemble the real-recall/var call sub-block                   ################################\n"
                  echo "$scriptdir/realrecal.sh $RealignOutputDir $realrecaloutputfile $chr $aligned_bam $RGparms ${region[$chromosomecounter]} ${realparms[$chromosomecounter]} ${recalparms[$chromosomecounter]} $runfile $flag $RealignOutputDir/logs/log.realrecal.$sample.$chr.in $RealignOutputDir/logs/log.realrecal.$sample.$chr.ou $email $RealignOutputDir/logs/realrecal.${sample}.${chr}" > $RealignOutputDir/logs/realrecal.${sample}.${chr}
           
                  if [ $skipvcall == "NO" ]
                  then
		      echo -e "\n#######################   assemble the vcallgatk call sub-block                     ################################\n"
		      echo "$scriptdir/vcallgatk.sh $VcallOutputDir  $RealignOutputDir $realrecaloutputfile $chr ${region[$chromosomecounter]} $runfile $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.ou $email $VcallOutputDir/logs/vcallgatk.${sample}.${chr}" >> $VcallOutputDir/logs/vcallgatk.${sample}.${chr}
                  fi

		  ((  chromosomecounter++ )) 
	      done # done going through chromosomes 
	      
	      
              (( samplecounter++ ))
          fi # done processing non-empty lines
      done <  $TheInputFile
      # end loop over samples

   # at the end of this set of nested loops, the variables chromosomecounter and samplecounter
   # reflect, respectively, number_of_chromosomes+1 and number_of_samples+1,
   # which is exactly the number of nodes required for anisimov launcher 


      echo -e "\n####################################################################################################"
      echo -e "\n##########                      main loop ends here                                      ###########"
      echo -e "\n####################################################################################################"


                #############################################################################################################
                ###################################                             #############################################
      echo -e "\n\n\n ###################################   now schedule these jobs   ###################################### \n\n\n"
                ###################################                             #############################################
                #############################################################################################################

   #set +x; echo -e "\n ### update autodocumentation script ### \n"; set -x;
   #echo -e "# @begin RealignRecalibrate_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
   #echo -e "   # @in sample_chr @as aligned_bam_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
   #RealignedBAMTemplate="{SampleName}/realign/{chromosome}.realrecal.${SampleName}.output.bam"
   #echo -e "   # @out realrecal  @as  realigned_bam_per_chromosome @URI ${RealignedBAMTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
   #echo -e "# @end RealignRecalibrate_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh

   #echo -e "# @begin VariantCalling_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
   #echo -e "   # @in  realrecal  @as  realigned_bam_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
   #VariantTemplate=${RealignedBAMTemplate}.raw.all.vcf
   #echo -e "   # @out vcf @as output_variants @URI ${VariantTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
   #echo -e "# @end VariantCalling_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh

   #echo -e "# @out vcf @as output_variants @URI sample_name/variant/chr_name.vcf" >> $outputdir/WorkflowAutodocumentationScript.sh
   #WorkflowName=`basename $outputdir`
   #echo -e "# @end $WorkflowName" >> $outputdir/WorkflowAutodocumentationScript.sh


   case $run_method in
   "APRUN")
      echo -e "\n\n ###### run_method=$run_method. generating and scheduling qsubs... ###### \n\n"

      truncate -s 0 $RealignOutputLogs/VERIFYXSAMPLEpbs
      truncate -s 0 $RealignOutputLogs/REALRECALpbs
      truncate -s 0 $VcallOutputLogs/VCALGATKpbs

      echo -e "\n\n ###### loop1 by sample  ###### \n\n"
      while read SampleLine
      do

          if [ `expr ${#SampleLine}` -gt 1 ]
          then
              ## processing non-empty line
              if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
              then
                   echo -e "this line has five fields, we just need the first field to create folders"
	           sample=$( echo "$SampleLine" | cut -f 1 )
	      else
                   echo -e "this line has one field, we will use filename as samplename"
	           sample=$( echo "$SampleLine" )
	      fi
              RealignOutputDir=$outputdir/$sample/realign
              VcallOutputDir=$outputdir/${sample}/variant 

	      echo -e "\n\n ###### loop2 by chromosome  ###### \n\n"              
              for chr in $indices
              do

                 qsub_realrecal=$RealignOutputDir/logs/qsub.realrecal.${sample}.${chr}

                 ###############################
                 echo -e "\n################# constructing qsub for realrecal\n"
                 # appending the generic header to the qsub
                 cat $outputdir/qsubGenericHeader > $qsub_realrecal
                 echo "#PBS -N ${pipeid}_realrecal.${sample}.${chr}" >> $qsub_realrecal
                 echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal
                 echo "#PBS -o $RealignOutputDir/logs/log.realrecal.${sample}.${chr}.ou" >> $qsub_realrecal
                 echo "#PBS -e $RealignOutputDir/logs/log.realrecal.${sample}.${chr}.in" >> $qsub_realrecal
                 echo "#PBS -W depend=afterok:$split_job" >> $qsub_realrecal
                 echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_realrecal

                 echo "aprun -n 1 -N 1 -d $thr /bin/bash $RealignOutputDir/logs/realrecal.${sample}.${chr}" >> $qsub_realrecal

                 realrecal_job=`qsub $qsub_realrecal`
                 `qhold -h u $realrecal_job`
                 #`qrls -h u $split_job` 
                 echo $realrecal_job >> $RealignOutputLogs/REALRECALpbs



                 if [ $skipvcall == "NO" ]
                 then
                     qsub_vcall=$VcallOutputDir/logs/qsub.vcallgatk.${sample}.${chr}



                     ###############################
                     echo -e "\n################# constructing qsub for vcall\n"
                     # appending the generic header to the qsub
                     cat $outputdir/qsubGenericHeader > $qsub_vcall
                     echo "#PBS -N ${pipeid}_vcall_${SampleName}.${chr}" >> $qsub_vcall
                     echo "#PBS -l walltime=$pbscpu" >> $qsub_vcall
                     echo "#PBS -o $VcallOutputDir/logs/log.vcall.${sample}.${chr}.ou" >> $qsub_vcall
                     echo "#PBS -e $VcallOutputDir/logs/log.vcall.${sample}.${chr}.in" >> $qsub_vcall
                     echo "#PBS -W depend=afterok:$realrecal_job" >> $qsub_vcall
                     echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_vcall
   
                     echo "aprun -n 1 -d $thr /bin/bash $VcallOutputDir/logs/vcallgatk.${SampleName}.${chr}" >> $qsub_vcall
   
                     vcall_job=`qsub $qsub_vcall`
                     #`qrls -h u $realrecal_job`
                     echo $vcall_job >> $VcallOutputLogs/VCALGATKpbs

                 fi
              done # going through chromosomes for this sample
         fi  # non-empty line of file
      done <  $TheInputFile
      # end loop over samples
   ;;
   "QSUB")
      # will add later
   ;;
   "LAUNCHER")

      echo -e "\n\n ###### run_method=$run_method. populate list of jobs first  and launching them later via launcher... ###### \n\n"

      echo -e "\n\n ###### loop1 by sample to create  ONE list of jobs sample with all chromosomes inside it  ###### \n\n"

      while read SampleLine
      do
           if [ `expr ${#SampleLine}` -gt 1 ]
           then
               ## processing non-empty line
               if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
               then
                     echo -e "this line has five fields, we just need the first field to create folders"
	             sample=$( echo "$SampleLine" | cut -f 1 )
	       else
                     echo -e "this line has one field, we will use filename as samplename"
	             sample=$( echo "$SampleLine" )
	       fi
	         
	       ## defining paths and reset joblists
	         
               RealignOutputDir=$outputdir/$sample/realign
               VcallOutputDir=$outputdir/${sample}/variant  
                          
	       echo -e "\n\n ###### loop2 by chromosome to populate the joblist by writing ONE line for each job call to a sample-chr pair ###### \n\n"         
                 
               for chr in $indices
               do
         
                     # creating a qsub out of the job file
                     # need to prepend "nohup" and append log file name, so that logs are properly created when Anisimov launches these jobs 

                     realrecal_log=$RealignOutputDir/logs/log.realrecal.${sample}.$chr.in

		     truncate -s 0 $RealignOutputDir/logs/jobfile.realrecal.${sample}.${chr}
                     truncate -s 0 $RealignOutputDir/logs/log.realrecal.${sample}.$chr.in

                     awk -v awkvar_realrecallog=$realrecal_log '{print "nohup "$0" > "awkvar_realrecallog}' $RealignOutputDir/logs/realrecal.${sample}.${chr} > $RealignOutputDir/logs/jobfile.realrecal.${sample}.${chr}
                     echo "$RealignOutputDir/logs/ jobfile.realrecal.${sample}.${chr}" >> $RealignOutputLogs/realrecal.${sample}.AnisimovJoblist

                     if [ $skipvcall == "NO" ]
                     then
                         vcall_log=$VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in

                         truncate -s 0 $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in
                         truncate -s 0 $VcallOutputDir/logs/jobfile.vcallgatk.${sample}.${chr}

                         awk -v awkvar_vcalllog=$vcall_log '{print "nohup "$0" > "awkvar_vcalllog}' $VcallOutputDir/logs/vcallgatk.${sample}.${chr} > $VcallOutputDir/logs/jobfile.vcallgatk.${sample}.${chr}
                         echo "$VcallOutputDir/logs/ jobfile.vcallgatk.${sample}.${chr}" >> $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist
                     fi
               done # end loop2 by chromosome


	       echo -e "\n\n ######                         outside loop2, joblist should be populated now                              ###### \n"         
               echo -e "\n\n ###### putting together the other pieces of the qsub file and then scheduling Anisimov Launcher joblists   ###### \n\n\n"
               # appending the generic header to the qsub
               
               qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecal.${sample}.AnisimovLauncher
               cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov

               if [ $skipvcall == "NO" ]
               then
                     qsub_vcallgatk_anisimov=$VcallOutputLogs/qsub.vcalgatk.${sample}.AnisimovLauncher
                     cat $outputdir/qsubGenericHeader > $qsub_vcallgatk_anisimov
               fi

               ###############
               ############### constructing the qsub for realrecal
               echo "#PBS -N ${pipeid}_realrecal_${sample}" >> $qsub_realrecal_anisimov
               echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal_anisimov
               echo "#PBS -o $RealignOutputLogs/log.realrecal.${sample}.ou" >> $qsub_realrecal_anisimov
               echo -e "#PBS -e $RealignOutputLogs/log.realrecal.${sample}.in\n" >> $qsub_realrecal_anisimov
               # the dependency on split_bam_by_chromosome job will be added when it is scheduled in the loop below

               # realrecal and vcall actually use multithreaded processes,
               # so we will give each chromosome its own node  +1 for the launcher
               # chromosmecounter is the variable that has this value already setup
               echo -e "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_realrecal_anisimov
               # the actual command
               echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realrecal.${sample}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realrecal.${sample}.AnisimovLauncher.log" >> $qsub_realrecal_anisimov

               if [ $skipvcall == "NO" ]
               then
                     ###############
                     ############### constructing the qsub for vcallgatk
                     echo "#PBS -N ${pipeid}_vcallgatk_${sample}" >> $qsub_vcallgatk_anisimov
                     echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk_anisimov
                     echo "#PBS -o $VcallOutputLogs/log.vcallgatk.${sample}.ou" >> $qsub_vcallgatk_anisimov
                     echo -e "#PBS -e $VcallOutputLogs/log.vcallgatk.${sample}.in\n" >> $qsub_vcallgatk_anisimov
                     # the dependency on realrecal job will be added when it is scheduled in the loop below
                     echo -e "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_vcallgatk_anisimov
                     # the actual command
                    echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $VcallOutputLogs/vcallgatk.${sample}.AnisimovJoblist /bin/bash > $VcallOutputLogs/vcallgatk.${sample}.AnisimovLauncher.log" >> $qsub_vcallgatk_anisimov
               fi

           fi # end non-empty line
      done <  $TheInputFile
      # end loop over samples


      ###############
      ############### 
      echo -e "\n ######                        Outside loop1.                                                  ######"
      echo -e "\n ###### Now. Going through samples again to schedule them in the right order                   ######"
      echo -e "\n ###### all realrecal jobs followed by verifysample and finally  vcallgatk                     ######"


      # reset the lists
      truncate -s 0 $RealignOutputLogs/VERIFYXSAMPLEpbs
      truncate -s 0 $RealignOutputLogs/REALRECALpbs
      if [ $skipvcall == "NO" ]
      then
         truncate -s 0 $VcallOutputLogs/VCALGATKpbs
      fi
      cd $RealignOutputLogs # so that whatever temp folders and pbs notifications would go there


      while read SampleLine
      do
          if [ `expr ${#SampleLine}` -gt 1 ]
          then
               ## processing non-empty line
               if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
               then
                     echo -e "this line has five fields, we just need the first field to create folders"
	             sample=$( echo "$SampleLine" | cut -f 1 )
	       else
                     echo -e "this line has one field, we will use filename as samplename"
	             sample=$( echo "$SampleLine" )
	       fi

               RealignOutputDir=$outputdir/$sample/realign
               VcallOutputDir=$outputdir/${sample}/variant  
               RealignLog=$RealignOutputDir/logs
               
               ### launching realrecal jobs
               realrecal_job=`qsub $RealignOutputLogs/qsub.realrecal.${sample}.AnisimovLauncher`
               `qhold -h u $realrecal_job`
               echo $realrecal_job >> $RealignOutputLogs/REALRECALpbs

      
               echo "####################################################################################################"
               echo "###############     constructing the qsub for verifysample   and launching it   ####################"
               echo "####################################################################################################"
  

               qsub_verifySample_anisimov=$RealignOutputLogs/qsub.verifySample.${sample}.AnisimovLauncher
               cat $outputdir/qsubGenericHeader > $qsub_verifySample_anisimov
	       echo "#PBS -N ${pipeid}_verifySample_${sample}" >> $qsub_verifySample_anisimov
	       echo "#PBS -l walltime=$pbscpu" >> $qsub_verifySample_anisimov
	       echo "#PBS -o $RealignLog/log.verifySample.${sample}.ou" >> $qsub_verifySample_anisimov
	       echo "#PBS -e $RealignLog/log.verifySample.${sample}.in" >> $qsub_verifySample_anisimov
	       echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_verifySample_anisimov
               echo "#PBS -W depend=afterok:$realrecal_job" >> $qsub_verifySample_anisimov
               echo "aprun -n 1 -N 1 -d $thr $scriptdir/verifySample.sh $runfile $sample $RealignOutputDir  ${sample}.verified $outputdir $RealignLog/log.verifySample.$sample.in $RealignLog/log.verifySample.$sample.ou $email $qsub_verifySample_anisimov" >> $qsub_verifySample_anisimov
               verifysample_job=`qsub $RealignOutputLogs/qsub.verifySample.${sample}.AnisimovLauncher` 
               echo $verifysample_job>> $RealignOutputLogs/VERIFYXSAMPLEpbs


               # now launchng vcallgatk
      
               if [ $skipvcall == "NO" ]
               then
                     sed -i "2i #PBS -W depend=afterok:$verifysample_job" $VcallOutputLogs/qsub.vcalgatk.${sample}.AnisimovLauncher
                     vcallgatk_job=`qsub $VcallOutputLogs/qsub.vcalgatk.${sample}.AnisimovLauncher`
                     echo $vcallgatk_job >> $VcallOutputLogs/VCALGATKpbs
               fi
          fi # done processing non-empty line    
      done <  $TheInputFile

   ;;
   esac


   echo -e "\n\n\n ####################   wrap up and produce summary table #########################\n\n\n"

   if [ $skipvcall == "NO" ]
   then
      summarydependids=$( cat $VcallOutputLogs/VCALGATKpbs | sed "s/\..*//" | tr "\n" ":" )
   else
      summarydependids=$( cat $RealignOutputLogs/VERIFYXSAMPLEpbs | sed "s/\..*//" | tr "\n" ":" )
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
   echo "$scriptdir/summary.sh $runfile $email exitok $reportticket"  >> $qsub_summary
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
       echo "$scriptdir/summary.sh $runfile $email exitnotok $reportticket"  >> $qsub_summary
       `chmod a+r $qsub_summary`
       badjobid=`qsub $qsub_summary`
       echo $badjobid >> $TopOutputLogs/SUMMARYpbs
   fi


   # release all jobs now

   realrecalids=$( cat $RealignOutputLogs/REALRECALpbs | sed "s/\..*//" | tr "\n" " " )
   `qrls -h u $realrecalids`

   if [ $skipvcall == "NO" ]
   then
      vcalids=$( cat $VcallOutputLogs/VCALGATKpbs | sed "s/\..*//" | tr "\n" " " )
      `qrls -h u $vcalids`
   fi


fi
