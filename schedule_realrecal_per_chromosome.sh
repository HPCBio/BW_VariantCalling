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
   echo -e "\n\n############# schedule_realrecal_per_chromosome: schedule jobs to split by chromosome, realrecal, and if necessary variant calling  ###############\n\n" >&2

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
      echo -e "Program $0 stopped. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
      exit 1;
   fi



   set +x; echo -e "\n\n" >&2;
   # wrapping commends in echo, so that the output logs would be easier to read: they will have more structure
   echo "####################################################################################################" >&2
   echo "#####################################                       ########################################" >&2
   echo "##################################### PARSING RUN INFO FILE ########################################" >&2
   echo "#####################################                       ########################################" >&2
   echo "####################################################################################################" >&2
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
   kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
   targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
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

   set +x; echo -e "\n"
   echo "#################################################################################" >&2
   echo "checking type of analysis " >&2
   echo "#################################################################################" >&2
   echo -e "\n"; set -x

   if [ $analysis == "MULTIPLEXED" ]
   then
      MSG="ANALYSIS=$analysis Program=$scriptfile Invalid pipeline program for this type of analysis. This program is for the NON-MULTIPLEXED case only"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi

   set +x; echo -e "\n"
   echo "#################################################################################" >&2
   echo "checking input type" >&2
   echo "#################################################################################" >&2;
   echo -e "\n"; set -x

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
         echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
         exit 1;
      fi
   fi


   set +x; echo -e "\n"
   echo "#################################################################################" >&2
   echo "checking cleanup options" >&2
   echo "#################################################################################" >&2
   echo -e "\n"; set -x

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


   set +x; echo -e "\n"
   echo "#################################################################################" >&2
   echo "skip/include variant calling module" >&2
   echo "#################################################################################" >&2
   echo -e "\n"; set -x

   if [ $skipvcall != "1" -a $skipvcall != "0" -a $skipvcall != "YES" -a $skipvcall != "NO" ]
   then
      MSG="Invalid value for SKIPVCALL=$skipvcall"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
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

   set +x; echo -e "\n"
   echo "#################################################################################" >&2
   echo "java for gatk" >&2
   echo "#################################################################################" >&2
   echo -e "\n"; set -x

   if [ -z $javamodule ]
   then
      MSG="Value for JAVAMODULE must be specified in configuration file"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi

   set +x; echo -e "\n"
   echo "#################################################################################" >&2
   echo "directories for output and for tools" >&2
   echo "#################################################################################" >&2
   echo -e "\n"; set -x

   if [ ! -d $outputdir ]
   then
      MSG="$outputdir ROOT directory for this run of the pipeline not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi

   if [ ! -d $picardir ]
   then
      MSG="$picardir picard directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi
   if [ ! -d $samdir ]
   then
      MSG="$samdir samtools directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi
   if [ ! -d $gatk ]
   then
      MSG="$gatk GATK directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi

   if [ ! -d $refdir ]
   then
      MSG="$refdir reference genome directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi
   if [ ! -s $refdir/$ref ]
   then
      MSG="$ref reference genome not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
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
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi


   set +x; echo -e "\n\n"
   echo "####################################################################################################" >&2
   echo "###########################      parameters ok. Now      creating log folders        ###############" >&2
   echo "####################################################################################################" >&2
   echo -e "\n"; set -x
  

   TopOutputLogs=$outputdir/logs
   RealignOutputLogs=$outputdir/logs/realign
   VcallOutputLogs=$outputdir/logs/variant
   if [ ! -d $TopOutputLogs ]
   then
      MSG="$TopOutputLogs realign directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
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








   set +x; echo -e "\n\n"
   echo "####################################################################################################" >&2
   echo "###########################   generating regions, intervals, known/knownSites        ###############" >&2
   echo "####################################################################################################" >&2
   echo -e "\n"
   echo "###### region is the array with snps per chr which will be used by vcallgatk-unifiedGenotyper/haplotypeCaller" >&2
   echo "###### realparms is the array with indels per chr which will be used for  gatk-IndelRealigner" >&2
   echo "###### recalparms is the array with indels per chr which will be used for  gatk-Recalibration" >&2
   echo -e "\n"; set -x

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

   set +x; echo -e "\n\n"
   echo "####################################################################################################" >&2
   echo "####################### done  generating regions, intervals, known/knownSites        ###############" >&2
   echo "####################################################################################################" >&2
   echo -e "\n"; set -x







   set +x; echo -e "\n\n"
   echo "#############################################################################################" >&2
   echo "##############################                              #################################" >&2
   echo "##############################  CREATE OUTPUT  DIRECTORIES  #################################" >&2
   echo "##############################                              #################################" >&2
   echo "#############################################################################################" >&2
   echo -e "\n"; set -x


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
              #this line has five fields, we just need the first field to create folders
	      sample=$( echo "$SampleLine" | cut -f 1 )
	  else
              #this line has one field, we need to parse the sample name out of the filename
	      sample=$( echo "$SampleLine" )
	  fi
	  
	  set +x; echo -e "\n ########## first, checking that alignment information exists  #########" >&2; set -x

	  alignedfile=`find $outputdir/$sample/align -name "*.wdups.sorted.bam"`
	  alignedfilehdr=`find $outputdir/$sample/align -name "*.wdups.sorted.bam.header"`
	  alignedfilebai=`find $outputdir/$sample/align -name "*.wdups.sorted.bam.bai"`

	  if [ -s $alignedfile -a -s $alignedfilehdr -a -s $alignedfilebai ]
	  then
              set +x; echo -e "\n # alignment files for this sample $sample were found at $outputdir/${sample}/align" >&2; set -x
          else
              MSG="No aligned bam file(s) found at $outputdir/${sample}/align"
              echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
              exit 1;              
	  fi

          set +x; echo -e "\n ########## now we can create the folders for the rest of the analysis  #########" >&2; set -x
	  if [ -d $outputdir/${sample} ]
	  then
              set +x; echo -e "\n # creating output folders for sample=$sample" >&2; set -x

              mkdir -p $outputdir/${sample}/realign/logs
              mkdir -p $outputdir/${sample}/variant/logs
          fi
      fi  # end processing non-empty lines
   done  <  $TheInputFile
   # end loop over samples









   set +x; echo -e "\n\n"
   echo "#############################################################################################################" >&2
   echo "########################                                                             ########################" >&2
   echo "########################                  main loops starts here                     ########################" >&2
   echo "########################   outer loops by sample; inner loops by chromosome/region   ########################" >&2
   echo "########################                                                             ########################" >&2
   echo "#############################################################################################################" >&2
   echo -e "\n"; set -x


   set +x; echo -e "\n\n"
   echo "####################################################################################################" >&2
   echo "##########                       loop by sample starts here                              ###########" >&2
   echo "####################################################################################################" >&2
   echo -e "\n"; set -x

      samplecounter=1 
      while read SampleLine
      do
          if [ `expr ${#SampleLine}` -gt 1 ]
          then
              ## processing non-empty line
	      echo `date`
	      set +x; echo -e "\n ######     processing this line: $SampleLine #####\n" >&2; set -x

              if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
              then
                   ###### this line has five fields, we just need the first field to create folders
	           sample=$( echo "$SampleLine" | cut -f 1 )
	      else
                   ###### this line has one field, we will use filename as samplename
	           sample=$( echo "$SampleLine" )
	      fi

              set +x; echo -e "\n\n"
              echo "####################################################################################################" >&2
              echo "######################################  get aligned bam files  #####################################" >&2
              echo "####################################################################################################" >&2
              echo "## we already checked that these files exist when we were creating folders in the previous block  ##" >&2
              echo "####################################################################################################" >&2
              echo -e "\n"; set -x
              

              cd $outputdir/${sample}/align
	      aligned_bam=`find ./ -name "*.wdups.sorted.bam"`
	      
              # now check that there is only one bam file
              aligned_bam=$outputdir/${sample}/align/${aligned_bam}
              aligned_bam=( $aligned_bam ) # recast variable as an array and count the number of members
              if [ ${#aligned_bam[@]} -ne 1 ]
              then
                 MSG="more than one bam file found in $outputdir/${sample}/align"
                 echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                 exit 1;
              fi

              aligned_bam="${aligned_bam[*]}" # recast variable back as string	      
              header=${aligned_bam}.header
 
              set +x; echo -e "\n\n"
              echo "####################################################################################################" >&2
              echo "##########                   now putting together the RG line for RGparms                ###########" >&2
	      echo "####################################################################################################" >&2
              echo -e "\n"; set -x
 

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

 
              set +x; echo -e "\n\n"
              echo "####################################################################################################" >&2
              echo "##########                   now defining paths to output folders                        ###########" >&2
	      echo "####################################################################################################" >&2
              echo -e "\n"; set -x
              
              RealignOutputDir=$outputdir/${sample}/realign
              VcallOutputDir=$outputdir/${sample}/variant

              set +x; echo -e "\n\n"
              echo "####################################################################################################" >&2
              echo "##########                   loop by chromosome starts here                              ###########" >&2
	      echo "####################################################################################################" >&2
              echo -e "\n"; set -x

	      chromosomecounter=1
	      for chr in $indices
	      do
                  set +x; echo -e "\n\n"
		  echo "#########################################################################################################\n" >&2
		  echo "######    generating real-recal, vcallgatk  calls for chr=${chr}                                   ######\n" >&2
		  echo "#########################################################################################################\n" >&2  
		  echo `date`
                  echo -e "\n"; set -x

                  #######################   assemble the real-recall/var call sub-block 
                  echo "$scriptdir/realrecal.sh $RealignOutputDir $chr.realrecal.$sample.output.bam $chr $aligned_bam $RGparms ${region[$chromosomecounter]} ${realparms[$chromosomecounter]} ${recalparms[$chromosomecounter]} $runfile $flag $RealignOutputDir/logs/log.realrecal.$sample.$chr.in $RealignOutputDir/logs/log.realrecal.$sample.$chr.ou $email $RealignOutputDir/logs/realrecal.${sample}.${chr}" > $RealignOutputDir/logs/realrecal.${sample}.${chr}
           
                  if [ $skipvcall == "NO" ]
                  then
		      #######################   assemble the vcallgatk call sub-block 
		      echo "$scriptdir/vcallgatk.sh $VcallOutputDir  $RealignOutputDir ${chr}.realrecal.${sample}.output.bam $chr ${region[$chromosomecounter]} $runfile $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in $VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.ou $email $VcallOutputDir/logs/vcallgatk.${sample}.${chr}" >> $VcallOutputDir/logs/vcallgatk.${sample}.${chr}
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


      set +x; echo -e "\n\n"
      echo "####################################################################################################" >&2
      echo "##########                      main loop ends here                                      ###########" >&2
      echo "####################################################################################################" >&2
      echo -e "\n"; set -x



      set +x; echo -e "\n\n"
      echo "#############################################################################################################" >&2
      echo "###################################                             #############################################" >&2
      echo "###################################   now schedule these jobs   #############################################" >&2
      echo "###################################                             #############################################" >&2
      echo "#############################################################################################################" >&2
      echo -e "\n"; set -x

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
      set +x; echo -e "\n ###### run_method=$run_method. generating and scheduling qsubs... ###### \n" >&2; set -x

      truncate -s 0 $RealignOutputLogs/SPLITBYCHROMOSOMEpbs
      truncate -s 0 $RealignOutputLogs/REALRECALpbs
      truncate -s 0 $VcallOutputLogs/VCALGATKpbs

      set +x; echo -e "\n ###### loop by sample  ###### \n" >&2; set -x
      while read SampleLine
      do

          if [ `expr ${#SampleLine}` -gt 1 ]
          then
              ## processing non-empty line
              if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
              then
                   #this line has five fields, we just need the first field to create folders
	           sample=$( echo "$SampleLine" | cut -f 1 )
	      else
                   #this line has one field, we will use filename as samplename
	           sample=$( echo "$SampleLine" )
	      fi
              RealignOutputDir=$outputdir/$sample/realign
              VcallOutputDir=$outputdir/${sample}/variant 

	      set +x; echo -e "\n ###### loop by chromosome  ###### \n" >&2; set -x
              for chr in $indices
              do

                 qsub_realrecal=$RealignOutputDir/logs/qsub.realrecal.${sample}.${chr}

                 ###############################
                 set +x; echo -e "\n################# constructing qsub for realrecal\n" >&2; set -x
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
                     set +x; echo -e "\n################# constructing qsub for vcall\n" >&2; set -x
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

      set +x; echo -e "\n"
      echo "###### run_method=$run_method. populate list of jobs first  and launching them later via launcher... ######" >&2
      echo "###### loop1 by chromosome to create  ONE list of jobs per chromosome  ######" >&2
      echo -e "\n"; set -x

      for chr in $indices
      do
         # clear out the joblists

         truncate -s 0 $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist
         if [ $skipvcall == "NO" ]
         then
            truncate -s 0 $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist
         fi
         

	 set +x; echo -e "\n ###### loop by sample to populate the list by writing ONE line for each job call to a sample-chr pair ###### \n" >&2; set -x
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
         
         
                # creating a qsub out of the job file
                # need to prepend "nohup" and append log file name, so that logs are properly created when Anisimov launches these jobs 

                realrecal_log=$RealignOutputDir/logs/log.realrecal.${sample}.$chr.in
                awk -v awkvar_realrecallog=$realrecal_log '{print "nohup "$0" > "awkvar_realrecallog}' $RealignOutputDir/logs/realrecal.${sample}.${chr} > $RealignOutputDir/logs/jobfile.realrecal.${sample}.${chr}
                echo "$RealignOutputDir/logs/ jobfile.realrecal.${sample}.${chr}" >> $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist

                if [ $skipvcall == "NO" ]
                then
                     vcall_log=$VcallOutputDir/logs/log.vcallgatk.${sample}.${chr}.in
                     awk -v awkvar_vcalllog=$vcall_log '{print "nohup "$0" > "awkvar_vcalllog}' $VcallOutputDir/logs/vcallgatk.${sample}.${chr} > $VcallOutputDir/logs/jobfile.vcallgatk.${sample}.${chr}
                     echo "$VcallOutputDir/logs/ jobfile.vcallgatk.${sample}.${chr}" >> $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist
                fi
              fi # non-empty line in file
         done <  $TheInputFile
         # end loop over samples

         set +x; echo -e "\n";
	 echo "###### finished creating the joblist for Launcher             ######" >&2
         echo "###### now putting together the other pieces of the qsub file ######" >&2
         echo "###### and then scheduling Anisimov Launcher joblists         ######" >&2
         echo -e "\n"; set -x

         qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher
         # appending the generic header to the qsub

         cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov

         if [ $skipvcall == "NO" ]
         then
            qsub_vcallgatk_anisimov=$VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher
            cat $outputdir/qsubGenericHeader > $qsub_vcallgatk_anisimov

         fi




         ###############
         ############### constructing the qsub for realrecal
         echo "#PBS -N ${pipeid}_realrecal_${chr}" >> $qsub_realrecal_anisimov
         echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal_anisimov
         echo "#PBS -o $RealignOutputLogs/log.realrecal.${chr}.ou" >> $qsub_realrecal_anisimov
         echo -e "#PBS -e $RealignOutputLogs/log.realrecal.${chr}.in\n" >> $qsub_realrecal_anisimov
         # the dependency on split_bam_by_chromosome job will be added when it is scheduled in the loop below

         # realrecal and vcall actually use multithreaded processes,
         # so we will give each sample its own node
         # +1 for the launcher
         # samplecounter is already more than actual number of samples by 1
         echo -e "#PBS -l nodes=$samplecounter:ppn=$thr\n" >> $qsub_realrecal_anisimov
         # the actual command
         echo "aprun -n $samplecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realrecal.${chr}.AnisimovLauncher.log" >> $qsub_realrecal_anisimov



         if [ $skipvcall == "NO" ]
         then
            ###############
            ############### constructing the qsub for vcallgatk
            echo "#PBS -N ${pipeid}_vcallgatk_${chr}" >> $qsub_vcallgatk_anisimov
            echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk_anisimov
            echo "#PBS -o $RealignOutputLogs/log.vcallgatk.${chr}.ou" >> $qsub_vcallgatk_anisimov
            echo -e "#PBS -e $RealignOutputLogs/log.vcallgatk.${chr}.in\n" >> $qsub_vcallgatk_anisimov
            # the dependency on realrecal job will be added when it is scheduled in the loop below

            # realrecal and vcall actually use multithreaded processes, 
            # so we will give each sample its own node
            # +1 for the launcher
            # samplecounter is already more than actual number of samples by 1
            echo -e "#PBS -l nodes=$samplecounter:ppn=$thr\n" >> $qsub_vcallgatk_anisimov
            # the actual command
            echo "aprun -n $samplecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist /bin/bash > $VcallOutputLogs/vcallgatk.${chr}.AnisimovLauncher.log" >> $qsub_vcallgatk_anisimov
        fi

      done # done going through chromosomes



      ###############
      ############### 
      set +x; echo -e "\n\n"
      echo "###### Now going through chromosomes again to schedule them in the right order " >&2
      echo "###### all split jobs are before all realrecal jobs, and all realrecal jobs are before vcallgatk jobs " >&2
      echo "###### this will efficiently work with a 25 job limit on queued state " >&2
      echo -e "\n"; set -x;

      # reset the lists
      
      truncate -s 0 $RealignOutputLogs/REALRECALpbs
      if [ $skipvcall == "NO" ]
      then
         truncate -s 0 $VcallOutputLogs/VCALGATKpbs
      fi
      cd $RealignOutputLogs # so that whatever temp fioles and pbs notifications would go there

      # now realrecal
      for chr in $indices
      do
         realrecal_job=`qsub $RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher`
         # not putting these on hold, as they will be helkd back by the dependency on respective split jobs
         # Add dependency on realrecal job to vcallgatk job
         if [ $skipvcall == "NO" ]
         then
             sed -i "2i #PBS -W depend=afterok:$realrecal_job" $VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher
	 fi
         echo $realrecal_job >> $RealignOutputLogs/REALRECALpbs
      done
      if [ $skipvcall == "NO" ]
      then
         # finally submit the vcallgatk
         cd $VcallOutputLogs # so that whatever temp files and pbs notifications would go there
         for chr in $indices
         do
            vcallgatk_job=`qsub $VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher`
            echo $vcallgatk_job >> $VcallOutputLogs/VCALGATKpbs
         done 
      fi


   ;;
   esac


   set +x; echo -e "\n\n\n ####################   wrap up and produce summary table #########################\n\n\n" >&2; set -x

   if [ $skipvcall == "NO" ]
   then
      summarydependids=$( cat $VcallOutputLogs/VCALGATKpbs | sed "s/\..*//" | tr "\n" ":" )
   else
      summarydependids=$( cat $RealignOutputLogs/REALRECALpbs | sed "s/\..*//" | tr "\n" ":" )
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
