#!/bin/bash
#
#  script to realign and recalibrate the aligned file(s)
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 7 ]
then
   MSG="parameter mismatch."
   echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
   exit 1;
else
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
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi



   echo "####################################################################################################"
   echo "####################################################################################################"
   echo "#####      SCRIPT=scriptfile        PARSING PARAMETERS AND SANITY CHECK                       ######"
   echo "####################################################################################################"
   echo "####################################################################################################"


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
   indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
   realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
   sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
   chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
   indices=$( echo $chrindex | sed 's/:/ /g' )
   javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
   shareREFGENOMEmode="YES"

   echo "#################################################################################"
   echo -e cheching type of analysis
   echo "#################################################################################"

   if [ $analysis != "TESTING" ]
   then
      MSG="ANALYSIS=$analysis Program=$scriptfile Invalid option for this type of analysis. This program is for the TESTING case only"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
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
         echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
         exit 1;
      fi
   fi

   echo -e "#################################################################################"
   echo -e "####  checking sample info file. For this test, it is a two-column file      ####"
   echo -e "####  col1: sample_id        col2: full-path-to-deduplicated.bam             ####"
   echo -e "#################################################################################"

   if [ ! -s $sampleinfo ]
   then
      MSG="SAMPLEINFORMATION=$sampleinfo file not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"    
      exit 1;
   fi

   echo "#################################################################################"
   echo -e java for gatk
   echo "#################################################################################"

   if [ -z $javamodule ]
   then
      MSG="Value for JAVAMODULE must be specified in configuration file"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi

   echo "#################################################################################"
   echo -e directories for output and for tools
   echo "#################################################################################"

   if [ ! -d $outputdir ]
   then
      MSG="$outputdir ROOT directory for this run of the pipeline not found. Creating it"
      mkdir -p $outputdir  
   fi

   if [ ! -d $picardir ]
   then
      MSG="$picardir picard directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi
   if [ ! -d $samdir ]
   then
      MSG="$samdir samtools directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi
   if [ ! -d $gatk ]
   then
      MSG="$gatk GATK directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi

   if [ ! -d $refdir ]
   then
      MSG="$refdir reference genome directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi
   if [ ! -s $refdir/$ref ]
   then
      MSG="$ref reference genome not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi
   
   if [ ! -d $refdir/$indeldir ]
   then
      MSG="$indeldir indel directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi


   echo "####################################################################################################"
   echo "####################################################################################################"
   echo "#####      parameters ok. Now creating log folders  and populating them                       ######"
   echo "####################################################################################################"
   echo "####################################################################################################"
  

   TopOutputLogs=$outputdir/logs
   RealignOutputLogs=$outputdir/logs/realign

   if [ ! -d $TopOutputLogs ]
   then
      MSG="$TopOutputLogs directory not found. Creating it"
      mkdir -p $TopOutputLogs
   fi
   if [ ! -d $RealignOutputLogs ]
   then
      mkdir -p $RealignOutputLogs
   fi


   pipeid=$( cat $TopOutputLogs/CONFIGUREpbs )

   if [ $run_method == "LAUNCHER" ]
   then
      truncate -s 0 $RealignOutputLogs/realrecal.AnisimovJoblist
   else
      MSG="RUNMETHOD=$run_method Invalid option for this type of analysis. This program is for the LAUNCHER case only"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;   
   fi

   cp $runfile $outputdir/runfile.txt
   runfile=$outputdir/runfile.txt   
   cp $sampleinfo $outputdir/SAMPLEBAMNAMES.list
   TheInputFile=$outputdir/SAMPLEBAMNAMES.list 
   
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


   echo -e "\n####################################################################################################"              
   echo -e "\n################ checking that aligned files exist. Creating output folders          ###############"
   echo -e "\n"############### ALso creating RGPAMS array with the RG line for all samples         ###############"     
   echo -e "\n####################################################################################################"

   i=1   # this time this index will iterate over samples

   while read SampleLine
   do
      if [ `expr ${#SampleLine}` -lt 1 ]
      then
          echo -e "######################################################################"
          echo -e "##########      skipping empty line                        ###########"
          echo -e "######################################################################"          
      else    
          echo -e "######################################################################"
          echo -e "##########      parsing line: $SampleLine                  ###########"
          echo -e "######################################################################" 
          
          sample=$( cat $SampleLine | cut -f 1 )
          alignedfile=$( cat $SampleLine | cut -f 2 )
          
          echo -e "######################################################################"
	  echo -e "########## let us checking that alignment info exists        #########"
          echo -e "######################################################################"

	  if [ -s $alignedfile ]
	  then
              echo -e "alignment file for this sample $sample were found at $alignedfile "
          else
              MSG="No aligned bam file found at $alignedfile"
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
              exit 1;              
	  fi
	  
          echo -e "######################################################################"
          echo -e "##### now we can create the folders for the rest of the analysis  ####"
          echo -e "######################################################################"
	  if [ -d $outputdir/${sample} ]
	  then
              echo -e "creating output folders for sample=$sample"
              mkdir -p $outputdir/${sample}/realign/logs
          fi

          # resetting logs
	  rm $outputdir/${sample}/realign/logs/*

          echo -e "######################################################################"
          echo -e "##### now forming the RG line for this sample                     ####"
          echo -e "######################################################################"
 

          RGline=${alignedfile}.RGline
          if [ `expr ${#RGline}` -lt 1 -o ! -s $RGline ]
          then
                 echo -e "RGparms line needs to be recreated from scratch" 
                 RGPARMS[$i]=$( grep "^@RG" ${alignedfile}.header | sed 's/@RG//' | tr ":" "=" | tr " " ":" | tr "\t" ":" | sed "s/ //g" )
          else 
                 RGPARMS[$i]=$( cat $RGline )
          fi	  
          (( i++ ))
          
      fi  # end processing non-empty lines
   done  <  $TheInputFile    # end loop over samples

   echo -e "\n####################################################################################################"
   echo -e "\n##########                      MAIN BLOCK START HERE                                     ###########"
   echo -e "\n####################################################################################################"

   echo -e "\n####################################################################################################"
   echo -e "\n##########            OUTER LOOP by chromosome  starts here                              ###########"
   echo -e "\n####################################################################################################"
	      
   chromosomecounter=1
   for chr in $indices
   do
       echo -e "\n ################################################################################################\n"
       echo -e "\n ######    generating real-recal, vcallgatk  calls for chr=${chr}                          ######\n"
       echo -e "\n ################################################################################################\n"      
       echo `date`

       truncate -s 0 $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist
		  
       echo -e "\n#################################################################################################"
       echo -e "\n##########        INNER LOOP by sample  starts here                                      ########"
       echo -e "\n#################################################################################################"			  
       samplecounter=1 
       while read SampleLine
       do
          if [ `expr ${#SampleLine}` -lt 1 ]
          then
               echo -e "######################################################################"
               echo -e "##########      skipping empty line                        ###########"
               echo -e "######################################################################"                      
	  else
               echo -e "######################################################################"
               echo -e "##########      parsing line: $SampleLine                  ###########"
               echo -e "######################################################################" 
          
               sample=$( cat $SampleLine | cut -f 1 )
               alignedfile=$( cat $SampleLine | cut -f 2 )	

               echo -e "######################################################################"
               echo -e "### forming the realrecal command and populating anisimov launcher  ##"
               echo -e "######################################################################" 
               
               RealignOutputDir=$outputdir/${sample}/realign
               realrecal_outputfile=${chr}.${sample}.realrecal.bam
               realrecal_cmd=$RealignOutputDir/logs/realrecal.${sample}.${chr}  

               if [ $shareREFGENOMEmode  == "YES" ]
               then
                    echo "$scriptdir/realrecal_sharedRefGenomeMode.sh $RealignOutputDir $realrecal_outputfile $chr $alignedfile ${RGPARMS[$samplecounter]} ${region[$chromosomecounter]} ${realparms[$chromosomecounter]} ${recalparms[$chromosomecounter]} $runfile $flag $RealignOutputDir/logs/log.realrecal.$sample.$chr.in $RealignOutputDir/logs/log.realrecal.$sample.$chr.ou $email $RealignOutputDir/logs/realrecal.${sample}.${chr}" > $RealignOutputDir/logs/realrecal.${sample}.${chr}
               else
                    echo "$scriptdir/realrecal_NONsharedRefGenomeMode.sh $RealignOutputDir $realrecal_outputfile $chr $alignedfile ${RGPARMS[$samplecounter]} ${region[$chromosomecounter]} ${realparms[$chromosomecounter]} ${recalparms[$chromosomecounter]} $runfile $flag $RealignOutputDir/logs/log.realrecal.$sample.$chr.in $RealignOutputDir/logs/log.realrecal.$sample.$chr.ou $email $RealignOutputDir/logs/realrecal.${sample}.${chr}" > $RealignOutputDir/logs/realrecal.${sample}.${chr}
               fi
               awk -v awkvar_realrecallog=$realrecal_cmd '{print "nohup "$0" > "awkvar_realrecallog}' $RealignOutputDir/logs/realrecal.${SampleName}.${chr} > $RealignOutputDir/logs/jobfile.realrecal.${SampleName}.${chr}
               echo "$RealignOutputDir/logs/ jobfile.realrecal.${SampleName}.${chr}" >> $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist
               
               
               (( samplecounter++ ))
          fi # dome processing non-empty line
       done # done INNER LOOP going through samples       
       ((  chromosomecounter++ ))               
   done <  $TheInputFile   # done OUTER LOOP going through chromosomes

   # at the end of this set of nested loops, the variables chromosomecounter and samplecounter
   # reflect, respectively, number_of_chromosomes+1 and number_of_samples+1,
   # which is exactly the number of nodes required for anisimov launcher 


   echo -e "\n####################################################################################################"
   echo -e "\n##########                      MAIN BLOCK ENDS HERE                                     ###########"
   echo -e "\n####################################################################################################"


   echo -e "\n####################################################################################################"
   echo -e "\n##########      NEXT LOOP TO FORM AND LAUNCH QSUB JOBS FOR ANISIMOV LAUNCHERS            ###########"
   echo -e "\n##########      ONE LAUNCHER PER CHROMOSOME. EACH LAUNCHER HAS THE SAME NUMBER OF JOBS   ###########"  
   echo -e "\n####################################################################################################"

   truncate -s 0 $RealignOutputLogs/REALRECAL.pbs

   for chr in $indices
   do
 
       qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher   
       cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov

       echo "#PBS -N ${pipeid}_realrecal_${chr} " >> $qsub_realrecal_anisimov
       echo "#PBS -l walltime=$pbscpu " >> $qsub_realrecal_anisimov
       echo "#PBS -o $RealignOutputLogs/log.realrecal.${chr}.ou " >> $qsub_realrecal_anisimov
       echo "#PBS -e $RealignOutputLogs/log.realrecal.${chr}.in " >> $qsub_realrecal_anisimov
       echo "#PBS -l nodes=$samplecounter:ppn=$thr " >> $qsub_realrecal_anisimov 
       echo "aprun -n $samplecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realrecal.${chr}.AnisimovLauncher.log" >> $qsub_realrecal_anisimov
       jobid=`qsub $qsub_realrecal_anisimov`
       echo $jobid >> $RealignOutputLogs/REALRECAL.pbs

   done # loop over chromosomes 

   echo -e "\n####################################################################################################"
   echo -e "\n##########                      DONE.    EXITING NOW                                     ###########"
   echo -e "\n####################################################################################################"


