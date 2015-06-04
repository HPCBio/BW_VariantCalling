#!/bin/bash
# written in collaboration with Mayo Bioinformatics core group
#
#  script to realign and recalibrate the aligned file(s)
#redmine=hpcbio-redmine@igb.illinois.edu
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

   if [ -z $javamodule ]
   then
      MSG="Value for JAVAMODULE must be specified in configuration file"
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

   if [ -s $refdir/$dbSNP ]
   then
      realparms="-known:$refdir/$dbSNP"
      recalparms="--knownSites:$refdir/$dbSNP"
   else
      MSG="dbSNP not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi
   if [ -s $refdir/$kgenome ]
   then
      realparms=$realparms":-known:$refdir/$kgenome"
      recalparms=$recalparms":--knownSites:$refdir/$kgenome"
   fi
   #if [ ! -d $RealignOutputLogs ]
   #then
   #   MSG="$RealignOutputLogs realignlog directory not found"
   #   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
   #   exit 1;
   #fi


   RealignOutputDir=$outputdir/realign
   if [ ! -d $RealignOutputDir ]
   then
      MSG="$RealignOutputDir realign directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

   if [ $skipvcall == "NO" ]
   then
      VcallOutputDir=$outputdir/variant
      if [ ! -d $VcallOutputDir ]
      then
         echo "$VcallOutputDir variant directory not found, creating it"
         mkdir $VcallOutputDir
      fi
   fi


   TopOutputLogs=$outputdir/logs
   if [ ! -d $TopOutputLogs ]
   then
      MSG="$TopOutputLogs realign directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
      exit 1;
   fi

   RealignOutputLogs=$outputdir/logs/realign
   if [ ! -d $RealignOutputLogs ]
   then
      MSG="$RealignOutputLogs realign directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi

   if [ $skipvcall == "NO" ]
   then
      VcallOutputLogs=$outputdir/logs/variant
      if [ ! -d $VcallOutputLogs ]
      then
         MSG="$VcallOutputLogs realign directory not found"
         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
         #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeli
         exit 1;
      fi
   fi


   if [ $run_method != "LAUNCHER" -a $run_method != "QSUB" -a $run_method != "APRUN" -a $run_method != "SERVER" ]
   then
      MSG="Invalid value for RUNMETHOD=$run_method"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
      exit 1;
   elif [ $schedule == "LAUNCHER" ]
      then
      truncate -s 0 $RealignOutputLogs/$realrecal.AnisimovJoblist
   fi


   pipeid=$( cat $TopOutputLogs/CONFIGUREpbs )






                ####################################################################################################
                #####################################                       ########################################
echo -e "\n\n\n #####################################  CREATE  DIRECTORIES  ######################################## \n\n\n"
                #####################################                       ########################################
                ####################################################################################################



   while read SampleName
   do
      if [ ! -d $RealignOutputDir/${SampleName} ]
      then
         mkdir -p $RealignOutputDir/${SampleName}
      else
         echo "$RealignOutputDir/${SampleName} already exists. resetting it"
         `rm -r $RealignOutputDir/${SampleName}/*`
      fi

      #if [ $multisample == "NO" ]
      #then

         # when running multiple INDEPENDENT samples, there will be lots of qsubs, logs and jobfiles
         # best keep them in the sample subfolder
         if [ ! -d $RealignOutputDir/${SampleName}/logs ]
         then
            mkdir $RealignOutputDir/${SampleName}/logs
         else
            echo "$RealignOutputDir/${SampleName}/logs already exists. resetting it"
            `rm -r $RealignOutputDir/${SampleName}/logs/*`
         fi

         if [ $skipvcall == "NO" ]
         then
            if [ ! -d $VcallOutputDir/${SampleName} ]
            then
               mkdir -p $VcallOutputDir/${SampleName}
            else 
               echo "$VcallOutputDir/${SampleName} already exists. resetting it"
               `rm -r $VcallOutputDir/${SampleName}/*`
            fi

            # when running multiple samples, there will be lots of qsubs, logs and jobfiles
            # best keep them in the sample subfolder
            if [ ! -d $VcallOutputDir/${SampleName}/logs ]
            then
               mkdir $VcallOutputDir/${SampleName}/logs
            else
               echo "$VcallOutputDir/${SampleName}/logs already exists. resetting it"
               `rm -r $VcallOutputDir/${SampleName}/logs/*`
            fi
         fi
      #fi
   done  <  $outputdir/SAMPLENAMES.list
   # end loop over samples








                ###########################################################################################################
                ###########################                                                          ######################"
echo -e "\n\n\n ###########################   generating regions and intervals files in BED format   ############### \n\n\n"
                ###########################                                                          ######################
                ###########################################################################################################




   echo `date`
   i=1
   for chr in $indices
   do
      #i=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
      if [ -d $targetkit ]
      then
         if [ `cat $targetkit/${chr}.bed | wc -l` -gt 0 ]
         then
            region[$i]="-L:$targetkit/${chr}.bed"
         else
            region[$i]="-L:$chr"
         fi
      else
         region[$i]="-L:$chr"
      fi
      (( i++ ))
   done
   echo `date`

   igv_files=""
   vcf_files=""
   sep=":"







                #############################################################################################################
                ###################################                                ##########################################
echo -e "\n\n\n ###################################       main loops start here    ###################################### 
                ########################   outer loops by chromosome; inner loops by sample   ######################## \n\n\n"
                ########################                                                      ###############################
                #############################################################################################################



   chromosomecounter=1
   for chr in $indices
   do
      echo "generating real-recal calls for chr=${chr} ..."
      echo `date`
      #inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )

      samplecounter=1 
      while read SampleName
      do
         echo "processing next sample"
         echo "realigning: $SampleName"

         echo -e "\n######################################  get aligned bam files  #####################################\n"

         cd $outputdir/align/${SampleName}
         # look for aligned, sorted bams with duplicates marked, inside the align folder for that particular sample
         # there should be only one such file
         aligned_bam=`find ./ -name "*.wdups.sorted.bam"`
         #  the ${# tells the length of the string
         if [ `expr ${#aligned_bam}` -lt 1 ]
         then
            MSG="No bam file(s) found to perform realign-recalibrate at $outputdir/align/${SampleName}"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            exit 1;
         fi

         header=${aligned_bam}.header

         if [ ! -s $header ]
         then
            MSG="$header file not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            exit 1;

         fi
         truncate -s 0 ${aligned_bam}.RGline
	 `grep "RG.*ID" $header | cut -f 2` >> ${aligned_bam}.RGline
         RGline=${aligned_bam}.RGline
         RGline=$( sed -i 's/ID/RGID/' "$RGline" )
         RGline=$( sed -i 's/ / RG/g' "$RGline" )
         RGparms=$( cat "$RGline" )


         # now check that there is only one file
         aligned_bam=$outputdir/align/${SampleName}/${aligned_bam}
         aligned_bam=( $aligned_bam ) # recast variable as an array and count the number of members
         if [ ${#aligned_bam[@]} -ne 1 ]
         then
            MSG="more than one bam file found in $outputdir/align/${SampleName}"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            exit 1;
         fi

         aligned_bam="${aligned_bam[*]}" # recast variable back as string



         echo -e "\n######################################  assemble the split_bam_by_chromosome sub-block #####################################\n"
         echo "$scriptdir/split_bam_by_chromosome.sh $runfile $picardir $samdir $javamodule $RealignOutputDir/$SampleName $aligned_bam ${SampleName}.$chr.bam ${SampleName}.$chr.sorted.bam $RGparms $flag $chr $RealignOutputDir/${SampleName}/logs/log.split.$SampleName.$chr.in $RealignOutputDir/${SampleName}/logs/log.split.$SampleName.$chr.ou $email $RealignOutputDir/${SampleName}/logs/split_bam_by_chromosome.${SampleName}.${chr} $RealignOutputLogs" > $RealignOutputDir/${SampleName}/logs/split_bam_by_chromosome.${SampleName}.${chr}

         # construct the list of all input files for the GATK real-recal procedure, multisample case
         #if [ $multisample == "YES" ]
         #then
         #   chrinfiles[$chromosomecounter]=${chrinfiles[$chromosomecounter]}":-I:$RealignOutputDir/${SampleName}/${SampleName}.$chr.sorted.bam"
# not used anymore?          chrinputfiles[$inx]=${chrinputfiles[$inx]}":INPUT=$outputdir/$SampleName.$chr.sorted.bam"
         #fi
         if [ $multisample == "NO" ]
         then         
            echo -e "\n#######################   assemble the real-recall/var call sub-block for INDEPENDENT SAMPLES    ################################\n"
            echo "$scriptdir/realrecal.sh $RealignOutputDir/${SampleName} $chr.realrecal.$SampleName.output.bam $chr -I:$RealignOutputDir/${SampleName}/${SampleName}.$chr.sorted.bam ${region[$chromosomecounter]} $realparms $recalparms $runfile $flag $RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.in $RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.ou $email $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr}" > $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr}
           
            if [ $skipvcall == "NO" ]
            then
               echo "$scriptdir/vcallgatk.sh $VcallOutputDir/${SampleName}  $RealignOutputDir/${SampleName} ${chr}.realrecal.${SampleName}.output.bam $chr ${region[$chromosomecounter]} $runfile $VcallOutputDir/${SampleName}/logs/log.vcallgatk.${SampleName}.${chr}.in $VcallOutputDir/${SampleName}/logs/log.vcallgatk.${SampleName}.${chr}.ou $email $VcallOutputDir/${SampleName}/logs/vcallgatk.${SampleName}.${chr}" >> $VcallOutputDir/${SampleName}/logs/vcallgatk.${SampleName}.${chr}
            fi
         fi

         (( samplecounter++ ))

      done <  $outputdir/SAMPLENAMES.list
      # end loop over samples
      echo `date`
      


#      if [ $multisample == "YES" ]
#      then      
#         echo "##"
#         echo "###################    assemble the split_bam_by_chromosome sub-block for MULTISAMPLE experiment    ############################"
#         echo "##" # this line makes output logs more readable
###############################################################################
###############################################################################
#########################   NOT CLEAR HOW TO HANDLE BAMFILE AND rg VARIABLES HERE IN THE MULTISAMPLE CASE
#######################3          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
#####################################################3333
         #echo "$scriptdir/realrecal.sh $RealignOutputDir $chr.realrecal.output.bam $chr ${chrinfiles[$chromosomecounter]} ${region[$chromosomecounter]} $realparms $recalparms $runfile $flag $RealignOutputLogs/log.realrecal.$bamfile.$chr.in $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou $email $RealignOutputLogs/qsub.realrecal.$bamfile.$chr" >> $qsub_realrecal
#      fi


      ((  chromosomecounter++ )) 
   done # done going through chromosomes 
   # at the end of this set of nested loops, the variables chromosomecounter and samplecounter
   # reflect, respectively, number_of_chromosomes+1 and number_of_samples+1,
   # which is exactly the number of nodes required for anisimov launcher 








#         qsub_split_bam_by_chromosome=$RealignOutputLogs/qsub.split.${SampleName}.$chr
#         echo "#PBS -N ${pipeid}_split_${SampleName}_$chr" >> $qsub_split_bam_by_chromosome
#         echo "#PBS -l walltime=$pbscpu" >> $qsub_split_bam_by_chromosome
#         echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_split_bam_by_chromosome
#         echo "#PBS -o $RealignOutputDir/${SampleName}/logs/log.split.${SampleName}.$chr.ou" >> $qsub_split_bam_by_chromosome
#         echo "#PBS -e $RealignOutputDir/${SampleName}/logs/log.split.${SampleName}.$chr.in" >> $qsub_split_bam_by_chromosome
# 
#         `chmod a+r $qsub_split_bam_by_chromosome`
#            splitjobid=`qsub $qsub_split_bam_by_chromosome`
#            # new line to avoid hiccup
#            #`qhold -h u $splitjobid`
#            echo $splitjobid >> $outputdir/logs/REALSORTEDpbs.$SampleName
#            echo $splitjobid >> $outputdir/logs/REALSORTED.$SampleName$chr
#
#            truncate -s 0 $realrecal.$SampleName.$chr.serialjobs
#            echo "$RealignOutputLogs realrecal.$SampleName.$chr.serialjobs" >> $RealignOutputLogs/$realrecal.AnisimovJoblist
#






#         if [ $multisample == "NO" ]
#         then


#      splitid=$( cat $outputdir/logs/REALSORTED.$bamfile_$chr | sed "s/\..*//g" | tr "\n" ":" )
#      outputfile=$chr.realrecal.$bamfile.output.bam
#      igv_files=${igv_files}":INPUT=${outputfile}"
#      echo "realign-recalibrate for interval:$chr..."
#      qsub_realrecal=$RealignOutputLogs/qsub.realrecal.$bamfile.$chr
#      echo "#PBS -V" > $qsub_realrecal
#      echo "#PBS -A $pbsprj" >> $qsub_realrecal
#      echo "#PBS -N ${pipeid}_realrecal_$bamfile.$chr" >> $qsub_realrecal
#      echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal
#      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_realrecal
#      echo "#PBS -o $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou" >> $qsub_realrecal
#      echo "#PBS -e $RealignOutputLogs/log.realrecal.$bamfile.$chr.in" >> $qsub_realrecal
#      echo "#PBS -q $pbsqueue" >> $qsub_realrecal
#      echo "#PBS -m a" >> $qsub_realrecal
#      echo "#PBS -M $email" >> $qsub_realrecal
#      echo "#PBS -W depend=afterok:$splitid" >> $qsub_realrecal
#
#      if [ $schedule eq "QSUB" ]
#      then
#         echo "aprun -n 1 -d $thr $scriptdir/realrecal.sh $????#####outputdir $outputfile $chr ${chrinfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $RealignOutputLogs/log.realrecal.$bamfile.$chr.in $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou $email $RealignOutputLogs/qsub.realrecal.$bamfile.$chr" >> $qsub_realrecal
#         `chmod a+r $qsub_realrecal`
#         #recaljobid=`qsub $qsub_realrecal`
#         # new line to avoid hiccup
#         #`qhold -h u $recaljobid`
#         echo $recaljobid >> $outputdir/logs/REALRECALpbs.$bamfile
#      else 
#         # this block constructs the list of jobs to perform in serial on that chromosome, on that sample
#         echo "nohup $scriptdir/realrecal.sh $????#####outputdir $outputfile $chr ${chrinfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $RealignOutputLogs/log.realrecal.$bamfile.$chr.in $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou $email $RealignOutputLogs/qsub.realrecal.$bamfile.$chr > $RealignOutputLogs/log.realrecal.$bamfile.$chr.in" >> $RealignOutputLogs/realrecal.$bamfile.$chr.serialjobs
#      fi

#      if [ $skipvcall == "NO" ]
#      then
#         vardir=$( echo $????#####outputdir | sed 's/realign/variant/' )
#         varlogdir=$( echo $RealignOutputLogs | sed 's/realign/variant/' )
#  
#         vcf_files=${vcf_files}":${outputfile}"
#         echo "variant calling call.."
#         qsub_vcallgatk=$varlogdir/qsub.vcallgatk.$bamfile.$chr
#         echo "#PBS -V" > $qsub_vcallgatk
#         echo "#PBS -A $pbsprj" >> $qsub_vcallgatk
#         echo "#PBS -N ${pipeid}_vcall_$chr" >> $qsub_vcallgatk
#         echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk
#         echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_vcallgatk
#         echo "#PBS -o $varlogdir/log.vcallgatk.$bamfile.$chr.ou" >> $qsub_vcallgatk
#         echo "#PBS -e $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $qsub_vcallgatk
#         echo "#PBS -q $pbsqueue" >> $qsub_vcallgatk
#         echo "#PBS -m a" >> $qsub_vcallgatk
#         echo "#PBS -M $email" >> $qsub_vcallgatk
#         echo "#PBS -W depend=afterok:$recaljobid" >> $qsub_vcallgatk
#         if [ $schedule eq "QSUB" ]
#         then
#            echo "aprun -n 1 -d $thr $scriptdir/vcallgatk.sh $vardir $????#####outputdir $outputfile $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$bamfile.$chr.in $varlogdir/log.vcallgatk.$bamfile.$chr.ou $email $varlogdir/qsub.vcallgatk.$bamfile.$chr" >> $qsub_vcallgatk
#            `chmod a+r $qsub_vcallgatk`
##            vcalljobid=`qsub $qsub_vcallgatk`
#            echo $vcalljobid >> $outputdir/logs/VCALLGATKpbs.$bamfile
#         else
#            # something is messed up with redirection here
#            echo "nohup $scriptdir/vcallgatk.sh $vardir $????#####outputdir $outputfile $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$bamfile.$chr.in $varlogdir/log.vcallgatk.$bamfile.$chr.ou $email $varlogdir/qsub.vcallgatk.$bamfile.$chr >> $qsub_vcallgatk > $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $RealignOutputLogs/realrecal.$bamfile.$chr.serialjobs
#         fi
#      else
#         echo "variant calling will not be run"
#      fi
#
#
#      echo `date`
#      echo "bottom of the loop over chromosomes"
#   done
#   echo `date`
#
#
#     
#








                #############################################################################################################
                ###################################                             #############################################
echo -e "\n\n\n ###################################   now schedule these jobs   ###################################### \n\n\n"
                ###################################                             #############################################
                #############################################################################################################





   case $run_method in
   "APRUN")
      echo -e "\n\nscheduling qsubs\n\n"
      truncate -s 0 $RealignOutputLogs/SPLITBYCHROMOSOMEpbs
      truncate -s 0 $RealignOutputLogs/REALRECALpbs
      truncate -s 0 $VcallOutputLogs/VCALGATKpbs

      while read SampleName
      do

         for chr in $indices
         do

            qsub_split=$RealignOutputDir/${SampleName}/logs/qsub.split_bam_by_chromosome.${SampleName}.${chr}
            qsub_realrecal=$RealignOutputDir/${SampleName}/logs/qsub.realrecal.${SampleName}.${chr}


            ###############################
            echo -e "\n################# constructing qsub for split\n"
            # appending the generic header to the qsub
            cat $outputdir/qsubGenericHeader > $qsub_split
            echo "#PBS -N ${pipeid}_split_bam.${SampleName}.${chr}" >> $qsub_split
            echo "#PBS -l walltime=$pbscpu" >> $qsub_split
            echo "#PBS -o $RealignOutputDir/${SampleName}/logs/log.split_bam_by_chromosome.${SampleName}.${chr}.ou" >> $qsub_split
            echo "#PBS -e $RealignOutputDir/${SampleName}/logs/log.split_bam_by_chromosome.${SampleName}.${chr}.in" >> $qsub_split
            echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_split

            echo "aprun -n 1 -N 1 -d $thr /bin/bash $RealignOutputDir/${SampleName}/logs/split_bam_by_chromosome.${SampleName}.${chr}" >> $qsub_split

            split_job=`qsub $qsub_split`
            `qhold -h u $split_job`
            echo $split_job >> $RealignOutputLogs/SPLITBYCHROMOSOMEpbs



            ###############################
            echo -e "\n################# constructing qsub for realrecal\n"
            # appending the generic header to the qsub
            cat $outputdir/qsubGenericHeader > $qsub_realrecal
            echo "#PBS -N ${pipeid}_realrecal.${SampleName}.${chr}" >> $qsub_realrecal
            echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal
            echo "#PBS -o $RealignOutputDir/${SampleName}/logs/log.realrecal.${SampleName}.${chr}.ou" >> $qsub_realrecal
            echo "#PBS -e $RealignOutputDir/${SampleName}/logs/log.realrecal.${SampleName}.${chr}.in" >> $qsub_realrecal
            echo "#PBS -W depend=afterok:$split_job" >> $qsub_realrecal
            echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_realrecal

            echo "aprun -n 1 -N 1 -d $thr /bin/bash $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr}" >> $qsub_realrecal

            realrecal_job=`qsub $qsub_realrecal`
            `qhold -h u $realrecal_job`
            #`qrls -h u $split_job` 
            echo $realrecal_job >> $RealignOutputLogs/REALRECALpbs



            if [ $skipvcall == "NO" ]
            then
               qsub_vcall=$VcallOutputDir/${SampleName}/logs/qsub.vcallgatk.${SampleName}.${chr}



               ###############################
               echo -e "\n################# constructing qsub for vcall\n"
               # appending the generic header to the qsub
               cat $outputdir/qsubGenericHeader > $qsub_vcall
               echo "#PBS -N ${pipeid}_vcall_${SampleName}.${chr}" >> $qsub_vcall
               echo "#PBS -l walltime=$pbscpu" >> $qsub_vcall
               echo "#PBS -o $VcallOutputDir/${SampleName}/logs/log.vcall.${SampleName}.${chr}.ou" >> $qsub_vcall
               echo "#PBS -e $VcallOutputDir/${SampleName}/logs/log.vcall.${SampleName}.${chr}.in" >> $qsub_vcall
               echo "#PBS -W depend=afterok:$realrecal_job" >> $qsub_vcall
               echo -e "#PBS -l nodes=1:ppn=$thr\n" >> $qsub_vcall
   
               echo "aprun -n 1 -d $thr /bin/bash $VcallOutputDir/${SampleName}/logs/vcallgatk.${SampleName}.${chr}" >> $qsub_vcall
   
               vcall_job=`qsub $qsub_vcall`
               #`qrls -h u $realrecal_job`
               echo $vcall_job >> $VcallOutputLogs/VCALGATKpbs
            #else
               #`qrls -h u $realrecal_job` 
            fi
         done # going through chromosomes for this sample
      done <  $outputdir/SAMPLENAMES.list
      # end loop over samples
   ;;
   "QSUB")
      # will add later
   ;;
   "LAUNCHER")
      for chr in $indices
      do
         # clear out the joblists
         truncate -s 0 $RealignOutputLogs/split_bam_by_chromosome.${chr}.AnisimovJoblist
         truncate -s 0 $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist
         if [ $skipvcall == "NO" ]
         then
            truncate -s 0 $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist
         fi
         while read SampleName
         do
            # creating a qsub out of the job file
            # need to prepend "nohup" and append log file name, so that logs are properly created when Anisimov launches these jobs 
            split_bam_by_chromosome_log=${RealignOutputDir}/${SampleName}/logs/log.split.$SampleName.$chr.in
            awk -v awkvar_split_bam_by_chromosomelog=$split_bam_by_chromosome_log '{print "nohup "$0" > "awkvar_split_bam_by_chromosomelog}' $RealignOutputDir/${SampleName}/logs/split_bam_by_chromosome.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/jobfile.split_bam_by_chromosome.${SampleName}.${chr}
            echo "$RealignOutputDir/${SampleName}/logs/ jobfile.split_bam_by_chromosome.${SampleName}.${chr}" >> $RealignOutputLogs/split_bam_by_chromosome.${chr}.AnisimovJoblist

            realrecal_log=$RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.in
            awk -v awkvar_realrecallog=$realrecal_log '{print "nohup "$0" > "awkvar_realrecallog}' $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/jobfile.realrecal.${SampleName}.${chr}
            echo "$RealignOutputDir/${SampleName}/logs/ jobfile.realrecal.${SampleName}.${chr}" >> $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist

            if [ $skipvcall == "NO" ]
            then
               vcall_log=$VcallOutputDir/${SampleName}/logs/log.vcallgatk.${SampleName}.${chr}.in
               awk -v awkvar_vcalllog=$vcall_log '{print "nohup "$0" > "awkvar_vcalllog}' $VcallOutputDir/${SampleName}/logs/vcallgatk.${SampleName}.${chr} > $VcallOutputDir/${SampleName}/logs/jobfile.vcallgatk.${SampleName}.${chr}
               echo "$VcallOutputDir/${SampleName}/logs/ jobfile.vcallgatk.${SampleName}.${chr}" >> $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist
            fi
         done <  $outputdir/SAMPLENAMES.list
         # end loop over samples

         echo -e "\nscheduling Anisimov Launcher joblists\n"
         qsub_split_bam_by_chromosome_anisimov=$RealignOutputLogs/qsub.split_bam_by_chromosome.${chr}.AnisimovLauncher
         qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher
         if [ $skipvcall == "NO" ]
         then
            qsub_vcallgatk_anisimov=$VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher
         fi

         # appending the generic header to the qsub
         cat $outputdir/qsubGenericHeader > $qsub_split_bam_by_chromosome_anisimov
         cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov
         cat $outputdir/qsubGenericHeader > $qsub_vcallgatk_anisimov



         ###############
         ############### constructing the qsub for split_bam_by_chromosome
         echo "#PBS -N ${pipeid}_split_bam_by_chromosome_${chr}" >> $qsub_split_bam_by_chromosome_anisimov
         echo "#PBS -l walltime=$pbscpu" >> $qsub_split_bam_by_chromosome_anisimov
         echo "#PBS -o $RealignOutputLogs/log.split_bam_by_chromosome.${chr}.ou" >> $qsub_split_bam_by_chromosome_anisimov
         echo -e "#PBS -e $RealignOutputLogs/log.split_bam_by_chromosome.${chr}.in\n" >> $qsub_split_bam_by_chromosome_anisimov

         # split_bam_by_chromosome only contains single-threaded commands
         # -n should be equal to the number of Anisimov jobs + 1
         # -d 2 is there to space the jobs every other integer core: bunches of fp operations are being performed
         NumberOfProcesses=$(( samplecounter )) # already more than actual number of samples by 1
         NumberOfProcPerNode=16
         if [ $NumberOfProcesses -lt 17 ]
         then
            NumberOfNodes=1
            NumberOfProcPerNode=$(( samplecounter ))
         else
            NumberOfNodes=$(( NumberOfProcesses/NumberOfProcPerNode ))
            if [ `expr $NumberOfProcesses % $NumberOfProcPerNode` -gt 0 ]
            then
               (( NumberOfNodes++ )) # there is a remainder in that division, and we give those remaining jobs an extra node
            fi
         fi
         echo -e "#PBS -l nodes=$NumberOfNodes:ppn=$thr\n" >> $qsub_split_bam_by_chromosome_anisimov
         echo "aprun -n $NumberOfProcesses -N $NumberOfProcPerNode -d 2 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/split_bam_by_chromosome.${chr}.AnisimovJoblist /bin/bash > $RealignOutputLogs/split_bam_by_chromosome.${chr}.AnisimovLauncher.log" >> $qsub_split_bam_by_chromosome_anisimov
         # send an email about failures if any; must check for failures first         
         # echo "cat $RealignOutputLogs/FAILEDmessages | mail -s '[Task #3820]' "$redmine,$email"" >> $qsub_split_bam_by_chromosome_anisimov



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
      echo -e "\n going through chromosomes again to schedule all split before all realrecal, and all realrecal before vcallgatk, in order to efficiently work with a 25 job limit on queued state\n"
      # reset the list of SPLITBYCHROMOSOME pbs ids
      truncate -s 0 $RealignOutputLogs/SPLITBYCHROMOSOMEpbs
      truncate -s 0 $RealignOutputLogs/REALRECALpbs
      if [ $skipvcall == "NO" ]
      then
         truncate -s 0 $VcallOutputLogs/VCALGATKpbs
      fi
      cd $RealignOutputLogs # so that whatever temp fioles and pbs notifications would go there
      # first split_bam_by_chromosome
      for chr in $indices
      do
         split_bam_by_chromosome_job=`qsub $RealignOutputLogs/qsub.split_bam_by_chromosome.${chr}.AnisimovLauncher`
         `qhold -h u $split_bam_by_chromosome_job` #I am not going to allow these to run right away,
         echo $split_bam_by_chromosome_job >> $RealignOutputLogs/SPLITBYCHROMOSOMEpbs # will read these in to release them later
         # add dependency on split_bam_by_chromosome job to realrecal job
         sed -i "2i #PBS -W depend=afterok:$split_bam_by_chromosome_job" $RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher
      done
      # now realrecal
      for chr in $indices
      do
         realrecal_job=`qsub $RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher`
         # I am not going to put these on hold, as they will be helkd back by the dependency on respective split jobs
         # Add dependency on realrecal job to vcallgatk job
         sed -i "2i #PBS -W depend=afterok:$realrecal_job" $VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher
         echo $realrecal_job >> $RealignOutputLogs/REALRECALpbs
      done
      if [ $skipvcall == "NO" ]
      then
         # finally submit the vcallgatk
         cd $VcallOutputLogs # so that whatever temp fioles and pbs notifications would go there
         for chr in $indices
         do
            vcallgatk_job=`qsub $VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher`
            echo $vcallgatk_job >> $VcallOutputLogs/VCALGATKpbs
         done 
      fi

      # now release the split_bam_by_chromosome jobs
      split_bam_by_chromosome_ids=$( cat $RealignOutputLogs/SPLITBYCHROMOSOMEpbs | sed "s/\..*//" | tr "\n" " " )
      #qrls -h u $split_bam_by_chromosome_ids
   ;;
   "SERVER")
      while read SampleName
      do
         nohup $RealignOutputDir/${SampleName}/logs/jobfile.split_bam_by_chromosome.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/log.split.$SampleName.$chr.in
         nohup $RealignOutputDir/${SampleName}/logs/jobfile.realrecal.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.in 

         if [ $skipvcall == "NO" ]
         then
            nohup $VcallOutputDir/${SampleName}/logs/jobfile.vcallgatk.${SampleName}.${chr} > $VcallOutputDir/${SampleName}/logs/log.vcallgatk.${SampleName}.${chr}.in
         fi
      done  <  $outputdir/SAMPLENAMES.list
      # end loop over samples
   ;;
   esac


   echo "wrap up and produce summary table"

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


   # release all jobs now
   splitids=$( cat $RealignOutputLogs/SPLITBYCHROMOSOMEpbs | sed "s/\..*//" | tr "\n" " " )
   `qrls -h u $splitids`
   realrecalids=$( cat $RealignOutputLogs/REALRECALpbs | sed "s/\..*//" | tr "\n" " " )
   `qrls -h u $realrecalids`

   if [ $skipvcall == "NO" ]
   then
      vcalids=$( cat $VcallOutputLogs/VCALGATKpbs | sed "s/\..*//" | tr "\n" " " )
      `qrls -h u $vcalids`
   fi


fi
