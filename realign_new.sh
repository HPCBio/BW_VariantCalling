#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#
#  script to realign and recalibrate the aligned file(s)
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 7 ]
then
   MSG="parameter mismatch."
   echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
   #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
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
   fi
   if [ -s $refdir/$kgenome ]
   then
      realparms=$realparms":-known:$refdir/$kgenome"
      recalparms=$recalparms":--knownSites:$refdir/$kgenome"
   fi
   if [ ! -d $RealignOutputLogs ]
   then
      MSG="$RealignOutputLogs realignlog directory not found"
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi


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
   elif [ $schedule eq "LAUNCHER" ]
      then
      truncate -s 0 $RealignOutputLogs/$realrecal.AnisimovJoblist
   fi


   pipeid=$( cat $TopOutputLogs/MAINpbs )






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

      if [ $multisample == "NO" ]
      then

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
      fi
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
         # now check that there is only one file
         aligned_bam=( $aligned_bam ) # recast variable as an array and count the number of members
         if [ ${#aligned_bam[@]} -ne 1 ]
         then
            MSG="more than one bam file found in $outputdir/align/${SampleName}"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline
            exit 1;
         fi

         aligned_bam="${aligned_bam[*]}" # recast variable back as string
         aligned_bam=$outputdir/align/${SampleName}/${aligned_bam}
         sID=$SampleName
         sPU=$SampleName
         sSM=$SampleName
         RGparms=$( echo "RGID=${sID}:RGLB=${sLB}:RGPU=${sPU}:RGSM=${sSM}:RGPL=${sPL}:RGCN=${sCN}" )



         echo -e "\n######################################  assemble the sortnode sub-block #####################################\n"
         echo "$scriptdir/sortnode.sh $picardir $samdir $javamodule $RealignOutputDir/$SampleName $aligned_bam ${SampleName}.$chr.bam ${SampleName}.$chr.sorted.bam $RGparms $flag $chr $RealignOutputDir/${SampleName}/logs/log.sort.$SampleName.$chr.in $RealignOutputDir/${SampleName}/logs/log.sort.$SampleName.$chr.ou $email $RealignOutputDir/${SampleName}/logs/sortnode.${SampleName}.${chr}" > $RealignOutputDir/${SampleName}/logs/sortnode.${SampleName}.${chr}

         # construct the list of all input files for the GATK real-recal procedure, multisample case
         if [ $multisample == "YES" ]
         then
            chrinfiles[$chromosomecounter]=${chrinfiles[$chromosomecounter]}":-I:$RealignOutputDir/${SampleName}/${SampleName}.$chr.sorted.bam"
# not used anymore?          chrinputfiles[$inx]=${chrinputfiles[$inx]}":INPUT=$outputdir/$SampleName.$chr.sorted.bam"
         else
            echo -e "\n#######################   assemble the real-recall/var call sub-block for INDEPENDENT SAMPLES    ################################\n"
            echo "$scriptdir/realrecalold.sh $RealignOutputDir/${SampleName} $chr.realrecal.$SampleName.output.bam $chr -I:$RealignOutputDir/${SampleName}/${SampleName}.$chr.sorted.bam ${region[$chromosomecounter]} $realparms $recalparms $runfile $flag $RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.in $RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.ou $email $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr}" > $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr}
           
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
#         echo "###################    assemble the sortnode sub-block for MULTISAMPLE experiment    ############################"
#         echo "##" # this line makes output logs more readable
###############################################################################
###############################################################################
#########################   NOT CLEAR HOW TO HANDLE BAMFILE AND rg VARIABLES HERE IN THE MULTISAMPLE CASE
#######################3          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
#####################################################3333
         #echo "$scriptdir/realrecalold.sh $RealignOutputDir $chr.realrecal.output.bam $chr ${chrinfiles[$chromosomecounter]} ${region[$chromosomecounter]} $realparms $recalparms $runfile $flag $RealignOutputLogs/log.realrecal.$bamfile.$chr.in $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou $email $RealignOutputLogs/qsub.realrecal.$bamfile.$chr" >> $qsub_realrecalold
#      fi


      ((  chromosomecounter++ )) 
   done # done going through chromosomes 
   # at the end of this set of nested loops, the variables chromosomecounter and samplecounter
   # reflect, respectively, number_of_chromosomes+1 and number_of_samples+1,
   # which is exactly the number of nodes required for anisimov launcher 








#         qsub_sortnode=$RealignOutputLogs/qsub.sort.${SampleName}.$chr
#         echo "#PBS -N ${pipeid}_sort_${SampleName}_$chr" >> $qsub_sortnode
#         echo "#PBS -l walltime=$pbscpu" >> $qsub_sortnode
#         echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_sortnode
#         echo "#PBS -o $RealignOutputDir/${SampleName}/logs/log.sort.${SampleName}.$chr.ou" >> $qsub_sortnode
#         echo "#PBS -e $RealignOutputDir/${SampleName}/logs/log.sort.${SampleName}.$chr.in" >> $qsub_sortnode
# 
#         `chmod a+r $qsub_sortnode`
#            sortjobid=`qsub $qsub_sortnode`
#            # new line to avoid hiccup
#            #`qhold -h u $sortjobid`
#            echo $sortjobid >> $outputdir/logs/REALSORTEDpbs.$SampleName
#            echo $sortjobid >> $outputdir/logs/REALSORTED.$SampleName$chr
#
#            truncate -s 0 $realrecal.$SampleName.$chr.serialjobs
#            echo "$RealignOutputLogs realrecal.$SampleName.$chr.serialjobs" >> $RealignOutputLogs/$realrecal.AnisimovJoblist
#






#         if [ $multisample == "NO" ]
#         then


#      sortid=$( cat $outputdir/logs/REALSORTED.$bamfile_$chr | sed "s/\..*//g" | tr "\n" ":" )
#      outputfile=$chr.realrecal.$bamfile.output.bam
#      igv_files=${igv_files}":INPUT=${outputfile}"
#      echo "realign-recalibrate for interval:$chr..."
#      qsub_realrecalold=$RealignOutputLogs/qsub.realrecal.$bamfile.$chr
#      echo "#PBS -V" > $qsub_realrecalold
#      echo "#PBS -A $pbsprj" >> $qsub_realrecalold
#      echo "#PBS -N ${pipeid}_realrecal_$bamfile.$chr" >> $qsub_realrecalold
#      echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecalold
#      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_realrecalold
#      echo "#PBS -o $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou" >> $qsub_realrecalold
#      echo "#PBS -e $RealignOutputLogs/log.realrecal.$bamfile.$chr.in" >> $qsub_realrecalold
#      echo "#PBS -q $pbsqueue" >> $qsub_realrecalold
#      echo "#PBS -m ae" >> $qsub_realrecalold
#      echo "#PBS -M $email" >> $qsub_realrecalold
#      echo "#PBS -W depend=afterok:$sortid" >> $qsub_realrecalold
#
#      if [ $schedule eq "QSUB" ]
#      then
#         echo "aprun -n 1 -d $thr $scriptdir/realrecalold.sh $????#####outputdir $outputfile $chr ${chrinfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $RealignOutputLogs/log.realrecal.$bamfile.$chr.in $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou $email $RealignOutputLogs/qsub.realrecal.$bamfile.$chr" >> $qsub_realrecalold
#         `chmod a+r $qsub_realrecalold`
#         #recaljobid=`qsub $qsub_realrecalold`
#         # new line to avoid hiccup
#         #`qhold -h u $recaljobid`
#         echo $recaljobid >> $outputdir/logs/REALRECALpbs.$bamfile
#      else 
#         # this block constructs the list of jobs to perform in serial on that chromosome, on that sample
#         echo "nohup $scriptdir/realrecalold.sh $????#####outputdir $outputfile $chr ${chrinfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $RealignOutputLogs/log.realrecal.$bamfile.$chr.in $RealignOutputLogs/log.realrecal.$bamfile.$chr.ou $email $RealignOutputLogs/qsub.realrecal.$bamfile.$chr > $RealignOutputLogs/log.realrecal.$bamfile.$chr.in" >> $RealignOutputLogs/realrecal.$bamfile.$chr.serialjobs
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
#         echo "#PBS -m ae" >> $qsub_vcallgatk
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
      while read SampleName
      do
         # using sed to edit the first and only line in each jobfile to add the relevant scheduling commands
         sed "1!b;s/^/aprun -n 1 -d $thr /" $RealignOutputDir/${SampleName}/logs/sortnode.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/qsub.sortnode.${SampleName}.${chr}
         sed "1!b;s/^/aprun -n 1 -d $thr /" $RealignOutputDir/${SampleName}/logs/realrecal.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/qsub.realrecal.${SampleName}.${chr}
         if [ $skipvcall == "NO" ]
         then
            sed "1!b;s/^/aprun -n 1 -d $thr /" $VcallOutputDir/${SampleName}/logs/vcallgatk.${SampleName}.${chr} > $VcallOutputDir/${SampleName}/logs/qsub.vcallgatk.${SampleName}.${chr}
         fi
      done <  $outputdir/SAMPLENAMES.list
      # end loop over samples
   ;;
   "QSUB")

   ;;
   "LAUNCHER")
      for chr in $indices
      do
         # clear out the joblists
         truncate -s 0 $RealignOutputLogs/sortnode.${chr}.AnisimovJoblist
         truncate -s 0 $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist
         truncate -s 0 $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist
         while read SampleName
         do
            # creating a qsub out of the job file
            # need to prepend "nohup" and append lof file name, so that logs are properly created when Anisimov launches these jobs 
            sortnode_log=${RealignOutputDir}/${SampleName}/logs/log.sort.$SampleName.$chr.in
            awk -v awkvar_sortnodelog=$sortnode_log '{print "nohup "$0" > "awkvar_sortnodelog}' $RealignOutputDir/${SampleName}/logs/sortnode.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/jobfile.sortnode.${SampleName}.${chr}
            echo "$RealignOutputDir/${SampleName}/logs/ jobfile.sortnode.${SampleName}.${chr}" >> $RealignOutputLogs/sortnode.${chr}.AnisimovJoblist

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
         qsub_sortnode_anisimov=$RealignOutputLogs/qsub.sortnode.${chr}.AnisimovLauncher
         qsub_realrecal_anisimov=$RealignOutputLogs/qsub.realrecal.${chr}.AnisimovLauncher
         qsub_vcallgatk_anisimov=$VcallOutputLogs/qsub.vcalgatk.${chr}.AnisimovLauncher

         # appending the generic header to the qsub
         cat $outputdir/qsubGenericHeader > $qsub_sortnode_anisimov
         cat $outputdir/qsubGenericHeader > $qsub_realrecal_anisimov
         cat $outputdir/qsubGenericHeader > $qsub_vcallgatk_anisimov



         ###############
         ############### constructing the qsub for sortnode
         echo "#PBS -N ${pipeid}_sortnode_${chr}" >> $qsub_sortnode_anisimov
         echo "#PBS -l walltime=$pbscpu" >> $qsub_sortnode_anisimov
         echo "#PBS -o $RealignOutputLogs/log.sortnode.${chr}.ou" >> $qsub_sortnode_anisimov
         echo -e "#PBS -e $RealignOutputLogs/log.sortnode.${chr}.in\n" >> $qsub_sortnode_anisimov

         # sortnode only contains single-threaded commands
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
         echo -e "#PBS -l nodes=$NumberOfNodes:ppn=$thr\n" >> $qsub_sortnode_anisimov
         echo "aprun -n $NumberOfProcesses -N $NumberOfProcPerNode -d 2 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/sortnode.${chr}.AnisimovJoblist /bin/bash > $RealignOutputLogs/sortnode.${chr}.AnisimovLauncher.log" >> $qsub_sortnode_anisimov

         sortnode_job=`qsub $qsub_sortnode_anisimov`
         `qhold -h u $sortnode_job`



         ###############
         ############### constructing the qsub for realrecal
         echo "#PBS -N ${pipeid}_realrecal_${chr}" >> $qsub_realrecal_anisimov
         echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecal_anisimov
         echo "#PBS -o $RealignOutputLogs/log.realrecal.${chr}.ou" >> $qsub_realrecal_anisimov
         echo -e "#PBS -e $RealignOutputLogs/log.realrecal.${chr}.in\n" >> $qsub_realrecal_anisimov
         # add dependency to realrecal job
         echo -e "#PBS -W depend=afterok:$sortnode_job\n" >> $qsub_realrecal_anisimov
         # realrecal and vcall actually use multithreaded processes,
         # so we will give each chromosome its own node
         # +1 for the launcher
         echo -e "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_realrecal_anisimov
         # the actual command
         echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $RealignOutputLogs/realrecal.${chr}.AnisimovJoblist /bin/bash > $RealignOutputLogs/realrecal.${chr}.AnisimovLauncher.log" >> $qsub_realrecal_anisimov

         realrecal_job=`qsub $qsub_realrecal_anisimov`
         `qhold -h u $realrecal_job`
         `qrls -h u $sortnode_job`



         ###############
         ############### constructing the qsub for vcallgatk
         echo "#PBS -N ${pipeid}_vcallgatk_${chr}" >> $qsub_vcallgatk_anisimov
         echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk_anisimov
         echo "#PBS -o $RealignOutputLogs/log.vcallgatk.${chr}.ou" >> $qsub_vcallgatk_anisimov
         echo -e "#PBS -e $RealignOutputLogs/log.vcallgatk.${chr}.in\n" >> $qsub_vcallgatk_anisimov
         # add dependency to realrecal job
         echo -e "#PBS -W depend=afterok:$realrecal_job\n" >> $qsub_vcallgatk_anisimov
         # realrecal and vcall actually use multithreaded processes, 
         # so we will give each chromosome its own node
         # +1 for the launcher
         echo -e "#PBS -l nodes=$chromosomecounter:ppn=$thr\n" >> $qsub_vcallgatk_anisimov
         # the actual command
         echo "aprun -n $chromosomecounter -N 1 -d 32 ~anisimov/scheduler/scheduler.x $VcallOutputLogs/vcallgatk.${chr}.AnisimovJoblist /bin/bash > $VcallOutputLogs/vcallgatk.${chr}.AnisimovLauncher.log" >> $qsub_vcallgatk_anisimov

         vcallgatk_job=`qsub $qsub_vcallgatk_anisimov`
         `qrls -h u $realrecal_job`

      done # done going through chromosomes
   ;;
   "SERVER")
      while read SampleName
      do
         nohup $RealignOutputDir/${SampleName}/logs/jobfile.sortnode.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/log.sort.$SampleName.$chr.in
         nohup $RealignOutputDir/${SampleName}/logs/jobfile.realrecal.${SampleName}.${chr} > $RealignOutputDir/${SampleName}/logs/log.realrecal.$SampleName.$chr.in 

         if [ $skipvcall == "NO" ]
         then
            nohup $VcallOutputDir/${SampleName}/logs/jobfile.vcallgatk.${SampleName}.${chr} > $VcallOutputDir/${SampleName}/logs/log.vcallgatk.${SampleName}.${chr}.in
         fi
      done  <  $outputdir/SAMPLENAMES.list
      # end loop over samples
   ;;
   esac



#      echo "scheduling Anisimov Launcher joblist"
#      qsub_aisimov=$varlogdir/qsub.vcallgatk.$bamfile.$chr
#      echo "#PBS -V" > $qsub_aisimov
#      echo "#PBS -A $pbsprj" >> $qsub_aisimov
#      echo "#PBS -N ${pipeid}_vcall_$chr" >> $qsub_aisimov
#      echo "#PBS -l walltime=$pbscpu" >> $qsub_aisimov
#      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_aisimov
#      echo "#PBS -o $varlogdir/log.vcallgatk.$bamfile.$chr.ou" >> $qsub_aisimov
#      echo "#PBS -e $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $qsub_aisimov
#      echo "#PBS -q $pbsqueue" >> $qsub_aisimov
#      echo "#PBS -m ae" >> $qsub_aisimov
#      echo "#PBS -M $email" >> $qsub_aisimov
#
#      echo "aprun -n $numsamples -N 1 -d $thr ~anisimov/scheduler/scheduler.x $RealignOutputLogs/$realrecal.AnisimovJoblist /bin/bash > $RealignOutputLogs/RealRecalVCallAnisimov.joblist.log" >> $qsub_anisimov
#      $RealignOutputLogs/$realrecal.AnisimovJoblist
 
##     qsub $qsub_anisimov
#
#   fi
#
#   echo "clean up and produce summary table"
#   if [ $skipvcall == "NO" ]
#   then
#      listjobids=$( cat $outputdir/logs/VCALLGATKpbs.$bamfile | sed "s/\..*//g" | tr "\n" ":" )
#   else
#      listjobids=$( cat $outputdir/logs/REALRECALpbs.$bamfile | sed "s/\..*//g" | tr "\n" ":" )
#   fi
#
#   # preparing output files for realignment and/or variant calling
#
#   qsub_igvbam=$outputdir/logs/qsub.igvbam.$bamfile
#   echo "#PBS -V" > $qsub_igvbam
#   echo "#PBS -A $pbsprj" >> $qsub_igvbam
#   echo "#PBS -N ${pipeid}_igvbam.$bamfile" >> $qsub_igvbam
#   echo "#PBS -l walltime=$pbscpu" >> $qsub_igvbam
#   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_igvbam
#   echo "#PBS -o $outputdir/logs/log.igvbam.$bamfile.ou" >> $qsub_igvbam
#   echo "#PBS -e $outputdir/logs/log.igvbam.$bamfile.in" >> $qsub_igvbam
#   echo "#PBS -q $pbsqueue" >> $qsub_igvbam
#   echo "#PBS -m ae" >> $qsub_igvbam
#   echo "#PBS -M $email" >> $qsub_igvbam
#   echo "#PBS -W depend=afterok:$listjobids" >> $qsub_igvbam
#   echo "aprun -n 1 -d $thr $scriptdir/igvbam.sh $????#####outputdir $igv_files $runfile $outputdir/logs/log.igvbam.$bamfile.in $outputroot/logs/log.igvbam.$bamfile.ou $email $outputdir/logs/qsub.igvbam.$bamfile"  >> $qsub_igvbam
#   `chmod a+r $qsub_igvbam`
##   igvjobid=`qsub $qsub_igvbam`
#   echo $igvjobid >> $outputdir/logs/IGVBAMpbs.$bamfile
#
#   if [ $cleanupflag == "YES" ]
#   then
#      qsub_cleanup=$outputdir/logs/qsub.cleanup.realn.$bamfile
#      echo "#PBS -V" > $qsub_cleanup
#      echo "#PBS -A $pbsprj" >> $qsub_cleanup
#      echo "#PBS -N ${pipeid}_cleanup_realn.$bamfile" >> $qsub_cleanup
#      echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
#      echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
#      echo "#PBS -o $outputdir/logs/log.cleanup.realn.$bamfile.ou" >> $qsub_cleanup
#      echo "#PBS -e $outputdir/logs/log.cleanup.realn.$bamfile.in" >> $qsub_cleanup
#      echo "#PBS -q $pbsqueue" >> $qsub_cleanup
#      echo "#PBS -m ae" >> $qsub_cleanup
#      echo "#PBS -M $email" >> $qsub_cleanup
#      echo "#PBS -W depend=afterok:$igvjobid" >> $qsub_cleanup
#      echo "$scriptdir/cleanup.sh $outputdir $analysis $outputdir/logs/log.cleanup.realn.$bamfile.in $outputdir/logs/log.cleanup.realn.$bamfile.ou $email $outputdir/logs/qsub.cleanup.realn.$bamfile"  >> $qsub_cleanup
#      `chmod a+r $qsub_cleanup`
##      cleanjobid=`qsub $qsub_cleanup`
#      echo $cleanjobid >> $outputdir/logs/CLEANUPpbs.$bamfile
#   fi
#
#   if [ $skipvcall == "NO" ]
#   then
#      qsub_mergevcf=$varlogdir/qsub.mergevcf.$bamfile
#      echo "#PBS -V" > $qsub_mergevcf
#      echo "#PBS -A $pbsprj" >> $qsub_mergevcf
#      echo "#PBS -N ${pipeid}_mergevcf.$bamfile" >> $qsub_mergevcf
#      echo "#PBS -l walltime=$pbscpu" >> $qsub_mergevcf
#      echo "#PBS -l nodes=1:ppn=1" >> $qsub_mergevcf
#      echo "#PBS -o $varlogdir/log.mergevcf.$bamfile.ou" >> $qsub_mergevcf
#      echo "#PBS -e $varlogdir/log.mergevcf.$bamfile.in" >> $qsub_mergevcf
#      echo "#PBS -q $pbsqueue" >> $qsub_mergevcf
#      echo "#PBS -m ae" >> $qsub_mergevcf
#      echo "#PBS -M $email" >> $qsub_mergevcf
#      echo "#PBS -W depend=afterok:$listjobids" >> $qsub_mergevcf
#      echo "$scriptdir/mergevcf.sh $vardir $vcf_files $tabixdir $vcftoolsdir $varlogdir/log.mergevcf.$bamfile.in $varlogdir/log.mergevcf.$bamfile.ou $email $varlogdir/qsub.mergevcf.$bamfile"  >> $qsub_mergevcf
#      `chmod a+r $qsub_mergevcf`
##      mergevcfjobid=`qsub $qsub_mergevcf`
#      echo $mergevcfjobid >> $outputdir/logs/MERGEVCFpbs
#
#      if [ $cleanupflag == "YES" ]
#      then
#         qsub_cleanup=$outputdir/logs/qsub.cleanup.vcall.$bamfile
#         echo "#PBS -V" > $qsub_cleanup
#         echo "#PBS -A $pbsprj" >> $qsub_cleanup
#         echo "#PBS -N ${pipeid}_cleanup_vcall.$bamfile" >> $qsub_cleanup
#         echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
#         echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
#         echo "#PBS -o $outputdir/logs/log.cleanup.vcall.$bamfile.ou" >> $qsub_cleanup
#         echo "#PBS -e $outputdir/logs/log.cleanup.vcall.$bamfile.in" >> $qsub_cleanup
#         echo "#PBS -q $pbsqueue" >> $qsub_cleanup
#         echo "#PBS -m ae" >> $qsub_cleanup
#         echo "#PBS -M $email" >> $qsub_cleanup
#         echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub_cleanup
#         echo "$scriptdir/cleanup.sh $outputdir VCALL_ONLY $outputdir/logs/log.cleanup.vcall.$bamfile.in $outputdir/logs/log.cleanup.vcall.$bamfile.ou $email $outputdir/logs/qsub.cleanup.vcall.$bamfile"  >> $qsub_cleanup
#         `chmod a+r $qsub_cleanup`
##         cleanjobid=`qsub $qsub_cleanup`
#         echo $cleanjobid >> $outputdir/logs/CLEANUPpbs.$bamfile
#      fi
#   fi
#
#  # new lines to avoid hiccups and missing jobs in summary
#
#  heldjobs_realsorted=$( cat $outputdir/logs/REALSORTEDpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
#  heldjobs_realrecal=$( cat $outputdir/logs/REALRECALpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
#  heldjobs_vcall=$( cat $outputdir/logs/VCALLGATKpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
#  cleanupjobs=$( cat $outputdir/logs/CLEANUPpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
#
#  #`qrls -h u $heldjobs_realsorted`
#  #`qrls -h u $heldjobs_realrecal`
#  #`qrls -h u $heldjobs_vcall`
#  `sleep 10s`
#
#  qsub_summary=$outputdir/logs/qsub.summary.realn.allok.$bamfile
#  echo "#PBS -V" > $qsub_summary
#  echo "#PBS -A $pbsprj" >> $qsub_summary
#  echo "#PBS -N ${pipeid}_summaryok.$bamfile" >> $qsub_summary
#  echo "#PBS -l walltime=$pbscpu" >> $qsub_summary
#  echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
#  echo "#PBS -o $outputdir/logs/log.summary.realn.$bamfile.ou" >> $qsub_summary
#  echo "#PBS -e $outputdir/logs/log.summary.realn.$bamfile.in" >> $qsub_summary
#  echo "#PBS -q $pbsqueue" >> $qsub_summary
#  echo "#PBS -m ae" >> $qsub_summary
#  echo "#PBS -M $email" >> $qsub_summary
#  if [ `expr ${#cleanupjobs}` -gt 0 ]
#  then 
#     echo "#PBS -W depend=afterok:$cleanupjobs" >> $qsub_summary
#  else
#     if [ $skipvcall == "YES" ]
#     then
#        echo "#PBS -W depend=afterok:$igvjobid" >> $qsub_summary
#     else
#        echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub_summary
#     fi
#  fi
#  echo "$scriptdir/summary.sh $outputdir $email exitok"  >> $qsub_summary
#  `chmod a+r $qsub_summary`
##  lastjobid=`qsub $qsub_summary`
#  echo $lastjobid >> $outputdir/logs/SUMMARYpbs.$bamfile
#
#  if [ `expr ${#lastjobid}` -lt 1 ]
#  then
#     echo "at least one job aborted"
#     qsub_summary=$outputdir/logs/qsub.summary.realn.afterany.$bamfile
#     echo "#PBS -V" > $qsub_summary
#     echo "#PBS -A $pbsprj" >> $qsub_summary
#     echo "#PBS -N ${pipeid}_summary_afterany.$bamfile" >> $qsub_summary
#     echo "#PBS -l walltime=$pbscpu" >> $qsub_summary
#     echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
#     echo "#PBS -o $outputdir/logs/log.summary.realn.afterany.$bamfile.ou" >> $qsub_summary
#     echo "#PBS -e $outputdir/logs/log.summary.realn.afterany.$bamfile.in" >> $qsub_summary
#     echo "#PBS -q $pbsqueue" >> $qsub_summary
#     echo "#PBS -m ae" >> $qsub_summary
#     echo "#PBS -M $email" >> $qsub_summary
#     echo "#PBS -W depend=afterany:$listjobids" >> $qsub_summary
#     echo "$scriptdir/summary.sh $outputtopdir $email exitnotok"  >> $qsub_summary
#     `chmod a+r $qsub_summary`
##     badjobid=`qsub $qsub_summary`
#     echo $badjobid >> $outputdir/logs/SUMMARYpbs.$bamfile
#  fi
#  `chmod -R 770 $outputroordir/logs`
#  echo `date`
#fi



# howe to implement the case switrching for schweuler


#            if [ $run_method == "LAUNCHER" ]
#            then
#               echo "$scriptdir/realign_new.sh $RealignOutputDir/$SampleName $RealignOutputLogs $listfiles $runfile $realrecalflag $RealignOutputLogs/log.${SampleName}.realign_new.in $RealignOutputLogs/log.${SampleName}.realign_new.ou $email $RealignOutputLogs/qsub.${SampleName}.realign_new" > $RealignOutputLogs/$realignnew_${SampleName}_LauncherJob
#               echo "$RealignOutputLogs $realignnew_${SampleName}_LauncherJob" >> $RealignOutputLogs/$realignnew_LauncherJobList
#            else
#               # in this case does not matter whether with or without aprun: the realign_new should run on a MOM node anyway
#               echo "$scriptdir/realign_new.sh $RealignOutputDir/$SampleName $RealignOutputLogs $listfiles $runfile $realrecalflag $RealignOutputLogs/log.${SampleName}.realign_new.in $RealignOutputLogs/log.${SampleName}.realign_new.ou $email $RealignOutputLogs/qsub.${SampleName}.realign_new" >> $RealignOutputLogs/qsub.${SampleName}.realign_new
#            fi
#            #`chmod a+r $qsub_realignnew`
#            # realrecaljob=`qsub $qsub_realignnew`
#            # `qhold -h u $realrecaljob`
#            # echo $realrecaljob >> $TopOutputLogs/RECALLpbs
#         fi # finish checking for empty lines in SampleNames file
#
#         (( samplecounter++ ))
#      done <  $outputdir/SAMPLENAMES.list
#      # end loop over input samples



      ### schedule the realign_new job(s)
#      case $run_method in
#         APRUN|QSUB)
#            # schedule a separate qsub for realign_new for each sample
#            while read SampleName
#            do
#
#               qsub_realignnew=$RealignOutputLogs/qsub.${SampleName}.realign_new
##
#               # constructing the qsub
#               # appending the generic header to the
#               cat $????#####outputdir/qsubGenericHeader $qsub_realignnew > /tmp/qsub.${SampleName}.realign_new && mv /tmp/qsub.${SampleName}.realign_new $qsub_realignnew
#
#                sed -i "#PBS -N ${pipeid}_realign_new.${SampleName}" >> $qsub_realignnew
#               echo "#PBS -l walltime=01:00:00" >> $qsub_realignnew # allowing an hour for realign_new:
#               # should be more than enough (it only takes ~5 minutes), and increases job priority
#               echo "#PBS -l nodes=1:ppn=1" >> $qsub_realignnew
#               echo "#PBS -o $RealignOutputLogs/log.${SampleName}.realign_new.ou" >> $qsub_realignnew
#               echo "#PBS -e $RealignOutputLogs/log.${SampleName}.realign_new.in" >> $qsub_realignnew
#               if [ `expr ${#JOBSmayo}` -gt 0 ]
#               then
#                  echo "#PBS -W depend=afterok:$JOBSmayo" >> $qsub_realignnew
#   #           else
#   #              echo "#PBS -W depend=afterok:$JOBSncsa" >> $qsub_realignnew
#               fi
#               echo -e "\n" >> $qsub_realignnew
#            done <  $outputdir/SAMPLENAMES.list
#            # end loop over input samples
#
#         ;;
#         LAUNCHER)
#            commands
#         ;;
#         SERVER)
#            commands
#         ;;
#      esac



#aprun -n 1 -d $thr sortnode
#
#         qsub_sortnode=$RealignOutputLogs/qsub.sort.${SampleName}.$chr
#         echo "#PBS -N ${pipeid}_sort_${SampleName}_$chr" >> $qsub_sortnode
#         echo "#PBS -l walltime=$pbscpu" >> $qsub_sortnode
#         echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_sortnode
#         echo "#PBS -o $RealignOutputDir/${SampleName}/logs/log.sort.${SampleName}.$chr.ou" >> $qsub_sortnode
#         echo "#PBS -e $RealignOutputDir/${SampleName}/logs/log.sort.${SampleName}.$chr.in" >> $qsub_sortnode
#
#         `chmod a+r $qsub_sortnode`
#            sortjobid=`qsub $qsub_sortnode`
#            # new line to avoid hiccup
#            #`qhold -h u $sortjobid`
#            echo $sortjobid >> $outputdir/logs/REALSORTEDpbs.$SampleName
#            echo $sortjobid >> $outputdir/logs/REALSORTED.$SampleName$chr
#
#            truncate -s 0 $realrecal.$SampleName.$chr.serialjobs
#            echo "$RealignOutputLogs realrecal.$SampleName.$chr.serialjobs" >> $RealignOutputLogs/$realrecal.AnisimovJoblist

fi
