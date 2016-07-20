#!/bin/bash
#
# start_realrecal_block.sh
# Second module in the GGPS analysis pipeline
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
   MSG="parameter mismatch. "
   echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
   exit 1;
fi
   echo -e "\n\n############# START REALRECAL BLOCK: decide if the aligned bams were made locally or externally; create the downstream schedule_realrecal_per_chromosome accordingly  ###############\n\n" >&2
   umask 0027
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


   sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
   outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
   thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
   refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
   scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
   refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
   input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
   paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
   rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
   multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   region=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
   resortbam=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2  )
   indices=$( echo $region | sed 's/:/ /g' )
   picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
   samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
   run_method=$( cat $runfile | grep -w RUNMETHOD | cut -d '=' -f2 )

   if [ ! -d $outputdir ]
   then
       mkdir -p $outputdir
   fi 
   
   if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
   then
      pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
      pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
   elif [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
   then
      pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
      pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
   else
      MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi

   if [ $resortbam != "1" -a $resortbam != "0" -a $resortbam != "YES" -a $resortbam != "NO" ]
   then
      MSG="Invalid value for RESORTBAM=$resortbam"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   else
      if [ $resortbam == "1" ]
      then
         $resortbam="YES"
      fi
      if [ $resortbam == "0" ]
      then
         $resortbam="NO"
      fi
   fi

   if [ ! -d $scriptdir ]
   then
      MSG="$scriptdir script directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi
   if [ ! -d $refdir ]
   then
      MSG="$refdir directory of reference genome  not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi
   if [ ! -d $picardir ]
   then
      MSG="$picardir picard directory  not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi


   # check that sampledir exists
   #if [ ! -d $sampledir ]
   #then
   #   MSG="INPUTDIR=$sampledir input directory not found"
   #   echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   #   exit 1;
   #fi

   numsamples=`wc -l $outputdir/SAMPLENAMES.list | cut -d ' ' -f 1`
   if [ $numsamples -lt 1 ]
   then
      MSG="No samples found in INPUTDIR=$sampledir."
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
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



   set +x; echo -e "\n\n" >&2;
   echo "####################################################################################################" >&2
   echo "#####################################                       ########################################" >&2
   echo "#####################################  CREATE  DIRECTORIES  ########################################" >&2
   echo "#####################################                       ########################################" >&2
   echo "####################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;



   #RealignOutputDir=$outputdir/realign
   #if [ -d $RealignOutputDir ]
   #then
   #   echo "$RealignOutputDir already exists; resetting it"
   #   `rm -r $RealignOutputDir/*`
   #else
   #   mkdir -p $RealignOutputDir
   #fi


   TopOutputLogs=$outputdir/logs
   if [ -d $TopOutputLogs ]
   then
      echo "$TopOutputLogs already exists"
      pbsids=""
   else
      mkdir -p $TopOutputLogs
   fi
   pipeid=$( cat $TopOutputLogs/pbs.CONFIGURE )



   RealignOutputLogs=$TopOutputLogs/realign
   if [ ! -d $RealignOutputLogs ]
   then
      mkdir $RealignOutputLogs
   fi
   #`chmod -R 770 $RealignOutputLogs/`
   # where messages about failures will go
   truncate -s 0 $RealignOutputLogs/FAILEDmessages
   if [ ! -d $RealignOutputLogs/FAILEDjobs ]
   then
      mkdir $RealignOutputLogs/FAILEDjobs
   else
      rm -r $RealignOutputLogs/FAILEDjobs/*
   fi
   #`chmod -R 770 $RealignOutputLogs/FAILEDjobs`




   set +x; echo -e "\n\n" >&2;
   echo "#############################################################################################################" >&2
   echo "#####################################                                ########################################" >&2
   echo "#####################################  THE FOLLWOING BLOCK IS FOR    ########################################" >&2
   echo "#####################################  BAM files aligned elsewhere   ########################################" >&2
   echo "#####################################                                ########################################" >&2
   echo "#############################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;

   listfiles="";
   sep=":";
   JOBSmayo=""  

   if [ $inputformat == "BAM" ]
   then
      set +x; echo -e "\n # alignment was done elsewhere; starting the workflow with realignment.\n" >&2; set -x
      
      while read SampleName
      do
         set +x; echo -e "\n # processing next sample\n" >&2; set -x
         if [ `expr ${#SampleName}` -lt 7 ]
         then
            set +x; echo -e "\n # skipping empty line\n" >&2; set -x
         else
            set +x; echo -e "\n # processing $SampleName. The line should have two fields <samplename> <bamfile>\n" >&2; set -x

            prefix=$( echo $sampledetail | cut -d ' ' -f1 )
            inbamfile=$( echo $sampledetail | cut -d ' ' -f2 )

            if [ ! -s $inbamfile ]
            then
               MSG="parsing $outputdir/SAMPLENAMES_multiplexed.list file failed. realignment failed to start"
               echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit 1;
            fi
            
            outputalign=$outputdir/$prefix/align
            outputlogs=$TopOutputLogs/$prefix/logs
            tmpbamfile=$inbamfile
            sortedplain=${prefix}.wrg.sorted.bam
            sorted=${prefix}.wdups.sorted.bam
            sortednodups=${prefix}.nodups.sorted.bam

            if [ ! -d $outputalign ]
            then
               mkdir -p $outputalign
               if [ ! -d $outputlogs ]
               then
                  mkdir -p $outputlogs
               else
                  `rm -r $outputlogs/*`
               fi
            fi

            if [ $resortbam == "YES" ]
            then

                set +x; echo -e "\n # $inbamfile needs to be resorted\n" >&2; set -x
            
            	qsub_sortbammayo=$outputlogs/qsub.sortbammayo.$prefix

            	echo "#PBS -A $pbsprj" >> $qsub_sortbammayo
            	echo "#PBS -N ${pipeid}_sortbamayo_${prefix}" >> $qsub_sortbammayo
            	echo "#PBS -l walltime=$pbscpu" >> $qsub_sortbammayo
            	echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_sortbammayo
            	echo "#PBS -o $outputlogs/log.sortbammayo.${prefix}.ou" >> $qsub_sortbammayo
            	echo "#PBS -e $outputlogs/log.sortbammayo.${prefix}.in" >> $qsub_sortbammayo
            	echo "#PBS -q $pbsqueue" >> $qsub_sortbammayo
            	echo "#PBS -m ae" >> $qsub_sortbammayo
            	echo "#PBS -M $email" >> $qsub_sortbammayo
            	echo "aprun -n 1 -d $thr $scriptdir/sortbammayo.sh $outputalign $tmpbamfile $sortedplain $sorted $sortednodups $runfile $outputlogs/log.sortbammayo.${prefix}.in $outputlogs/log.sortbammayo.${prefix}.ou $email $outputlogs/qsub.sortbammayo.${prefix}" >> $qsub_sortbammayo
            	#`chmod a+r $qsub_sortbammayo`
            	sortid=`qsub $qsub_sortbammayo`
            	#`qhold -h u $sortid`
            	echo $sortid >> $TopOutputLogs/pbs.REALSORTEDMAYO
            else
                set +x; echo -e "\n # $inbamfile DOES NOT need to be resorted. create symlinks\n" >&2; set -x
                ### TODO: header and index files need to be generated afresh
                
                cd $outputalign 
		ln -s $inbamfile $sortedplain
		ln -s $inbamfile $sorted
		ln -s $inbamfile $sortednodups            
            
            fi # end of resortbam if stmt
         fi # end of if statement checking for empty line in the SampleName file
      done <  $outputdir/SAMPLENAMES_multiplexed.list
      # end loop over input samples

      cp $TopOutputLogs/pbs.REALSORTEDMAYO $TopOutputLogs/ALN_MAYO_jobids
      JOBSmayo=$( cat $TopOutputLogs/ALN_MAYO_jobids | sed "s/\..*//g" | tr "\n" ":" | sed "s/::/:/g" )
   fi



   set +x; echo -e "\n\n" >&2;
   echo "#############################################################################################################" >&2
   echo "#####################################  END OF BLOCK FOR              ########################################" >&2
   echo "#####################################  BAM files aligned elsewhere   ########################################" >&2
   echo "#############################################################################################################" >&2
   echo -e "\n\n" >&2; set -x;

   # grab job ids for align and for preprocessing done in this module
   alignids=$( cat $TopOutputLogs/pbs.ALIGNED | sed "s/\..*//" | tr "\n" " " )
   mergeids=$( cat $TopOutputLogs/pbs.MERGED | sed "s/\..*//" | tr "\n" " " )
   sortedmayoids=$( cat $TopOutputLogs/pbs.REALSORTEDMAYO | sed "s/\..*//" | tr "\n" " " )


   set +x; echo -e "\n\n" >&2;
   echo "######################################################################################" >&2
   echo "#############   NOW THAT THE INPUT HAVE BEEN CHECKED AND RESORTED,  ##################" >&2
   echo "#############   WE CAN PROCEED TO SCHEDULE REAL/RECAL ETC           ##################" >&2
   echo "######################################################################################" >&2
   echo -e "\n\n" >&2; set -x;


   # the schedule_realrecal_per_chromosome should run on a MOM node, so submitting with qsub without aprun

   qsub_schedule_realrecal_per_chromosome=$RealignOutputLogs/qsub.schedule_realrecal_per_chromosome
   cat $outputdir/qsubGenericHeader > $qsub_schedule_realrecal_per_chromosome
   echo "#PBS -N ${pipeid}_schedule_realrecal_per_chromosome" >> $qsub_schedule_realrecal_per_chromosome
   echo "#PBS -l walltime=01:00:00" >> $qsub_schedule_realrecal_per_chromosome
   echo "#PBS -l nodes=1:ppn=1" >> $qsub_schedule_realrecal_per_chromosome
   echo "#PBS -o $RealignOutputLogs/log.schedule_realrecal_per_chromosome.ou" >> $qsub_schedule_realrecal_per_chromosome
   echo "#PBS -e $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in" >> $qsub_schedule_realrecal_per_chromosome

   # inserting dependencies 
   if [ `expr ${#JOBSmayo}` -gt 0 ]
   then
        echo "1i #PBS -W depend=afterok:$JOBSmayo" >> $qsub_schedule_realrecal_per_chromosome
   fi

   # choosing the script to run in the qsub script
   if [ $analysis != "MULTIPLEXED" ]
   then
        echo "$scriptdir/schedule_realrecal_per_chromosome.sh $outputdir $runfile $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in $RealignOutputLogs/log.schedule_realrecal_per_chromosome.ou $email $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome" >> $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome
   else
        echo "$scriptdir/schedule_realrecal_multiplexed.sh $outputdir $runfile $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in $RealignOutputLogs/log.schedule_realrecal_per_chromosome.ou $email $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome" >> $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome
   fi

   #`chmod a+r $qsub_schedule_realrecal_per_chromosome`               
   realrecaljob=`qsub $qsub_schedule_realrecal_per_chromosome`
   # `qhold -h u $realrecaljob` 
   echo $realrecaljob >> $TopOutputLogs/pbs.RECALL


   set +x; echo -e "\n ###################### now making PBS log files read accessible to the group #################################\n" >&2; set -x
   echo `date`
   #`chmod -R 770 $outputdir/`
   #`chmod -R 770 $TopOutputLogs/`


   find $outputdir -name logs -type d | awk '{print "chmod -R g=rwx "$1}' | sh -x

   echo `date`
