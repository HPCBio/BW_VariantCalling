#!/bin/bash
#
# start_realrecal_block.sh
# Second module in the GGPS analysis pipeline
##redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
   MSG="parameter mismatch. "
   echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
   exit 1;
fi
   echo -e "\n\n############# START REALRECAL BLOCK: decide if the aligned bams were made locally or externally; create the downstream schedule_realrecal_per_chromosome accordingly  ###############\n\n" >&2

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
   analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   realrecalflag=$( cat $runfile | grep -w REALIGNORDER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
   paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
   rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
   multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
   region=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
   resortflag=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
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

   if [ $resortflag != "1" -a $resortflag != "0" -a $resortflag != "YES" -a $resortflag != "NO" ]
   then
      MSG="Invalid value for RESORTBAM=$resortflag"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   else
      if [ $resortflag == "1" ]
      then
         $resortflag="YES"
      fi
      if [ $resortflag == "0" ]
      then
         $resortflag="NO"
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
   if [ ! -d $sampledir ]
   then
      MSG="INPUTDIR=$sampledir input directory not found"
      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
      exit 1;
   fi

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


   if [ $realrecalflag != "1" -a $realrecalflag != "0" ]
   then
      echo "realign-recalibration order flag is not set properly. Default value [1] will be assiged to it"
      realrecalflag="1"
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
   `chmod -R 770 $RealignOutputLogs/`
   # where messages about failures will go
   truncate -s 0 $RealignOutputLogs/FAILEDmessages
   if [ ! -d $RealignOutputLogs/FAILEDjobs ]
   then
      mkdir $RealignOutputLogs/FAILEDjobs
   else
      rm -r $RealignOutputLogs/FAILEDjobs/*
   fi
   `chmod -R 770 $RealignOutputLogs/FAILEDjobs`




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

   # finding all aligned BAMs to be realigned-recalibrated

   if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" -o $analysis == "MULTIPLEXED" ]
   then
      set +x; echo -e "\n # alignment was done inhouse. no need to resort\n" >&2; set -x
   else
      set +x; echo -e "\n # we need to check entries in samplefileinfo before launching other realignment analyses\n" >&2; set -x
      while read sampledetail
      do
         set +x; echo -e "\n # processing next line in file ...\n" >&2; set -x
         if [ `expr ${#sampledetail}` -lt 7 ]
         then
            set +x; echo -e "\n # skip empty line\n" >&2; set -x
         else
            set +x; echo -e "\n # preprocessing for realignment $sampledetail\n" >&2; set -x
            bamfile=$( echo $sampledetail )
            #bamfile=$( echo $sampledetail | cut -d '=' -f2 )
            #sampletag=$( echo $sampledetail | cut -d '=' -f1 | cut -d ':' -f2 )
            if [ `expr ${#PrefixName}` -lt 1 ]
            then
               MSG="parsing SAMPLEFILENAMES file failed. realignment failed to start."
               echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit 1;
            fi
         
            if [ ! -s $bamfile ]
            then
               MSG="parsing SAMPLEFILENAMES file failed. realignment failed to start"
               echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit 1;
            fi
         fi
      done <  $outputdir/SAMPLENAMES.list
      # end loop over input samples

   fi   

   if [ $resortflag == "YES" -a $analysis == "REALIGN_ONLY" ]
   then
      set +x; echo -e "\n # alignment was done separately; starting the workflow with realignment. Need to resort bam files. Checking input files\n" >&2; set -x
      # loop over samples, by reading the SampleName file we constructed above: $TopOutputLogs/SAMPLENAMES.list
      while read SampleName
      do
         set +x; echo -e "\n # processing next sample\n" >&2; set -x
         if [ `expr ${#SampleName}` -lt 7 ]
         then
            set +x; echo -e "\n # skipping empty line\n" >&2; set -x
         else
            set +x; echo -e "\n # realigning: $SampleName\n" >&2; set -x

            prefix=`basename $SampleName .wrg.sorted.bam`
            outputalign=$outputdir/align/$prefix
            outputlogs=$TopOutputLogs/align/$prefix
            tmpbamfile=$SampleName
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

            qsub_sortbammayo=$outputlogs/qsub.sortbammayo.$prefix
            echo "#PBS -V" > $qsub_sortbammayo
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
            `chmod a+r $qsub_sortbammayo`
            sortid=`qsub $qsub_sortbammayo`
            #`qhold -h u $sortid`
            echo $sortid >> $TopOutputLogs/pbs.REALSORTEDMAYO

            listfiles=$outputalign/$sorted${sep}${listfiles}
         fi # end of if statement checking for empty line in the SampleName file
      done <  $outputdir/SAMPLENAMES.list
      # end loop over input samples

      cp $TopOutputLogs/pbs.REALSORTEDMAYO $TopOutputLogs/ALN_MAYO_jobids
      JOBSmayo=$( cat $TopOutputLogs/ALN_MAYO_jobids | sed "s/\..*//g" | tr "\n" ":" | sed "s/::/:/g" )
   fi

   if [ $resortflag == "NO" -a $analysis == "REALIGN_ONLY" ]
   then
      set +x; echo -e "\n # starting the workflow with realignment, alignment step was done separartely. BAM files will not be resorted\n" >&2; set -x
      if [ $revertsam == "0" -o $revertsam == "NO" ]
      then
         set +x; echo -e "\n # nput is aligned bam that is suitable for realignment and recalibration... no need for preprocessing\n" >&2; set -x
         while read sampledetail
         do
            bam=$( echo $sampledetail )
            listfiles=${bam}${sep}${listfiles}
         done <  $outputdir/SAMPLENAMES.list
         # end loop over input samples
      else
         MSG="invalid value for preprocessing this kind of input: aligned bam. set RESORTBAM=YES and rerun the pipeline"
         echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
         exit 1;
      fi
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

   if [ $analysis != "MULTIPLEXED" ]
   then
       echo "$scriptdir/schedule_realrecal_per_chromosome.sh $outputdir $runfile $realrecalflag $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in $RealignOutputLogs/log.schedule_realrecal_per_chromosome.ou $email $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome" > $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome
   else
       echo "$scriptdir/schedule_realrecal_multiplexed.sh $outputdir $runfile $realrecalflag $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in $RealignOutputLogs/log.schedule_realrecal_per_chromosome.ou $email $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome" > $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome
   fi

   set +x; echo -e "\n # schedule the schedule_realrecal_per_chromosome job(s)\n" >&2; set -x
   case $run_method in
   "APRUN"|"QSUB"|"LAUNCHER")
      # schedule a single qsub for all samples
      qsub_schedule_realrecal_per_chromosome=$RealignOutputLogs/qsub.schedule_realrecal_per_chromosome

      ################## constructing the qsub by editing the file in place
      sed -i "1i \ " $qsub_schedule_realrecal_per_chromosome # this will separate the header from command

      # inserting dependencies; 
      # does not  depend on NCSA alignment jobs, as start_realrecal_block.sh depends on them
      if [ `expr ${#JOBSmayo}` -gt 0 ]
      then
         sed -i "1i #PBS -W depend=afterok:$JOBSmayo" $qsub_schedule_realrecal_per_chromosome
      fi

      # appending the PBS options to the file in reverse order
      sed -i "1i #PBS -l walltime=03:00:00" $qsub_schedule_realrecal_per_chromosome # allowing an hour for schedule_realrecal_per_chromosome: 
      # should be more than enough (it only takes ~5 minutes), and increases job priority
      sed -i "1i #PBS -l nodes=1:ppn=1" $qsub_schedule_realrecal_per_chromosome
      sed -i "1i #PBS -o $RealignOutputLogs/log.schedule_realrecal_per_chromosome.ou" $qsub_schedule_realrecal_per_chromosome
      sed -i "1i #PBS -e $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in" $qsub_schedule_realrecal_per_chromosome

      # appending the generic header to the qsub
      cat $outputdir/qsubGenericHeader $qsub_schedule_realrecal_per_chromosome > $qsub_schedule_realrecal_per_chromosome.tmp && mv $qsub_schedule_realrecal_per_chromosome.tmp $qsub_schedule_realrecal_per_chromosome

      sed -i "1i #PBS -N ${pipeid}_schedule_realrecal_per_chromosome" $qsub_schedule_realrecal_per_chromosome

      `chmod a+r $qsub_schedule_realrecal_per_chromosome`               
      realrecaljob=`qsub $qsub_schedule_realrecal_per_chromosome`
      # `qhold -h u $realrecaljob` 
      echo $realrecaljob >> $TopOutputLogs/pbs.RECALL
      ;;
   "SERVER")
      nohup $RealignOutputLogs/qsub.schedule_realrecal_per_chromosome > $RealignOutputLogs/log.schedule_realrecal_per_chromosome.in 
      ;;
   esac



   set +x; echo -e "\n ###################### done scheduling realign/recalibrate.#################################\n" >&2; set -x
   echo `date`
   `chmod -R 770 $outputdir/`
   `chmod -R 770 $TopOutputLogs/`

