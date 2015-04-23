#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=lmainzer@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 6 ]
then
        MSG="Parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO. Reason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO. Reason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        runmode=$2
        elog=$3
        olog=$4
        email=$5
        qsubfile=$6
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        resortbam=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 )
        bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        sampledir=$( cat $runfile | grep -w SAMPLEDIR | cut -d '=' -f2 )

        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi

        if [ -z $thr -o -z $outputdir -o -z $pbsprj -o -z $epilogue ]
        then
 		MSG="Invalid value specified for any of these paramaters in configuration file:\nPBSTHREADS=$thr\nOUTPUTDIR=$outputdir\nPBSPROJECTID=$pbsprj\nEPILOGUE=$epilogue"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
        fi

        if [ $resortbam != "1" -a $resortbam != "0" -a $resortbam != "YES" -a $resortbam != "NO" ]
        then
           MSG="Invalid value for RESORTBAM=$resortbam"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
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

        if [ $bamtofastqflag != "YES" -a $bamtofastqflag != "NO" -a $bamtofastqflag != "1" -a $bamtofastqflag != "0" ]
        then
            MSG="BM2FASTQFLAG=$bamtofastqflag  invalid value"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        else
            if [ $bamtofastqflag == "1" ]
            then
                $bamtofastqflag="YES"
            fi
            if [ $bamtofastqflag == "0" ]
            then
                $bamtofastqflag="NO"
            fi
        fi


        if [ $resortbam == "YES" -a $bamtofastqflag == "YES" ]
        then
            MSG="Incompatible values for the pair RESORTBAM=$resortbam and BAM2FASTQFLAG=$bam2fqflag in the configuration file."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi

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
            mkdir -p $outputdir/logs
	elif [ $runmode != "batch" ]
        then
	    echo "resetting logs"
	    `rm -r $outputdir/logs/*`
        fi
	`chmod -R 770 $outputdir`
        `chmod 750 $epilogue`
        TopOutputLogs=$outputdir/logs
        pipeid=$( cat $TopOutputLogs/CONFIGUREpbs )


   # construct a list of SampleNames, check that files actually exist
   if [ ! -e $outputdir/SAMPLENAMES.list ]
   then
      truncate -s 0 $outputdir/SAMPLENAMES.tmp.list
      for fastqfile in $sampledir/*
      do
         # strip path, which read (left/right), and extension from input files
         # and put that info into the SampleNames file
         SampleName=$( basename $fastqfile | sed 's/_read.\?\.[^.]*$//' )
         echo -e "$SampleName" >> $outputdir/SAMPLENAMES.tmp.list
      done
      # paired-ended fastq will produce duplicate lines in the SampleNames file, so remove the duplicates
      uniq  $outputdir/SAMPLENAMES.tmp.list >  $outputdir/SAMPLENAMES.list
      sed -i '/^\s*$/d' $outputdir/SAMPLENAMES.list # remove blank lines
      rm  $outputdir/SAMPLENAMES.tmp.list
   fi
   # check that this actually worked, 
   # because otherwise the bash script will just go on, as if there is no problem
   if [ ! -s $outputdir/SAMPLENAMES.list ]
   then
      MSG="$outputdir/SAMPLENAMES.list is empty"
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
      exit 1;
   fi


   numsamples=`wc -l $outputdir/SAMPLENAMES.list | cut -d ' ' -f 1`
   if [ $numsamples -lt 1 ]
   then
      MSG="No samples found in SAMPLEDIR=$sampledir."
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      exit 1;
   fi




# generate a qsub header so we would not have to repeat the same lines
   generic_qsub_header=$outputdir/qsubGenericHeader
   truncate -s 0 $generic_qsub_header
   echo "#!/bin/bash" > $generic_qsub_header
   echo "#PBS -V" >> $generic_qsub_header
   echo "#PBS -A $pbsprj" >> $generic_qsub_header
   echo "#PBS -q $pbsqueue" >> $generic_qsub_header
   echo "#PBS -m ae" >> $generic_qsub_header
   echo "#PBS -M $email" >> $generic_qsub_header
   # check that this actually worked,
   # because otherwise the bash script will just go on, as if there is no problem
   if [ ! -s $generic_qsub_header ]
   then
      MSG="$generic_qsub_header is empty"
      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
      exit 1;
   fi



        
        case=""
        if [ $analysis == "ALIGN" -o $analysis == "ALIGNMENT" ]
        then
            echo "Type of analysis to run: ALIGNMENT only" 
            
            qsub1=$TopOutputLogs/qsub.start_align_block
            echo "#PBS -V" > $qsub1
            echo "#PBS -A $pbsprj" >> $qsub1
            echo "#PBS -N ${pipeid}_START_ALIGN_BLOCK" >> $qsub1
            echo "#pbs -l epilogue=$epilogue" >> $qsub1
	    echo "#PBS -l walltime=00:30:00" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $TopOutputLogs/start_align_block.ou" >> $qsub1
	    echo "#PBS -e $TopOutputLogs/start_align_block.in" >> $qsub1
            echo "#PBS -q $pbsqueue" >> $qsub1
            echo "#PBS -m ae" >> $qsub1
            echo "#PBS -M $email" >> $qsub1
            echo "$scriptdir/start_align_block.sh $runfile $TopOutputLogs/start_align_block.in $TopOutputLogs/start_align_block.ou $email $TopOutputLogs/qsub.start_align_block" >> $qsub1
            `chmod a+r $qsub1`               
            `qsub $qsub1 >> $TopOutputLogs/ALIGNpbs`
            echo `date`
            case="alignonly"
	fi
	if [ $analysis == "REALIGNONLY" -o $analysis == "REALIGN_ONLY" ]
        then
	    echo "Type of analysis to run: REALIGNMENT only. bams provided"
	    qsub2=$TopOutputLogs/qsub.start_realrecal_block
	    echo "#PBS -V" > $qsub2
	    echo "#PBS -A $pbsprj" >> $qsub2
	    echo "#PBS -N ${pipeid}_START_REALRECAL_BLOCK"
	    echo "#pbs -l epilogue=$epilogue" >> $qsub2
	    echo "#PBS -l walltime=00:30:00" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	    echo "#PBS -o $TopOutputLogs/start_realrecal_block.ou" >> $qsub2
	    echo "#PBS -e $TopOutputLogs/start_realrecal_block.in" >> $qsub2
	    echo "#PBS -q $pbsqueue" >> $qsub2
	    echo "#PBS -m ae" >> $qsub2
	    echo "#PBS -M $email" >> $qsub2
	    echo "$scriptdir/start_realrecal_block.sh $runfile $TopOutputLogs/start_realrecal_block.in $TopOutputLogs/start_realrecal_block.ou $email $TopOutputLogs/qsub.start_realrecal_block" >> $qsub2
	    `chmod a+r $qsub2` 
	    `qsub $qsub2 >> $TopOutputLogs/REALRECALpbs`
	    echo `date`
            case="realignonly" 
        fi
        if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ]
        then
	    echo "Type of analysis to run: ALIGNMENT and REALIGNMENT"
	    qsub1=$TopOutputLogs/qsub.start_align_block
	    echo "#PBS -V" > $qsub1
	    echo "#PBS -A $pbsprj" >> $qsub1
	    echo "#PBS -N ${pipeid}_START_ALIGN_BLOCK" >> $qsub1
	    echo "#pbs -l epilogue=$epilogue" >> $qsub1
	    echo "#PBS -l walltime=00:30:00" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $TopOutputLogs/start_align_block.ou" >> $qsub1
	    echo "#PBS -e $TopOutputLogs/start_align_block.in" >> $qsub1
	    echo "#PBS -q $pbsqueue" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "$scriptdir/start_align_block.sh $runfile $TopOutputLogs/start_align_block.in $TopOutputLogs/start_align_block.ou $email $TopOutputLogs/qsub.main.aln" >> $qsub1
	    `chmod a+r $qsub1`               
	    `qsub $qsub1 >> $TopOutputLogs/ALIGNpbs`
	    echo `date`
            echo "Note: realign module will be scheduled after align module ends"
            case="align and realign"  
        fi
        if [ $analysis == "VCALL_ONLY" -o $analysis == "VCALL" ]
        then

            echo "variant calling only"
	    qsub3=$TopOutputLogs/qsub.start_varcall_block
	    echo "#PBS -V" > $qsub3
	    echo "#PBS -A $pbsprj" >> $qsub3
	    echo "#PBS -N ${pipeid}_START_VARCALL_BLOCK" >> $qsub3
	    echo "#PBS -l epilogue=$epilogue" >> $qsub3
	    echo "#PBS -l walltime=00:30:00" >> $qsub3
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub3
	    echo "#PBS -o $TopOutputLogs/log.start_varcall_block.ou" >> $qsub3
	    echo "#PBS -e $TopOutputLogs/log.start_varcall_block.in" >> $qsub3
	    echo "#PBS -q $pbsqueue" >> $qsub3
	    echo "#PBS -m ae" >> $qsub3
	    echo "#PBS -M $email" >> $qsub3
	    echo "$scriptdir/start_varcall_block.sh $runfile $TopOutputLogs/log.start_varcall_block.in $TopOutputLogs/log.start_varcall_block.ou $email $TopOutputLogs/qsub.start_varcall_block" >> $qsub3
	    `chmod a+r $qsub3`
	    vcalljobid=`qsub $qsub3`
	    echo $vcalljobid >> $TopOutputLogs/VARCALLpbs
            case="vcall_only"  
        fi
        if [ `expr ${#case}` -lt 1 ]
        then
	       MSG="Invalid value for parameter ANALYSIS=$analysis in configuration file."
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
               exit 1; 
       fi
fi