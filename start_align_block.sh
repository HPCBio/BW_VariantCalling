#!/bin/sh
#
# written in collaboration with Mayo Bioinformatics core group
# start_align_block.sh
# First module in the GGPS analysis pipeline
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
else 
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
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        ## parsing run info file

        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
	inputdir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
        nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        aligner=$( cat $runfile | grep -w ALIGNER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        igvdir=$( cat $runfile | grep -w IGVDIR | cut -d '=' -f2 )
        fastqcdir=$( cat $runfile | grep -w FASTQCDIR | cut -d '=' -f2 )
        fastqcflag=$( cat $runfile | grep -w FASTQCFLAG | cut -d '=' -f2 )
        fastqcparms=$( cat $runfile | grep -w FASTQCPARMS | cut -d '=' -f2 | tr " " "_" )
        bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        bamtofastqparms=$( cat $runfile | grep -w BAM2FASTQPARMS | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE  | cut -d '=' -f2 )
        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )

        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        sortool=$( cat $runfile | grep -w SORTMERGETOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

        if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
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

        if [ $aligner != "NOVOALIGN" -a $aligner != "BWA" -a $aligner != "BWA_MEM" ]
        then
            MSG="ALIGNER=$aligner  is not available at this site"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi

        if [ -z $epilogue ]
        then
           MSG="Value for EPILOGUE must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
        fi

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        #deprecated
        #if [ ! -s $samplefileinfo ]
        #then
        #   MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
        #   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
        #   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
        #   exit 1;
        #fi

        if [ ! -d $outputdir ]
        then
           mkdir -p $outputdir
        fi

        if [ ! -d $alignerdir ]
        then
           MSG="$alignerdir aligner directory not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        #resetting output directories, logs, files

        oualigndir=$outputdir/align
        output_logs=$outputdir/logs
        pipeid=$( cat $output_logs/CONFIGUREpbs )
        
        if [ -d $oualigndir ]
        then
           echo "$oualigndir is there; resetting it"
           `rm -r $oualigndir/*`
        else
           mkdir -p $oualigndir
        fi

        if [ -d $output_logs ]
        then
           echo "$output_logs is there; resetting it"
           #`rm -r $output_logs/*`
           pbsids=""
        else
           mkdir -p $output_logs
        fi

        #################################
        # QC preprocessing on the samplefileinfo file
        # it verifies that the samplefileinfo file contains the correct information
        #################################

       # WE DO NOT USE SAMPLEDETAIL FILE ANYMORE,INSTEAD GRABBING SAMPLE FILES FROM A GIVEN FOLDER
       # WHICH IS WHY THIS BLOCK DOES NOT WORK
       # AND AS A RESULT ONLY INPUT FASTQ CAN BE PROCESSED
       # MAY NEED TO RECUSCITATE THIS IN A MODIFIED FORM IN ORDER TO BE ABLE TO PROCESS BAMS
       # ---  LSM April 22, 2014

        while read sampledetail 
        do
          if [ `expr ${#sampledetail}` -lt 5 ]
            then
              echo "skipping empty line"
            else
              echo "reading next line: $sampledetail"
              sampleTag=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f1 )
              samplefile=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f2 )
              R1=$( echo $sampledetail | cut -d '=' -f2 | cut -d ' ' -f1 )
              R2=$( echo $sampledetail | cut -d ' ' -f2 )

              if [ `expr ${#inputformat}` -lt 1 ]
              then
		MSG="Invalid input format. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi

              if [ $inputformat != "BAM" -a $inputformat != "FASTQ" ]
              then
		MSG="Invalid input format. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi
 
              if [ `expr ${#sampleTag}` -lt 1  ]
              then
		MSG="Invalid samplename. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi
 
              if [ `expr ${#samplefile}` -lt 1  ]
              then
		MSG="inputfile was not specified. parsing of line in SAMPLEFILENAMES failed. alignment stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
              fi
	      if [ $inputformat == "BAM" ]
	      then
                  if [ ! -s $samplefile ]
                  then
		      MSG="$samplefile input file specified in SAMPLEFILENAMES was not found. alignment stopped"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		      exit 1;
		  else
		      totlines=`tail $samplefile | wc -l`
		      if [ $totlines -lt 1 ]
		      then
			  MSG="$samplefile input file specified in SAMPLEFILENAMES is empty. alignment stopped"
			  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
			  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			  exit 1;
                      fi
                  fi
              fi
	      if [ $inputformat == "FASTQ" ]
              then 
                  if [ ! -s $R1 ]
                  then
		      MSG="$R1 input file specified in SAMPLEFILENAMES was not found. alignment stopped"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		      exit 1;
		  else
		      if [ `tail $R1 | wc -l` -lt 1 ]
		      then
			  MSG="$R1 input file specified in SAMPLEFILENMAES is empty. alignment stopped"
			  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
			  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			  exit 1;
                      fi
                  fi
                  if [ $paired -eq 1 ]
                  then
                      if [ ! -s $R2 ]
                      then 
			  MSG="$R2 input file specified in SAMPLEFILENAMES was not found. alignment stopped"
			  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
			  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			  exit 1;
		      else
			  if [ `tail $R2 | wc -l` -lt 1 ]
			  then
			      MSG="$R2 input file specified in SAMPLEFILENAMES is empty. alignment stopped"
			      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
			      #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			      exit 1;
                          fi
                      fi
                  fi
              fi
          fi
          echo "bottom of the loop"
        done < $samplefileinfo

        #################################
        # PREPROCESSING BAMS
        #################################

        CONVERT=""
        case=""
        typeOfupdateconfig=""

	if [ $inputformat == "BAM" ]
        then
            echo "input files are BAMs; preprocessing is required before aligning files"
            newfqfiles=""
            newbamfiles=""
            sep=":"

            if [ `expr ${#inputdir}` -lt 1 ]
            then
		MSG="a path must be specified in INPUTDIR when input files are BAMS"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
                exit 1;
            else
		if [ ! -d $inputdir ]
		then
		    mkdir -p $inputdir
		fi
            fi
           
            while read sampledetail
            do
		if [ `expr ${#sampledetail}` -lt 7 ]
		then
                    echo "skipping empty line"
		else
                    echo "Preprocessing step with BAMS"
		    originalBAM=$( echo $sampledetail | cut -d ':' -f2 | cut -d '=' -f2 )
		    suffix=$( echo $RANDOM )
		    tmpconversion=$inputdir/$suffix
                    # weird behavior with this next folder

		    if [ ! -d $tmpconversion ]
		    then
			`mkdir $tmpconversion`
                        exitcode=$?
                    else
                        `rm -r $tmpconversion/*`
                        exitcode=$?
		    fi

                    if [ $exitcode -ne 0 ]
                    then
			MSG="$tmpconversion tempdir for bam conversion was not created"
			echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
			#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
			exit 1;
                    fi
                    newsuffix4fq=${suffix}.preprocessed
                    newsuffix4bam=${suffix}.preprocessed
                    newfqfiles=${newsuffix4fq}${sep}${newfqfiles}
                    newbamfiles=${newsuffix4bam}${sep}${newbamfiles}


		    if [ $bamtofastqflag == "NO" ] 
                    then
                        echo "bam2newbam preprocessing..."
			typeOfupdateconfig="bam2newbam"
                        
                        #########################
                        # TODO: new script goes here
                        # this script performs: namesort, remove singletons, revertsam (if so indicated)
                        # input is originalBAM, output is newsuffix4bam
                        #########################

			qsub1=$output_logs/qsub.convertbam2newbam.$suffix
			echo "#PBS -V" > $qsub1
			echo "#PBS -A $pbsprj" >> $qsub1
			echo "#PBS -N ${pipeid}_convertbam2newbam_${suffix}" >> $qsub1
			echo "#pbs -l epilogue=$epilogue" >> $qsub1
			echo "#PBS -l walltime=$pbscpu" >> $qsub1
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
			echo "#PBS -o $output_logs/log.convertbam2newbam.${suffix}.ou" >> $qsub1
			echo "#PBS -e $output_logs/log.convertbam2newbam.${suffix}.in" >> $qsub1
			echo "#PBS -q $pbsqueue" >> $qsub1
			echo "#PBS -m ae" >> $qsub1
			echo "#PBS -M $email" >> $qsub1
			echo "aprun -n 1 -d $thr $scriptdir/convertbam2newbam.sh $inputdir $tmpconversion $originalBAM $newsuffix4bam $runfile $output_logs/log.convertbam2newbam.${suffix}.in $output_logs/log.convertbam2newbam.${suffix}.ou $email $output_logs/qsub.convertbam2newbam.$suffix" >> $qsub1
			`chmod a+r $qsub1` 
			convertjob=`qsub $qsub1`
			`qhold -h u $convertjob` 
			echo $convertjob >> $output_logs/CONVERTBAMpbs
			echo `date`

                    else
                        echo "bam2fastq preprocessing..."
			typeOfupdateconfig="bam2fastq"

                        #########################                        
                        # TODO: add preprocessing stuff to convertbam.sh
                        # the new convertbam.sh  script must perform: 
                        # namesort, remove singletons, revertsam (if so indicated, bam2fastq
                        # input is originalBAM, output is newsuffix4fq
                        #########################

			qsub1=$output_logs/qsub.convertbam2fq.$suffix
			echo "#PBS -V" > $qsub1
			echo "#PBS -A $pbsprj" >> $qsub1
			echo "#PBS -N ${pipeid}_convertbam2fq_${suffix}" >> $qsub1
			echo "#pbs -l epilogue=$epilogue" >> $qsub1
			echo "#PBS -l walltime=$pbscpu" >> $qsub1
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
			echo "#PBS -o $output_logs/log.convertbam2fq.${suffix}.ou" >> $qsub1
			echo "#PBS -e $output_logs/log.convertbam2fq.${suffix}.in" >> $qsub1
			echo "#PBS -q $pbsqueue" >> $qsub1
			echo "#PBS -m ae" >> $qsub1
			echo "#PBS -M $email" >> $qsub1
			echo "aprun -n 1 -d $thr $scriptdir/convertbam2fastq.sh $inputdir $tmpconversion $originalBAM $newsuffix4fq $runfile $output_logs/log.convertbam2fq.${suffix}.in $output_logs/log.convertbam2fq.${suffix}.ou $email $output_logs/qsub.convertbam2fq.$suffix" >> $qsub1
			`chmod a+r $qsub1` 
			convertjob=`qsub $qsub1`
			`qhold -h u $convertjob` 
			echo $convertjob >> $output_logs/CONVERTBAMpbs
			echo `date`
                    fi
		fi
	    done < $samplefileinfo
            echo `date`

            ###########################
            ## updating config files
            ###########################
        
            CONVERTids=$( cat $output_logs/CONVERTBAMpbs | sed "s/\.[a-z]*//" | tr "\n" ":" )

	    if [ $input_typeOfupdateconfig == "bam2fastq" ]
            then
                ###############################
                # todo:
                # rename updateconfig script for newfastq files
                # rename updateconfig.sh to updateconfig.wnewfq.sh
                ###############################

		qsub2=$output_logs/qsub.updateconfig_wnewfq
		echo "#PBS -V" > $qsub2
		echo "#PBS -A $pbsprj" >> $qsub2
		echo "#PBS -N ${pipeid}_updateconfig_wnewfq" >> $qsub2
		echo "#pbs -l epilogue=$epilogue" >> $qsub2
		echo "#PBS -l walltime=$pbscpu" >> $qsub2
		echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		echo "#PBS -o $output_logs/log.updateconfig_wnewfq.ou" >> $qsub2
		echo "#PBS -e $output_logs/log.updateconfig_wnewfq.in" >> $qsub2
		echo "#PBS -q $pbsqueue" >> $qsub2
		echo "#PBS -m ae" >> $qsub2
		echo "#PBS -M $email" >> $qsub2
		echo "#PBS -W depend=afterok:$CONVERTids" >> $qsub2
		echo "aprun -n 1 -d 1 $scriptdir/updateconfig.wnewfq.sh $inputdir $newfqfiles $runfile $samplefileinfo $output_logs/log.updateconfig_wnewfq.in $output_logs/log.updateconfig_wnewfq.ou $email $output_logs/qsub.updateconfig_wnewfq" >> $qsub2
		`chmod a+r $qsub2`       
		updatejob=`qsub $qsub2` 
		echo $updatejob >> $output_logs/UPDATECONFIGpbs
            else
               if [ $input_typeOfupdateconfig == "bam2newbam" ]
               then
                   ###############################
                   # todo:
                   # create new script here for updateconfig for newbams
                   # updateconfig.wnewbams.sh
                   ###############################

		   qsub3=$output_logs/qsub.updateconfig_wnewbam
		   echo "#PBS -V" > $qsub3
		   echo "#PBS -A $pbsprj" >> $qsub3
		   echo "#PBS -N ${pipeid}_updateconfig_wnewbam" >> $qsub3
		   echo "#pbs -l epilogue=$epilogue" >> $qsub3
		   echo "#PBS -l walltime=$pbscpu" >> $qsub3
		   echo "#PBS -l nodes=1:ppn=1" >> $qsub3
		   echo "#PBS -o $output_logs/log.updateconfig_wnewbam.ou" >> $qsub3
		   echo "#PBS -e $output_logs/log.updateconfig_wnewbam.in" >> $qsub3
		   echo "#PBS -q $pbsqueue" >> $qsub3
		   echo "#PBS -m ae" >> $qsub3
		   echo "#PBS -M $email" >> $qsub3
		   echo "#PBS -W depend=afterok:$CONVERTids" >> $qsub3
		   echo "aprun -n 1 -d 1 $scriptdir/updateconfig.wnewbam.sh $inputdir $newbamfiles $runfile $samplefileinfo $output_logs/log.updateconfig_wnewbam.in $output_logs/log.updateconfig_wnewbam.ou $email $output_logs/qsub.updateconfig_wnewbam" >> $qsub3
		   `chmod a+r $qsub3`       
		   updatejob=`qsub $qsub3` 
		   echo $updatejob >> $output_logs/UPDATECONFIGpbs

               fi
            fi

            allconjobs=$( echo $CONVERTids | tr ":" " " )
            `qrls -h u $allconjobs`
	    echo `date`

        fi

        #################################
        # PREPROCESSING OF BAMS --- DONE
        #################################
        # ALIGNMENT PROPER CAN PROCEED NOW
        #################################

        if [ $bamtofastqflag == "YES" -a $inputformat == "BAM" ]
        then
	    qsub1=$output_logs/qsub.main.alnFQ.afterbam2fastq
	    echo "#PBS -V" > $qsub1
	    echo "#PBS -A $pbsprj" >> $qsub1
	    echo "#PBS -N ${pipeid}_alnFQ_afterbam2fastq" >> $qsub1
	    echo "#pbs -l epilogue=$epilogue" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $output_logs/alnFQ.afterbam2fastq.ou" >> $qsub1
	    echo "#PBS -e $output_logs/alnFQ.afterbam2fastq.in" >> $qsub1
	    echo "#PBS -q $pbsqueue" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "#PBS -W depend=afterok:$updatejob" >> $qsub1
	    echo "$scriptdir/alignfastq.sh $runfile $output_logs/alnFQ.afterbam2fastq.in $output_logs/alnFQ.afterbam2fastq.ou $email $output_logs/qsub.main.alnFQ.afterbam2fastq" >> $qsub1
	    `chmod a+r $qsub1`               
	    `qsub $qsub1 >> $output_logs/ALIGNpbs`
            case="bam2fastq"
	    echo `date`
	fi

        if [ $bamtofastqflag == "NO" -a $inputformat == "BAM" ]
        then
            echo "aligning bam files directly"


            ####################################
            #  TODO:
            #  aligbam.sh does not need to call revertsam
            ####################################

            qsub2=$output_logs/qsub.alignbams
            echo "#PBS -V" > $qsub2
            echo "#PBS -A $pbsprj" >> $qsub2
            echo "#PBS -N ${pipeid}_alignbam" >> $qsub2
            echo "#PBS -l epilogue=$epilogue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	    echo "#PBS -o $output_logs/log.alignbams.ou" >> $qsub2
	    echo "#PBS -e $output_logs/log.alignbams.in" >> $qsub2
            echo "#PBS -q $pbsqueue" >> $qsub2
            echo "#PBS -m ae" >> $qsub2
            echo "#PBS -M $email" >> $qsub2
	    echo "#PBS -W depend=afterok:$updatejob" >> $qsub2
            echo "$scriptdir/alignbam.sh $runfile $output_logs/log.alignbams.in $output_logs/log.alignbams.ou $email $output_logs/qsub.alignbams" >> $qsub2
            `chmod a+r $qsub2`               
            `qsub $qsub2 >> $output_logs/ALIGNpbs`
            case="alignbams"
            echo `date`
        fi

        if [ $bamtofastqflag == "NO" -a $inputformat == "FASTQ" ]
        then
            echo "aligning fastq files directly"
            qsub3=$output_logs/qsub.alignfastq
            echo "#PBS -V" > $qsub3
            echo "#PBS -A $pbsprj" >> $qsub3
            echo "#PBS -N ${pipeid}_alignfastq" >> $qsub3
            echo "#PBS -l epilogue=$epilogue" >> $qsub3
	    echo "#PBS -l walltime=$pbscpu" >> $qsub3
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub3
	    echo "#PBS -o $output_logs/log.alignfastq.ou" >> $qsub3
	    echo "#PBS -e $output_logs/log.alignfastq.in" >> $qsub3
            echo "#PBS -q $pbsqueue" >> $qsub3
            echo "#PBS -m ae" >> $qsub3
            echo "#PBS -M $email" >> $qsub3
            echo "$scriptdir/alignfastq.sh $runfile $output_logs/log.alignfastq.in $output_logs/log.alignfastq.ou $email $output_logs/qsub.alignfastq" >> $qsub3
            `chmod a+r $qsub3` 
            `qsub $qsub3 >> $output_logs/ALIGNpbs`
            case="alignfastq"
            echo `date`
        fi

        if [ `expr ${#case}` -lt 1 ]
        then
           MSG="Alignment module failed to launch. Incompatible values specified in config files bam2fastqflag=$bamtofastqflag input_format=$inputformat analysis=$analysis"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi         
fi
