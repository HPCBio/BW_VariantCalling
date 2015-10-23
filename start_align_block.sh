#!/bin/bash
#
# start_align_block.sh
# First module in the GGPS analysis pipeline
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
else 

        echo -e "\n\n############# START ALIGN BLOCK ###############\n\n" >&2
        umask 0037
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
           exit 1;
        fi

# wrapping commends in echo, so that the output logs would be easier to read: they will have more structure
set +x
echo -e "\n\n####################################################################################################" >&2
echo        "##################################### PARSING RUN INFO FILE ########################################" >&2
echo        "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "\n\n####################################################################################################\n\n" >&2; set -x;


        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
	sampledir=$( cat $runfile | grep -w INPUTDIR | cut -d '=' -f2 )
        nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        aligner=$( cat $runfile | grep -w ALIGNER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        bamtofastqparms=$( cat $runfile | grep -w BAM2FASTQPARMS | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE  | cut -d '=' -f2 )
        input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
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
                echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1;
            fi
        fi


        if [ $inputformat != "FASTQ" -a $inputformat != "BAM" ]
        then
            MSG="Incorrect value for INPUTFORMAT=$inputformat in the configuration file."
            echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi


        if [ $bamtofastqflag != "YES" -a $bamtofastqflag != "NO" -a $bamtofastqflag != "1" -a $bamtofastqflag != "0" ]
        then
            MSG="BAM2FASTQFLAG=$bamtofastqflag  invalid value"
            echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
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

        if [ $aligner != "NOVOALIGN" -a $aligner != "BWA_ALN" -a $aligner != "BWA_MEM" ]
        then
            MSG="ALIGNER=$aligner  is not available at this site"
            echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit 1;
        fi

        if [ -z $epilogue ]
        then
           MSG="Value for EPILOGUE must be specified in configuration file"
           echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        else
           `chmod 740 $epilogue`
        fi

        if [ ! -d $scriptdir ]
        then
           MSG="SCRIPTDIR=$scriptdir directory not found"
           echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi

        if [ ! -d $outputdir ]
        then
           mkdir -p $outputdir
        fi


        set +x; echo -e "\n\n#####################    resetting output directories, logs, files     #############################\n\n" >&2; set -x; 


        #AlignmentOutputFolder=$outputdir/align
        TopLogsFolder=$outputdir/logs
        
        #if [ -d $AlignmentOutputFolder ]
        #then
        #   echo "$AlignmentOutputFolder is there; resetting it"
        #   `rm -r $AlignmentOutputFolder/*`
        #else
        #   mkdir -p $AlignmentOutputFolder
        #fi

        pipeid=$( cat $TopLogsFolder/pbs.CONFIGURE )
        pbsids=""



        set +x; echo -e "\n\n" >&2
        echo "####################################################################" >&2
        echo "####################################################################" >&2
        echo "############      preprocessing BAM files         ##################" >&2
        echo "####################################################################" >&2
        echo "####################################################################" >&2
        echo -e "\n\n" >&2; set -x

        CONVERT=""
        case=""
        typeOfupdateconfig=""
        truncate -s 0 $TopLogsFolder/pbs.CONVERTBAM

	if [ $inputformat == "BAM" ]
        then
            set +x; echo -e "\n #input files are BAMs; preprocessing is required before aligning files \n" >&2; set -x;
            newfqfiles=""
            newbamfiles=""
            sep=":"

           set +x; echo -e "\n\n #### Begin loop over samples! #### \n\n" >&2; set -x; 
           while read SampleName
           do
              set +x; echo -e "\n # processing next sample\n" >&2; set -x;
              # this will evaluate the length of string; 
              if [ `expr ${#SampleName}` -lt 1 ]
              then
                echo "skipping empty line"
              else

                 set +x; echo -e "\n # processing: $SampleName \n" >&2; set -x;
                 # make sure this is actually BAM
                 if [ ! -f ${sampledir}/${SampleName}.bam ]
                 then
                    MSG="${sampledir}/${SampleName}.bam does not exist"
                    echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                    exit 1;
                 else
           
		    originalBAM=${sampledir}/${SampleName}.bam
		    suffix=$( echo $RANDOM )
		    tmpconversion=$sampledir/$suffix
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
                        echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
			exit 1;
                    fi
                    newsuffix4fq=${suffix}.preprocessed
                    newsuffix4bam=${suffix}.preprocessed
                    newfqfiles=${newsuffix4fq}${sep}${newfqfiles}
                    newbamfiles=${newsuffix4bam}${sep}${newbamfiles}


		    if [ $bamtofastqflag == "NO" ] 
                    then
                        set +x; echo -e "\n # bam2newbam preprocessing... \n" >&2; set -x;
			typeOfupdateconfig="bam2newbam"
                        
                        #########################
                        # TODO: new script goes here
                        # this script performs: namesort, remove singletons, revertsam (if so indicated)
                        # input is originalBAM, output is newsuffix4bam
                        #########################

			qsub1=$TopLogsFolder/qsub.convertbam2newbam.$suffix
			echo "#PBS -V" > $qsub1
			echo "#PBS -A $pbsprj" >> $qsub1
			echo "#PBS -N ${pipeid}_convertbam2newbam_${suffix}" >> $qsub1
			echo "#pbs -l epilogue=$epilogue" >> $qsub1
			echo "#PBS -l walltime=$pbscpu" >> $qsub1
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
			echo "#PBS -o $TopLogsFolder/log.convertbam2newbam.${suffix}.ou" >> $qsub1
			echo "#PBS -e $TopLogsFolder/log.convertbam2newbam.${suffix}.in" >> $qsub1
			echo "#PBS -q $pbsqueue" >> $qsub1
			echo "#PBS -m ae" >> $qsub1
			echo "#PBS -M $email" >> $qsub1
			echo "aprun -n 1 -d $thr $scriptdir/convertbam2newbam.sh $sampledir $tmpconversion $originalBAM $newsuffix4bam $runfile $TopLogsFolder/log.convertbam2newbam.${suffix}.in $TopLogsFolder/log.convertbam2newbam.${suffix}.ou $email $TopLogsFolder/qsub.convertbam2newbam.$suffix" >> $qsub1
                        echo -e "\n\n" >> $qsub1
                        echo "exitcode=\$?" >> $qsub1
                        echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub1
                        echo "   echo -e \"\n\n convertbam2newbam.sh failed with exit code = \$exitcode \n logfile=$TopLogsFolder/log.convertbam2newbam.${suffix}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub1
                        echo "fi" >> $qsub1
                        echo -e "\n\n exit 1" >> $qsub1

			#`chmod a+r $qsub1` 
			convertjob=`qsub $qsub1`
			`qhold -h u $convertjob` 
			echo $convertjob >> $TopLogsFolder/pbs.CONVERTBAM
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

			qsub1=$TopLogsFolder/qsub.convertbam2fq.$suffix
			echo "#PBS -V" > $qsub1
			echo "#PBS -A $pbsprj" >> $qsub1
			echo "#PBS -N ${pipeid}_convertbam2fq_${suffix}" >> $qsub1
			echo "#pbs -l epilogue=$epilogue" >> $qsub1
			echo "#PBS -l walltime=$pbscpu" >> $qsub1
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
			echo "#PBS -o $TopLogsFolder/log.convertbam2fq.${suffix}.ou" >> $qsub1
			echo "#PBS -e $TopLogsFolder/log.convertbam2fq.${suffix}.in" >> $qsub1
			echo "#PBS -q $pbsqueue" >> $qsub1
			echo "#PBS -m ae" >> $qsub1
			echo "#PBS -M $email" >> $qsub1
			echo "aprun -n 1 -d $thr $scriptdir/convertbam2fastq.sh $sampledir $tmpconversion $originalBAM $newsuffix4fq $runfile $TopLogsFolder/log.convertbam2fq.${suffix}.in $TopLogsFolder/log.convertbam2fq.${suffix}.ou $email $TopLogsFolder/qsub.convertbam2fq.$suffix" >> $qsub1

                        echo "exitcode=\$?" >> $qsub1
                        echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub1
                        echo "   echo -e \"\n\n convertbam2fastq.sh failed with exit code = \$exitcode \n logfile=$TopLogsFolder/log.convertbam2fq.${suffix}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub1
                        echo "fi" >> $qsub1
                        echo -e "\n\n exit 1" >> $qsub1

			#`chmod a+r $qsub1` 
			convertjob=`qsub $qsub1`
			`qhold -h u $convertjob` 
			echo $convertjob >> $TopLogsFolder/pbs.CONVERTBAM
			echo `date`
                    fi
		fi
             fi 
	    done < $outputdir/SAMPLENAMES.list
            set +x; echo -e "\n\n #### End loop over samples! #### \n\n" >&2; set -x; 
            echo `date`




            set +x; echo -e "\n\n ########################### " >&2
                    echo -e " ## updating config files " >&2
                    echo -e " ########################### \n\n" >&2; set -x;
        
            CONVERTids=$( cat $TopLogsFolder/pbs.CONVERTBAM | sed "s/\..*//" | tr "\n" ":" )

	    if [ $bamtofastqflag == "YES" ]
            then
                ###############################
                # todo:
                # rename updateconfig script for newfastq files
                # rename updateconfig.sh to updateconfig.wnewfq.sh
                ###############################

		qsub2=$TopLogsFolder/qsub.updateconfig_wnewfq
		echo "#PBS -V" > $qsub2
		echo "#PBS -A $pbsprj" >> $qsub2
		echo "#PBS -N ${pipeid}_updateconfig_wnewfq" >> $qsub2
		echo "#pbs -l epilogue=$epilogue" >> $qsub2
		echo "#PBS -l walltime=$pbscpu" >> $qsub2
		echo "#PBS -l nodes=1:ppn=1" >> $qsub2
		echo "#PBS -o $TopLogsFolder/log.updateconfig_wnewfq.ou" >> $qsub2
		echo "#PBS -e $TopLogsFolder/log.updateconfig_wnewfq.in" >> $qsub2
		echo "#PBS -q $pbsqueue" >> $qsub2
		echo "#PBS -m ae" >> $qsub2
		echo "#PBS -M $email" >> $qsub2
		echo "#PBS -W depend=afterok:$CONVERTids" >> $qsub2
		echo "aprun -n 1 -d 1 $scriptdir/updateconfig.wnewfq.sh $sampledir $newfqfiles $runfile $samplefileinfo $TopLogsFolder/log.updateconfig_wnewfq.in $TopLogsFolder/log.updateconfig_wnewfq.ou $email $TopLogsFolder/qsub.updateconfig_wnewfq" >> $qsub2

                echo "exitcode=\$?" >> $qsub2
                echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub2
                echo "   echo -e \"\n\n updateconfig.wnewfq.sh failed with exit code = \$exitcode \n logfile=$TopLogsFolder/log.updateconfig_wnewfq.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub2
                echo "fi" >> $qsub2
                echo -e "\n\n exit 1" >> $qsub2

		#`chmod a+r $qsub2`       
		updatejob=`qsub $qsub2` 
		echo $updatejob >> $TopLogsFolder/pbs.UPDATECONFIG
            else
                   ###############################
                   # todo:
                   # create new script here for updateconfig for newbams
                   # updateconfig.wnewbams.sh
                   ###############################

		   qsub3=$TopLogsFolder/qsub.updateconfig_wnewbam
		   echo "#PBS -V" > $qsub3
		   echo "#PBS -A $pbsprj" >> $qsub3
		   echo "#PBS -N ${pipeid}_updateconfig_wnewbam" >> $qsub3
		   echo "#pbs -l epilogue=$epilogue" >> $qsub3
		   echo "#PBS -l walltime=$pbscpu" >> $qsub3
		   echo "#PBS -l nodes=1:ppn=1" >> $qsub3
		   echo "#PBS -o $TopLogsFolder/log.updateconfig_wnewbam.ou" >> $qsub3
		   echo "#PBS -e $TopLogsFolder/log.updateconfig_wnewbam.in" >> $qsub3
		   echo "#PBS -q $pbsqueue" >> $qsub3
		   echo "#PBS -m ae" >> $qsub3
		   echo "#PBS -M $email" >> $qsub3
		   echo "#PBS -W depend=afterok:$CONVERTids" >> $qsub3
		   echo "aprun -n 1 -d 1 $scriptdir/updateconfig.wnewbam.sh $sampledir $newbamfiles $runfile $samplefileinfo $TopLogsFolder/log.updateconfig_wnewbam.in $TopLogsFolder/log.updateconfig_wnewbam.ou $email $TopLogsFolder/qsub.updateconfig_wnewbam" >> $qsub3

                   echo "exitcode=\$?" >> $qsub3
                   echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub3
                   echo "   echo -e \"\n\n updateconfig.wnewbam.sh failed with exit code = \$exitcode \n logfile=$TopLogsFolder/log.updateconfig_wnewbam.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub3
                   echo "fi" >> $qsub3
                   echo -e "\n\n exit 1" >> $qsub3


		   #`chmod a+r $qsub3`       
		   updatejob=`qsub $qsub3` 
		   echo $updatejob >> $TopLogsFolder/pbs.UPDATECONFIG
            fi

            allconjobs=$( echo $CONVERTids | tr ":" " " )
            `qrls -h u $allconjobs`
	    echo `date`

        fi

        set +x; echo -e "\n\n" >&2
        echo "####################################################################" >&2
        echo "############    done with preprocessing of BAMs   ################" >&2
        echo "############    alignment case selection is next    ################" >&2
        echo "####################################################################" >&2
        echo -e "\n\n" >&2; set -x;

        if [ $bamtofastqflag == "YES" -a $inputformat == "BAM" ]
        then
            set +x; echo -e "\n# input is BAM, convert to fastq\n" >&2;  set -x;
	    qsub1=$TopLogsFolder/qsub.main.alnFQ.afterbam2fastq
	    echo "#PBS -V" > $qsub1
	    echo "#PBS -A $pbsprj" >> $qsub1
	    echo "#PBS -N ${pipeid}_alnFQ_afterbam2fastq" >> $qsub1
	    echo "#pbs -l epilogue=$epilogue" >> $qsub1
	    echo "#PBS -l walltime=$pbscpu" >> $qsub1
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	    echo "#PBS -o $TopLogsFolder/alnFQ.afterbam2fastq.ou" >> $qsub1
	    echo "#PBS -e $TopLogsFolder/alnFQ.afterbam2fastq.in" >> $qsub1
	    echo "#PBS -q $pbsqueue" >> $qsub1
	    echo "#PBS -m ae" >> $qsub1
	    echo "#PBS -M $email" >> $qsub1
	    echo "#PBS -W depend=afterok:$updatejob" >> $qsub1
	    echo "$scriptdir/alignfastq.sh $runfile $TopLogsFolder/alnFQ.afterbam2fastq.in $TopLogsFolder/alnFQ.afterbam2fastq.ou $email $TopLogsFolder/qsub.main.alnFQ.afterbam2fastq" >> $qsub1

	    #`chmod a+r $qsub1`               
	    `qsub $qsub1 >> $TopLogsFolder/pbs.ALIGN`
            case="bam2fastq"
	    echo `date`
	fi

        if [ $bamtofastqflag == "NO" -a $inputformat == "BAM" ]
        then
            set +x; echo "# aligning bam files directly\n" >&2;  set -x;


            ####################################
            #  TODO:
            #  aligbam.sh does not need to call revertsam
            ####################################

            qsub2=$TopLogsFolder/qsub.alignbams
            echo "#PBS -V" > $qsub2
            echo "#PBS -A $pbsprj" >> $qsub2
            echo "#PBS -N ${pipeid}_alignbam" >> $qsub2
            echo "#PBS -l epilogue=$epilogue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	    echo "#PBS -o $TopLogsFolder/log.alignbams.ou" >> $qsub2
	    echo "#PBS -e $TopLogsFolder/log.alignbams.in" >> $qsub2
            echo "#PBS -q $pbsqueue" >> $qsub2
            echo "#PBS -m ae" >> $qsub2
            echo "#PBS -M $email" >> $qsub2
	    echo "#PBS -W depend=afterok:$updatejob" >> $qsub2
            echo "$scriptdir/alignbam.sh $runfile $TopLogsFolder/log.alignbams.in $TopLogsFolder/log.alignbams.ou $email $TopLogsFolder/qsub.alignbams" >> $qsub2
            #`chmod a+r $qsub2`               
            `qsub $qsub2 >> $TopLogsFolder/pbs.ALIGN`
            case="alignbams"
            echo `date`
        fi

        if [ $bamtofastqflag == "NO" -a $inputformat == "FASTQ" ]
        then
            set +x; echo -e "\n# aligning fastq files directly\n" >&2;  set -x;
            qsub3=$TopLogsFolder/qsub.alignfastq
            echo "#PBS -V" > $qsub3
            echo "#PBS -A $pbsprj" >> $qsub3
            echo "#PBS -N ${pipeid}_alignfastq" >> $qsub3
            echo "#PBS -l epilogue=$epilogue" >> $qsub3
	    echo "#PBS -l walltime=$pbscpu" >> $qsub3
	    echo "#PBS -l nodes=1:ppn=1" >> $qsub3
	    echo "#PBS -o $TopLogsFolder/log.alignfastq.ou" >> $qsub3
	    echo "#PBS -e $TopLogsFolder/log.alignfastq.in" >> $qsub3
            echo "#PBS -q $pbsqueue" >> $qsub3
            echo "#PBS -m ae" >> $qsub3
            echo "#PBS -M $email" >> $qsub3
            echo "$scriptdir/alignfastq.sh $runfile $TopLogsFolder/log.alignfastq.in $TopLogsFolder/log.alignfastq.ou $email $TopLogsFolder/qsub.alignfastq" >> $qsub3
            #`chmod a+r $qsub3` 
            `qsub $qsub3 >> $TopLogsFolder/pbs.ALIGN`
            case="alignfastq"
            echo `date`
        fi

        if [ `expr ${#case}` -lt 1 ]
        then
           MSG="Alignment module failed to launch. Incompatible values specified in config files bam2fastqflag=$bamtofastqflag inputformat=$inputformat analysis=$analysis"
           echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi         
fi

echo -e "now we need to make the PBS log files group readable"

find $outputdir -name logs -type d | awk '{print "chmod -R g+r "$1}' | sh -x
