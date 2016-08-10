#!/bin/bash
#
# start_align_block.sh
# First module in the GGPS analysis pipeline
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu

if [ $# != 5 ]
then
        MSG="Parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi 

echo -e "\n\n############# START ALIGN BLOCK ###############\n\n" >&2
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
nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 )
scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
aligner=$( cat $runfile | grep -w ALIGNER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
epilogue=$( cat $runfile | grep -w EPILOGUE  | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

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
	echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
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
fi

if [ $bamtofastqflag == "1" ]
then
    bamtofastqflag="YES"
fi
if [ $bamtofastqflag == "0" ]
then
    bamtofastqflag="NO"
fi


if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
then
    MSG="$outputdir/SAMPLENAMES_multiplexed.list configuration file not found"
    echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ $aligner != "NOVOALIGN" -a $aligner != "BWA_ALN" -a $aligner != "BWA_MEM" ]
then
    MSG="ALIGNER=$aligner  processing this case is not available yet"
    echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

#if [ -z $epilogue ]
#then
#   MSG="Value for EPILOGUE must be specified in configuration file"
#   echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
#   exit 1;
#else
#   `chmod 740 $epilogue`
#fi

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

TopLogsFolder=$outputdir/logs
pipeid=$( cat $TopLogsFolder/pbs.CONFIGURE )
pbsids=""


set +x; echo -e "\n\n#####################    Do we need to run Preprocessing steps? if FASTQ, then NO, else YES     #####\n\n" >&2; set -x; 

#ALL input files that will be processed with this pipeline MUST have the same input format


if [ $inputformat == "FASTQ" ]
then
	set +x; echo -e "\n\n" >&2
	echo "####################################################################" >&2
	echo "####################################################################" >&2
	echo "############SKIPPING preprocessing because INPUT is in FASTQ #######" >&2
	echo "####################################################################" >&2
	echo "####################################################################" >&2
	echo -e "\n\n" >&2; set -x
	
	qsub3=$TopLogsFolder/qsub.alignfastq

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

else
	set +x; echo -e "\n\n" >&2
	echo "####################################################################" >&2
	echo "####################################################################" >&2
	echo "############    preprocessing INPUT files before alignment   #######" >&2
	echo "####################################################################" >&2
	echo "####################################################################" >&2
	echo -e "\n\n" >&2; set -x

	CONVERT=""
	case=""
	typeOfupdateconfig=""
	truncate -s 0 $TopLogsFolder/pbs.CONVERTBAM


	set +x; echo -e "\n #input files are BAMs; preprocessing is required before aligning files \n" >&2; set -x;
	newfqfiles=""
	newbamfiles=""
	sep=":"


	while read SampleName
	do
              # this will evaluate the length of string; 
              if [ `expr ${#SampleName}` -lt 1 ]
              then
                    set +x; echo -e "\n\n #### skipping empty line #### \n\n" >&2; set -x;              

              else
              
                    set +x; echo -e "\n\n #### processing $SampleName #### \n\n" >&2; set -x;
                    set +x; echo -e "\n\n #### parsing the line. It should have two fields: samplename inputfile.bam #### \n\n" >&2; set -x;
                    
                    sample=$( echo $SampleName | cut -d ' ' -f 1 )
                    inputbam=$( echo $SampleName | cut -d ' ' -f 2 )
                    
                    if [ ! -s $inputbam ]
		    then
			MSG="$inputbam input BAM file not found. Parsing of $SampleName failed"
			echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
			exit 1;
		    fi                    


		    AlignmentOutputFolder=$outputdir/$sample/align
		    AlignmentOutputLog=$outputdir/$sample/align/logs
		    AlignedOutputFile=${sample}.wdups.sorted.bam
		    
		    if [ -d $AlignmentOutputFolder ]
		    then
			echo "$AlignmentOutputFolder is there; resetting it"
			`rm -r $AlignmentOutputFolder/*`
			mkdir -p $AlignmentOutputLog
		    else
		        mkdir -p $AlignmentOutputLog
		    fi


                    set +x; echo -e "\n\n #### TWO CASES: convert/not convert input.bam to fastq #### \n\n" >&2; set -x;
                    
		    if [ $bamtofastqflag == "NO" ] 
                    then
			set +x; echo -e "\n # bam2newbam preprocessing... \n" >&2; set -x;
			typeOfupdateconfig="bam2newbam"

			#########################
			# TODO: new script goes here
			# this script performs: namesort, remove singletons, revertsam (if so indicated)
			# input is originalBAM, output is newsuffix4bam
			#########################

			qsub1=$AlignmentOutputLog/qsub.convertbam2newbam.$sample
			echo "#PBS -A $pbsprj" >> $qsub1
			echo "#PBS -N ${pipeid}_convertbam2newbam_${sample}" >> $qsub1
			echo "#pbs -l epilogue=$epilogue" >> $qsub1
			echo "#PBS -l walltime=$pbscpu" >> $qsub1
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
			echo "#PBS -o $AlignmentOutputLog/log.convertbam2newbam.${sample}.ou" >> $qsub1
			echo "#PBS -e $AlignmentOutputLog/log.convertbam2newbam.${sample}.in" >> $qsub1
			echo "#PBS -q $pbsqueue" >> $qsub1
			echo "#PBS -m ae" >> $qsub1
			echo "#PBS -M $email" >> $qsub1
			echo "$scriptdir/convertbam2newbam.sh $AlignmentOutputFolder $inputbam $AlignedOutputFile $runfile $AlignmentOutputLog/log.convertbam2newbam.${sample}.in $AlignmentOutputLog/log.convertbam2newbam.${sample}.ou $email $AlignmentOutputLog/qsub.convertbam2newbam.$sample" >> $qsub1
			echo -e "\n\n" >> $qsub1
			echo "exitcode=\$?" >> $qsub1
			echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub1
			echo "   echo -e \"\n\n convertbam2newbam.sh failed with exit code = \$exitcode \n logfile=$AlignmentOutputLog/log.convertbam2newbam.${sample}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub1
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

			qsub1=$AlignmentOutputLog/qsub.convertbam2fq.$sample

			echo "#PBS -A $pbsprj" >> $qsub1
			echo "#PBS -N ${pipeid}_convertbam2fq_${sample}" >> $qsub1
			echo "#pbs -l epilogue=$epilogue" >> $qsub1
			echo "#PBS -l walltime=$pbscpu" >> $qsub1
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
			echo "#PBS -o $AlignmentOutputLog/log.convertbam2fq.${sample}.ou" >> $qsub1
			echo "#PBS -e $AlignmentOutputLog/log.convertbam2fq.${sample}.in" >> $qsub1
			echo "#PBS -q $pbsqueue" >> $qsub1
			echo "#PBS -m ae" >> $qsub1
			echo "#PBS -M $email" >> $qsub1
			echo "$scriptdir/convertbam2fastq.sh $AlignmentOutputFolder $inputbam $sample $runfile $AlignmentOutputLog/log.convertbam2fq.${sample}.in $AlignmentOutputLog/log.convertbam2fq.${sample}.ou $email $AlignmentOutputLog/qsub.convertbam2fq.$sample" >> $qsub1

                        echo "exitcode=\$?" >> $qsub1
                        echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub1
                        echo "   echo -e \"\n\n convertbam2fastq.sh failed with exit code = \$exitcode \n logfile=$AlignmentOutputLog/log.convertbam2fq.${sample}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub1
                        echo "fi" >> $qsub1
                        echo -e "\n\n exit 1" >> $qsub1

			#`chmod a+r $qsub1` 
			convertjob=`qsub $qsub1`
			`qhold -h u $convertjob` 
			echo $convertjob >> $TopLogsFolder/pbs.CONVERTBAM
			echo `date`
                    fi  #end if $bamtofastqflag
		
              fi  # end non-empty line
              
	done < $outputdir/SAMPLENAMES_multiplexed.list
	
	set +x; echo -e "\n\n #### End loop over samples to preprocess input BAMs! #### \n\n" >&2; set -x; 
	set +x; echo -e "\n\n #### WE NEED TO UPDATE THE CONFIGURATION FILES #### \n\n" >&2; set -x; 

	CONVERTids=$( cat $TopLogsFolder/pbs.CONVERTBAM | sed "s/\..*//" | tr "\n" ":" )

	if [ $bamtofastqflag == "YES" ]
	then
                ###############################
                # todo:
                # rename updateconfig script for newfastq files
                # rename updateconfig.sh to updateconfig.wnewfq.sh
                ###############################

		qsub2=$TopLogsFolder/qsub.updateconfig_wnewfq

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
		echo "$scriptdir/updateconfig.wnewfq.sh $sampledir $newfqfiles $runfile $samplefileinfo $TopLogsFolder/log.updateconfig_wnewfq.in $TopLogsFolder/log.updateconfig_wnewfq.ou $email $TopLogsFolder/qsub.updateconfig_wnewfq" >> $qsub2

                echo "exitcode=\$?" >> $qsub2
                echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsub2
                echo "   echo -e \"\n\n updateconfig.wnewfq.sh failed with exit code = \$exitcode \n logfile=$TopLogsFolder/log.updateconfig_wnewfq.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsub2
                echo "fi" >> $qsub2
                echo -e "\n\n exit 1" >> $qsub2

		#`chmod a+r $qsub2`       
		updatejob=`qsub $qsub2` 
		echo $updatejob >> $TopLogsFolder/pbs.UPDATECONFIG
	fi

	allconjobs=$( echo $CONVERTids | tr ":" " " )
	`qrls -h u $allconjobs`
	echo `date`


        set +x; echo -e "\n\n" >&2
        echo "####################################################################" >&2
        echo "############    done with update config             ################" >&2
        echo "############    alignment case selection is next    ################" >&2
        echo "####################################################################" >&2
        echo -e "\n\n" >&2; set -x;

        if [ $bamtofastqflag == "YES" -a $inputformat == "BAM" ]
        then
            set +x; echo -e "\n# input is BAM, convert to fastq\n" >&2;  set -x;
	    qsub1=$TopLogsFolder/qsub.main.alnFQ.afterbam2fastq

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


        if [ `expr ${#case}` -lt 1 ]
        then
           MSG="Alignment module failed to launch. Incompatible values specified in config files bam2fastqflag=$bamtofastqflag inputformat=$inputformat analysis=$analysis"
           echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
           exit 1;
        fi         

fi # end cases with inputformat

echo -e "now we need to make the PBS log files group readable"

find $outputdir -name logs -type d | awk '{print "chmod -R g=rwx "$1}' | sh -x
