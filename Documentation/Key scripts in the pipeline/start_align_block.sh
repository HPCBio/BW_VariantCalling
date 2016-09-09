#!/bin/bash
#
# start_align_block.sh
# First module in the GGPS analysis pipeline. This module launches qsubs to convert file formats -if needed- and to align data according to input type
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

set +x; echo -e "\n\n############# CHECKING PARAMETERS ###############\n\n"; set -x;

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

set +x; echo -e "\n\n\n############ checking input type: WGS or WES\n" >&2; set -x

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

set +x; echo -e "\n\n\n############ checking input format\n" >&2; set -x


if [ $inputformat != "FASTQ" -a $inputformat != "BAM" ]
then
    MSG="Incorrect value for INPUTFORMAT=$inputformat in the configuration file."
    echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

set +x; echo -e "\n\n\n############ checking input conversion flag\n" >&2; set -x


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

set +x; echo -e "\n\n\n############ checking sample configuration file\n" >&2; set -x

if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
then
    MSG="$outputdir/SAMPLENAMES_multiplexed.list configuration file not found"
    echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

set +x; echo -e "\n\n\n############ checking aligner\n" >&2; set -x

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

set +x; echo -e "\n\n\n############ checking workflow scripts directory\n" >&2; set -x;

if [ ! -d $scriptdir ]
then
   MSG="SCRIPTDIR=$scriptdir directory not found"
   echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi
set +x; echo -e "\n\n\n############ checking output directory\n" >&2; set -x;

if [ ! -d $outputdir ]
then
   mkdir -p $outputdir
fi


set +x; echo -e "\n\n">&2
echo "############################################################################################################" >&2
echo "#############  CREATE  DIRECTORY STRUCTURE FOR LOGS AND RESET FILES. TOP LEVEL ONLY       ##################" >&2
echo "############################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

TopLogsFolder=$outputdir/logs
pipeid=$( cat $TopLogsFolder/pbs.CONFIGURE )
pbsids=""


set +x; echo -e "\n\n">&2
echo "############################################################################################################" >&2
echo "#############  QSUBS FOR SCHEDULERS OF ALIGNMENT JOBS                                     ##################" >&2
echo "#############  CASE1: input format is FASTQ                                               ##################" >&2
echo "#############  CASE2: input format is BAM and conversion from bam to fastq IS NOT NEEDED  ##################" >&2
echo "#############  CASE3: input format is BAM and conversion from bam to fastq IS NEEDED      ##################" >&2
echo "############################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


if [ $inputformat == "FASTQ" ]
then
	set +x; echo -e "\n\n" >&2
	echo "############################################################################################################" >&2
	echo "############  CASE1: input format is FASTQ " >&2
	echo "############################################################################################################" >&2
	echo -e "\n\n" >&2; set -x
	
	qsubSchAlignFq=$TopLogsFolder/qsub.alignfastq

	echo "#PBS -A $pbsprj" >> $qsubSchAlignFq
	echo "#PBS -N ${pipeid}_alignfastq" >> $qsubSchAlignFq
	echo "#PBS -l epilogue=$epilogue" >> $qsubSchAlignFq
	echo "#PBS -l walltime=$pbscpu" >> $qsubSchAlignFq
	echo "#PBS -l nodes=1:ppn=1" >> $qsubSchAlignFq
	echo "#PBS -o $TopLogsFolder/log.alignfastq.ou" >> $qsubSchAlignFq
	echo "#PBS -e $TopLogsFolder/log.alignfastq.in" >> $qsubSchAlignFq
	echo "#PBS -q $pbsqueue" >> $qsubSchAlignFq
	echo "#PBS -m ae" >> $qsubSchAlignFq
	echo "#PBS -M $email" >> $qsubSchAlignFq
	echo "$scriptdir/schedule_alignfastq.sh $runfile $TopLogsFolder/log.alignfastq.in $TopLogsFolder/log.alignfastq.ou $email $qsubSchAlignFq" >> $qsubSchAlignFq
	#`chmod a+r $qsubSchAlignFq` 
	`qsub $qsubSchAlignFq >> $TopLogsFolder/pbs.ALIGN`
	case="alignfastq"
	echo `date`

elif [ $bamtofastqflag == "NO" -a $inputformat == "BAM" ]
then
	set +x; echo -e "\n\n" >&2
	echo "############################################################################################################" >&2
	echo "############  CASE2: input format is BAM and conversion from bam to fastq IS NOT NEEDED" >&2
	echo "############################################################################################################" >&2
	echo -e "\n\n" >&2; set -x

	##############################################################
	#this section needs editing
	##############################################################


	qsubSchAlignBam=$TopLogsFolder/qsub.alignbams

	echo "#PBS -A $pbsprj" >> $qsubSchAlignBam
	echo "#PBS -N ${pipeid}_alignbam" >> $qsubSchAlignBam
	echo "#PBS -l epilogue=$epilogue" >> $qsubSchAlignBam
	echo "#PBS -l walltime=$pbscpu" >> $qsubSchAlignBam
	echo "#PBS -l nodes=1:ppn=1" >> $qsubSchAlignBam
	echo "#PBS -o $TopLogsFolder/log.alignbams.ou" >> $qsubSchAlignBam
	echo "#PBS -e $TopLogsFolder/log.alignbams.in" >> $qsubSchAlignBam
	echo "#PBS -q $pbsqueue" >> $qsubSchAlignBam
	echo "#PBS -m ae" >> $qsubSchAlignBam
	echo "#PBS -M $email" >> $qsubSchAlignBam
	echo "#PBS -W depend=afterok:$updatejob" >> $qsubSchAlignBam
	echo "$scriptdir/schedule_alignbam.sh $runfile $TopLogsFolder/log.alignbams.in $TopLogsFolder/log.alignbams.ou $email $qsubSchAlignBam" >> $qsubSchAlignBam
	#`chmod a+r $qsubSchAlignBam`               
	`qsub $qsubSchAlignBam >> $TopLogsFolder/pbs.ALIGN`
	case="alignbams"
	echo `date`

elif [ $bamtofastqflag == "YES" -a $inputformat == "BAM" ]
then
	set +x; echo -e "\n\n" >&2
	echo "############################################################################################################" >&2
	echo "############  CASE3: input format is BAM and conversion from bam to fastq IS NEEDED" >&2
	echo "############################################################################################################" >&2
	echo -e "\n\n" >&2; set -x

	##############################################################
	#this section needs editing
	##############################################################
                

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

			set +x; echo -e "\n\n #### define output folder and output file  #### \n\n" >&2; set -x;

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
                    
			echo "bam2fastq preprocessing..."
			typeOfupdateconfig="bam2fastq"

			#########################                        
			# TODO: add preprocessing stuff to convertbam.sh
			# the new convertbam.sh  script must perform: 
			# namesort, remove singletons, revertsam (if so indicated, bam2fastq
			# input is originalBAM, output is newsuffix4fq
			#########################

			qsubConvert=$AlignmentOutputLog/qsub.convertbam2fq.$sample

			echo "#PBS -A $pbsprj" >> $qsubConvert
			echo "#PBS -N ${pipeid}_convertbam2fq_${sample}" >> $qsubConvert
			echo "#pbs -l epilogue=$epilogue" >> $qsubConvert
			echo "#PBS -l walltime=$pbscpu" >> $qsubConvert
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsubConvert
			echo "#PBS -o $AlignmentOutputLog/log.convertbam2fq.${sample}.ou" >> $qsubConvert
			echo "#PBS -e $AlignmentOutputLog/log.convertbam2fq.${sample}.in" >> $qsubConvert
			echo "#PBS -q $pbsqueue" >> $qsubConvert
			echo "#PBS -m ae" >> $qsubConvert
			echo "#PBS -M $email" >> $qsubConvert
			echo "$scriptdir/convertbam2fastq.sh $AlignmentOutputFolder $inputbam $sample $runfile $AlignmentOutputLog/log.convertbam2fq.${sample}.in $AlignmentOutputLog/log.convertbam2fq.${sample}.ou $email $qsubConvert" >> $qsubConvert

			echo "exitcode=\$?" >> $qsubConvert
			echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubConvert
			echo "   echo -e \"\n\n convertbam2fastq.sh failed with exit code = \$exitcode \n logfile=$AlignmentOutputLog/log.convertbam2fq.${sample}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubConvert
			echo "fi" >> $qsubConvert
			echo -e "\n\n exit 1" >> $qsubConvert

			#`chmod a+r $qsubConvert` 
			convertjob=`qsub $qsubConvert`
			`qhold -h u $convertjob` 
			echo $convertjob >> $TopLogsFolder/pbs.CONVERTBAM
			echo `date`

		
		fi  # end non-empty line
              
	done < $outputdir/SAMPLENAMES_multiplexed.list
	
	set +x; echo -e "\n\n #### End loop over samples to preprocess input BAMs! #### \n\n" >&2; set -x; 
	set +x; echo -e "\n\n #### WE NEED TO UPDATE THE CONFIGURATION FILES #### \n\n" >&2; set -x; 

	CONVERTids=$( cat $TopLogsFolder/pbs.CONVERTBAM | sed "s/\..*//" | tr "\n" ":" )


	qsubUpdate=$TopLogsFolder/qsub.updateconfig_wnewfq

	echo "#PBS -A $pbsprj" >> $qsubUpdate
	echo "#PBS -N ${pipeid}_updateconfig_wnewfq" >> $qsubUpdate
	echo "#pbs -l epilogue=$epilogue" >> $qsubUpdate
	echo "#PBS -l walltime=$pbscpu" >> $qsubUpdate
	echo "#PBS -l nodes=1:ppn=1" >> $qsubUpdate
	echo "#PBS -o $TopLogsFolder/log.updateconfig_wnewfq.ou" >> $qsubUpdate
	echo "#PBS -e $TopLogsFolder/log.updateconfig_wnewfq.in" >> $qsubUpdate
	echo "#PBS -q $pbsqueue" >> $qsubUpdate
	echo "#PBS -m ae" >> $qsubUpdate
	echo "#PBS -M $email" >> $qsubUpdate
	echo "#PBS -W depend=afterok:$CONVERTids" >> $qsubUpdate
	echo "$scriptdir/updateconfig.wnewfq.sh $sampledir $newfqfiles $runfile $samplefileinfo $TopLogsFolder/log.updateconfig_wnewfq.in $TopLogsFolder/log.updateconfig_wnewfq.ou $email $qsubUpdate" >> $qsubUpdate

	echo "exitcode=\$?" >> $qsubUpdate
	echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubUpdate
	echo "   echo -e \"\n\n updateconfig.wnewfq.sh failed with exit code = \$exitcode \n logfile=$TopLogsFolder/log.updateconfig_wnewfq.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubUpdate
	echo "fi" >> $qsubUpdate
	echo -e "\n\n exit 1" >> $qsubUpdate

	#`chmod a+r $qsubUpdate`       
	updatejob=`qsub $qsubUpdate` 
	echo $updatejob >> $TopLogsFolder/pbs.UPDATECONFIG

	allconjobs=$( echo $CONVERTids | tr ":" " " )
	`qrls -h u $allconjobs`
	echo `date`


	set +x; echo -e "\n# input is BAM, convert to fastq\n" >&2;  set -x;
	qsubAlignFq=$TopLogsFolder/qsub.alignfastq.afterbam2fastq

	echo "#PBS -A $pbsprj" >> $qsubAlignFq
	echo "#PBS -N ${pipeid}_alnFQ_afterbam2fastq" >> $qsubAlignFq
	echo "#pbs -l epilogue=$epilogue" >> $qsubAlignFq
	echo "#PBS -l walltime=$pbscpu" >> $qsubAlignFq
	echo "#PBS -l nodes=1:ppn=1" >> $qsubAlignFq
	echo "#PBS -o $TopLogsFolder/log.alnFQ.afterbam2fastq.ou" >> $qsubAlignFq
	echo "#PBS -e $TopLogsFolder/log.alnFQ.afterbam2fastq.in" >> $qsubAlignFq
	echo "#PBS -q $pbsqueue" >> $qsubAlignFq
	echo "#PBS -m ae" >> $qsubAlignFq
	echo "#PBS -M $email" >> $qsubAlignFq
	echo "#PBS -W depend=afterok:$updatejob" >> $qsubAlignFq
	echo "$scriptdir/schedule_alignfastq.sh $runfile $TopLogsFolder/log.alnFQ.afterbam2fastq.in $TopLogsFolder/log.alnFQ.afterbam2fastq.ou $email $qsubAlignFq" >> $qsubAlignFq
	
	#`chmod a+r $qsubAlignFq`               
	`qsub $qsubAlignFq >> $TopLogsFolder/pbs.ALIGN`
	case="bam2fastq"
	echo `date`


else
	MSG="Alignment module failed to launch. Incompatible values specified in config files bam2fastqflag=$bamtofastqflag inputformat=$inputformat analysis=$analysis"
	echo -e "Program $scriptfile stopped at line=$LINENO.\n$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
    
fi # end cases with inputformat

set +x; echo -e "now we need to make the PBS log files group readable" >&2; set -x; 

find $outputdir -name logs -type d | awk '{print "chmod -R g=rwx "$1}' | sh -x

echo `date`

set +x; echo -e "\n\n">&2
echo "############################################################################################################" >&2
echo "#############  DONE. EXITING NOW                                                         ##################" >&2
echo "############################################################################################################" >&2
echo -e "\n\n" >&2; set -x;
