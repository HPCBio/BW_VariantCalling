#!/bin/bash
########################### 
# program configure.sh is the script that generates the sample configuration files
# and selects the types of analyses that will be executed in this run
###########################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch."
        echo -e "Program $0 stopped. Reason=$MSG" # | mail -s "Variant Calling Workflow failure message" "$redmine"
        exit 1;
fi


echo -e "\n\n############# CONFIGURE VARIANT CALLING WORKFLOW ###############\n\n" >&2

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



set +x; echo -e "\n\n############# CHECKING PARAMETERS ###############\n\n" >&2; set -x;

if [ !  -s $runfile ]
then
   MSG="$runfile configuration file not found."
   echo -e "Program $0 stopped. Reason=$MSG" # | mail -s "Variant Calling Workflow failure message" "$redmine"
   exit 1;
fi


set +x
echo -e "\n\n####################################################################################################" >&2
echo        "##################################### PARSING RUN INFO FILE ########################################" >&2
echo        "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "\n\n####################################################################################################\n\n" >&2; set -x;

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
resortbam=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
bamtofastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
#sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
sampleinfo=$outputdir/sampleinfo.txt

set +x; echo -e "\n\n\n############ checking input type: WGS or WES\n" >&2; set -x

if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
	pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
	pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
elif [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
	pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
	pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
else
	MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n\n############ checking output directory\n" >&2; set -x;

if [ -z $thr -o -z $outputdir ]
then
	MSG="Invalid value specified for any of these paramaters in configuration file:\nPBSTHREADS=$thr\nOUTPUTDIR=$outputdir\nPBSPROJECTID=$pbsprj\nEPILOGUE=$epilogue"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n\n############ checking analysis\n" >&2; set -x;

if [ `expr ${#analysis}` -lt 1 ]
then
	MSG="Incorrect value for ANALYSIS=$analysis in the configuration file."
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi


set +x; echo -e "\n\n\n############ checking input format\n" >&2; set -x;

if [ $inputformat != "FASTQ" -a $inputformat != "BAM" ]
then
	MSG="Incorrect value for INPUTFORMAT=$inputformat in the configuration file."
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n\n############ checking resortbam option\n" >&2; set -x;


if [ $resortbam != "1" -a $resortbam != "0" -a $resortbam != "YES" -a $resortbam != "NO" ]
then
	MSG="Invalid value for RESORTBAM=$resortbam"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi
if [ $resortbam == "1" ]
then
	resortbam="YES"
fi
if [ $resortbam == "0" ]
then
	resortbam="NO"
fi

set +x; echo -e "\n\n\n############ checking bamtofastq option\n" >&2; set -x;

if [ $bamtofastqflag != "YES" -a $bamtofastqflag != "NO" -a $bamtofastqflag != "1" -a $bamtofastqflag != "0" ]
then
	MSG="BM2FASTQFLAG=$bamtofastqflag  invalid value"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
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

set +x; echo -e "\n\n\n############ checking multisamples option\n" >&2; set -x;

if [ $multisample != "YES" -a $multisample != "NO" -a $multisample != "1" -a $multisample != "0" ]
then
	MSG="MULTISAMPLE=$multisample  invalid value"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ $multisample == "1" ]
then
	multisample="YES"
fi
if [ $multisample == "0" ]
then
	multisample="NO"
fi

set +x; echo -e "\n\n\n############ checking clash in conversion options\n" >&2; set -x;

if [ $resortbam == "YES" -a $bamtofastqflag == "YES" ]
then
	MSG="Incompatible values for the pair RESORTBAM=$resortbam and BAM2FASTQFLAG=$bam2fqflag in the configuration file."
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\nChecking  samples comfiguration file\n" >&2; set -x;

if [ ! -s $sampleinfo ]
then
	MSG="SAMPLEINFORMATION=$sampleinfo invalid value. A tab delimited file with lanes and samples must be specified. "
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n\n############ checking workflow scripts directory\n" >&2; set -x;

if [ ! -d $scriptdir ]
then
	MSG="SCRIPTDIR=$scriptdir directory not found"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
#`chmod -R 770 $outputdir`
#`chmod 740 $epilogue`

set +x; echo -e "\n\n\n############ checking output directory\n" >&2; set -x;

TopOutputLogs=$outputdir/logs
pipeid=$( cat $TopOutputLogs/pbs.CONFIGURE )

if [ ! -d $outputdir ]
then
	set +x; echo -e "\n new output folder creation\n" >&2; set -x; 
	mkdir -p $outputdir
	mkdir -p $outputdir/logs
else
	set +x; echo -e "\n resetting output folder's logs\n" >&2; set -x; 
	#`rm -r $outputdir/logs/*`
	`rm $TopOutputLogs/pbs.ALIGNED`
	`rm $TopOutputLogs/pbs.CONVERTBAM`
	`rm $TopOutputLogs/pbs.FASTQC`
	`rm $TopOutputLogs/pbs.MARKED`
	`rm $TopOutputLogs/pbs.REALIGN`
	`rm $TopOutputLogs/pbs.RECALIBRATE`
	
fi



set +x; echo -e "\n\n">&2
echo "#################################################################################################################################################################################" >&2
echo "############# constructing files with list(s) of input files to analyze in this run of the pipeline" >&2
echo "############# CASE1: analysis is multiplexed. Input format is FASTQ. Info sheet is parsed and three files are generated. SAMPLENAMES_multiplexed.list has five columns" >&2
echo "############# CASE2: analysis is NOT multiplexed. Input format is FASTQ. Info sheet is parsed and three files are generated. SAMPLENAMES_multiplexed.list has three columns" >&2
echo "############# CASE3: analysis is NOT multiplexed. Input format is BAM. Info sheet is parsed and three files are generated. SAMPLENAMES_multiplexed.list has two columns" >&2 
echo "#################################################################################################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


set +x; echo -e "\n ### run perl script to create three configuration files ### \n"  >&2; set -x;
cd $outputdir
sampleinfo=`basename $sampleinfo`

### all options below call the script configSAMPLENAMES.pl with the same parameters 
### we use this nested if-then-else to document the cases

if [ $analysis == "MULTIPLEXED" -a $inputformat == "FASTQ" ]
then
	set +x; echo -e "\n\n ######## CASE1: ANALYSIS IS MULTIPLEXED. Input format is FASTQ ! ############\n\n" >&2; 
	echo -e " ###### produce SAMPLENAMES.list SAMPLENAMES_multiplexed.list SAMPLEGROUPS.list " >&2
	echo -e " ###### SAMPLENAMES_multiplexed.list has five columns: sampleid r1 r2 flowcell lib  " >&2            
	echo -e "\n ###### from  info sheet specified in runfile in line SAMPLEINFORMATION=$sampleinfo." >&2; set -x;

	perl $scriptdir/configSAMPLENAMES.pl $outputdir $sampleinfo $inputformat SAMPLENAMES.list SAMPLENAMES_multiplexed.list SAMPLEGROUPS.list

elif [ $analysis != "MULTIPLEXED" -a $inputformat == "FASTQ" ]
then
	echo +x; echo -e "\n ###### CASE2: ANALYSIS IS NOT MULTIPLEXED. Input format is FASTQ" >&2;
	echo -e "###### Parse info sheet in $sampleinfo  and produce three files as we did in the multiplexed case" >&2;
	echo -e "###### SAMPLENAMES_multiplexed.list has three columns: sampleid r1 r2  " >&2 	     	
	echo -e "###### At least one of them will be redundant SAMPLEGROUPS.list. No big deal, as long as we don't break the code" >&2; set -x;
	
	perl $scriptdir/configSAMPLENAMES.pl $outputdir $sampleinfo $inputformat SAMPLENAMES.list SAMPLENAMES_multiplexed.list SAMPLEGROUPS.list

elif [ $analysis != "MULTIPLEXED" -a $inputformat == "BAM" ]
then
	echo +x; echo -e "\n ###### CASE3: ANALYSIS IS NOT MULTIPLEXED AND Input format is BAM." >&2;
	echo -e "###### Parse info sheet in $sampleinfo  and produce three files as we did in the multiplexed case" >&2;
	echo -e " ###### SAMPLENAMES_multiplexed.list has two columns: sampleid bamfile  " >&2 	     	
	echo -e "###### At least one of them will be redundant SAMPLEGROUPS.list. No big deal, as long as we don't break the code" >&2; set -x;

	perl $scriptdir/configSAMPLENAMES.pl $outputdir $sampleinfo $inputformat SAMPLENAMES.list SAMPLENAMES_multiplexed.list SAMPLEGROUPS.list

else
	MSG="Configuration case is NOT covered. ANALYSIS=$analysis INPUTFORMAT=$inputformat SAMPLEINFORMATION=$sampleinfo"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;

fi

set +x; echo -e "\n ### testing that the perl script actually worked. it should always create three files ### \n"  >&2; set -x;

if [ ! -s $outputdir/SAMPLENAMES_multiplexed.list ]
then
	MSG="$outputdir/SAMPLENAMES_multiplexed.list is empty"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -s $outputdir/SAMPLENAMES.list ]
then
	MSG="$outputdir/SAMPLENAMES.list is empty"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -s $outputdir/SAMPLEGROUPS.list ]
then
	MSG="$outputdir/SAMPLEGROUPS.list is empty"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n ### update autodocumentation script ### \n"; set -x;
numinputs=`wc -l $outputdir/SAMPLENAMES_multiplexed.list | cut -d ' ' -f 1`
numsamplegroups=`wc -l $outputdir/SAMPLEGROUPS.list | cut -d ' ' -f 1`
echo -e "# @in fastq_inputs @URI SAMPLENAMES_multiplexed.list=${numinputs}_inputs_for_${numsamplegroups}_samples" >> $outputdir/WorkflowAutodocumentationScript.sh
    

set +x; echo -e "\n\n">&2
echo "###########################################################################################################################################################" >&2
echo "#############Two more configuration tasks: generate a QC_Results file and a qsub header so we would not have to repeat the same lines\n\n" >&2; set -x;
echo "###########################################################################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

generic_qsub_header=$outputdir/qsubGenericHeader
truncate -s 0 $outputdir/QC_Results.txt  # to report results of all QC tests
truncate -s 0 $generic_qsub_header       # to print the generic qsub header lines

echo "#!/bin/bash" > $generic_qsub_header
echo "#PBS -A $pbsprj" >> $generic_qsub_header
echo "#PBS -q $pbsqueue" >> $generic_qsub_header
echo "#PBS -m ae" >> $generic_qsub_header
echo "#PBS -M $email" >> $generic_qsub_header
set +x; echo -e "\n### check that this actually worked, " >&2
echo -e "### because otherwise the bash script will just go on, as if there is no problem \n"  >&2; set -x;
if [ ! -s $generic_qsub_header ]
then 
	MSG="$generic_qsub_header is empty"
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi



set +x; echo -e "\n\n">&2
echo "###########################################################################################################################################################" >&2
echo "#############  select analysis or analyses to run. possible values for ANALYSIS are "  >&2 
echo "#############  vcall_only realign_only align_only align_and_realign align_and_realign_and_vcall"  >&2 
echo "###########################################################################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


case=""
if [ $analysis == "REALIGNONLY" -o $analysis == "REALIGN_ONLY" ]
then
	echo "############# Type of analysis to run: REALIGNMENT only. bams must be provided as input############# "
	qsub2=$TopOutputLogs/qsub.START_REALRECAL_BLOCK
	echo "#!/bin/bash" > $qsub2
	echo "#PBS -A $pbsprj" >> $qsub2
	echo "#PBS -N ${pipeid}_START_REALRECAL_BLOCK"
	echo "#pbs -l epilogue=$epilogue" >> $qsub2
	echo "#PBS -l walltime=00:30:00" >> $qsub2
	echo "#PBS -l nodes=1:ppn=1" >> $qsub2
	echo "#PBS -o $TopOutputLogs/log.START_REALRECAL_BLOCK.ou" >> $qsub2
	echo "#PBS -e $TopOutputLogs/log.START_REALRECAL_BLOCK.in" >> $qsub2
	echo "#PBS -q $pbsqueue" >> $qsub2
	echo "#PBS -m ae" >> $qsub2
	echo "#PBS -M $email" >> $qsub2
	echo "$scriptdir/start_realrecal_block.sh $runfile $TopOutputLogs/log.START_REALRECAL_BLOCK.in $TopOutputLogs/log.START_REALRECAL_BLOCK.ou $email $TopOutputLogs/qsub.START_REALRECAL_BLOCK" >> $qsub2
	#`chmod a+r $qsub2` 
	`qsub $qsub2 >> $TopOutputLogs/pbs.REALRECAL`
	echo `date`
	case="realign_only" 
 fi
 
 if [ $analysis == "VCALL_ONLY" -o $analysis == "VCALL" ]
 then
 	echo "############# Type of analysis to run: VARIANT CALLING ONLY. cleaned bams must be provided as input############# "
 	qsub3=$TopOutputLogs/qsub.START_VARCALL_BLOCK
 	echo "#!/bin/bash" > $qsub3
 	echo "#PBS -A $pbsprj" >> $qsub3
 	echo "#PBS -N ${pipeid}_START_VARCALL_BLOCK" >> $qsub3
 	echo "#PBS -l epilogue=$epilogue" >> $qsub3
 	echo "#PBS -l walltime=00:30:00" >> $qsub3
 	echo "#PBS -l nodes=1:ppn=1" >> $qsub3
 	echo "#PBS -o $TopOutputLogs/log.START_VARCALL_BLOCK.ou" >> $qsub3
 	echo "#PBS -e $TopOutputLogs/log.START_VARCALL_BLOCK.in" >> $qsub3
 	echo "#PBS -q $pbsqueue" >> $qsub3
 	echo "#PBS -m ae" >> $qsub3
 	echo "#PBS -M $email" >> $qsub3
 	echo "$scriptdir/start_varcall_block.sh $runfile $TopOutputLogs/log.START_VARCALL_BLOCK.in $TopOutputLogs/log.START_VARCALL_BLOCK.ou $email $TopOutputLogs/qsub.START_VARCALL_BLOCK" >> $qsub3
 	#`chmod a+r $qsub3`
 	vcalljobid=`qsub $qsub3`
 	echo $vcalljobid >> $TopOutputLogs/pbs.VARCALL
 	case="vcall_only"  
 fi

 if [ $analysis == "ALIGN" -o $analysis == "ALIGNMENT" -o $analysis == "ALIGNONLY" ]
 then
	echo "############# Type of analysis to run: ALIGNMENT only. fastq or unaligned bams must be provided as input"      
	qsub1=$TopOutputLogs/qsub.START_ALIGN_BLOCK
	echo "#!/bin/bash" > $qsub1
	echo "#PBS -A $pbsprj" >> $qsub1
	echo "#PBS -N ${pipeid}_START_ALIGN_BLOCK" >> $qsub1
	echo "#pbs -l epilogue=$epilogue" >> $qsub1
	echo "#PBS -l walltime=00:30:00" >> $qsub1
	echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	echo "#PBS -o $TopOutputLogs/log.START_ALIGN_BLOCK.ou" >> $qsub1
	echo "#PBS -e $TopOutputLogs/log.START_ALIGN_BLOCK.in" >> $qsub1
	echo "#PBS -q $pbsqueue" >> $qsub1
	echo "#PBS -m ae" >> $qsub1
	echo "#PBS -M $email" >> $qsub1
	echo "$scriptdir/start_align_block.sh $runfile $TopOutputLogs/log.START_ALIGN_BLOCK.in $TopOutputLogs/log.START_ALIGN_BLOCK.ou $email $TopOutputLogs/qsub.START_ALIGN_BLOCK" >> $qsub1
	#`chmod a+r $qsub1`               
	`qsub $qsub1 >> $TopOutputLogs/pbs.ALIGN`
	echo `date`
	case="align_only"
 fi

 if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" -o $analysis == "MULTIPLEXED" ]
 then
	echo "Type of analysis to run: ALIGNMENT and REALIGNMENT and/or VARiANT CALLING. fastq or unaligned bams must be provided as input"
	qsub1=$TopOutputLogs/qsub.START_ALIGN_BLOCK
	echo "#!/bin/bash" > $qsub1
	echo "#PBS -A $pbsprj" >> $qsub1
	echo "#PBS -N ${pipeid}_START_ALIGN_BLOCK" >> $qsub1
	echo "#pbs -l epilogue=$epilogue" >> $qsub1
	echo "#PBS -l walltime=00:30:00" >> $qsub1
	echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	echo "#PBS -o $TopOutputLogs/log.START_ALIGN_BLOCK.ou" >> $qsub1
	echo "#PBS -e $TopOutputLogs/log.START_ALIGN_BLOCK.in" >> $qsub1
	echo "#PBS -q $pbsqueue" >> $qsub1
	echo "#PBS -m ae" >> $qsub1
	echo "#PBS -M $email" >> $qsub1
	echo "$scriptdir/start_align_block.sh $runfile $TopOutputLogs/log.START_ALIGN_BLOCK.in $TopOutputLogs/log.START_ALIGN_BLOCK.ou $email $TopOutputLogs/qsub.main.aln" >> $qsub1
	#`chmod a+r $qsub1`               
	`qsub $qsub1 >> $TopOutputLogs/pbs.ALIGN`
	echo `date`
	echo "Note: realign module will be scheduled after align module ends"
	case="align_and_realign"  
 fi
 
 if [ `expr ${#case}` -lt 1 ]
 then
	MSG="Invalid value for parameter ANALYSIS=$analysis in configuration file."
	echo -e "$MSG\n\nDetails:\n\n$LOGS" # | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1; 
 fi

