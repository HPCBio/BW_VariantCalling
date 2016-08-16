#!/bin/bash
# written in collaboration with Mayo Bioinformatics core group
#
# script to schedule qsub jobs for realign, recalibrate and/or variant-call the aligned file(s)
#########################################################################################################
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu

if [ $# != 6 ]
then
	MSG="parameter mismatch."
	echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
	exit 1;
fi

umask 0027
set -x
echo `date`
scriptfile=$0
outputdir=$1
runfile=$2
elog=$3
olog=$4
email=$5
qsubfile=$6
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

if [ ! -s $runfile ]
then
	MSG="$runfile configuration file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s 'Variant Calling Workflow failure message' "$redmine,$email"
	exit 1;
fi



set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;




reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
run_method=$( cat $runfile | grep -w RUNMETHOD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
launcherdir=$( cat $runfile | grep -w LAUNCHERDIR | cut -d '=' -f2 )
run_cmd=$( cat $runfile | grep -w LAUNCHERCMD | cut -d '=' -f2 )   #string value aprun|mpirun
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
indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
indelfile=$( cat $runfile | grep -w INDELFILE | cut -d '=' -f2 )
snpdir=$( cat $runfile | grep -w SNPDIR | cut -d '=' -f2 )
targetdir=$( cat $runfile | grep -w TARGETREGIONSDIR | cut -d '=' -f2 )
targetfile=$( cat $runfile | grep -w TARGETREGIONSFILE | cut -d '=' -f2 )
realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
omnisites=$( cat $runfile | grep -w OMNISITES | cut -d '=' -f2 ) 
recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
runverify=$( cat $runfile | grep -w RUNVERIFYBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 )
chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
indices=$( echo $chrindex | sed 's/:/ /g' )
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
splitRealign=$( cat $runfile | grep -w SPLIT2REALIGN | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
packageByChr=$( cat $runfile | grep -w SPLIT2REALIGNBYCHR | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
bash_cmd=`which bash`

set +x; echo -e "\n\n#############      checking analysis ###############\n\n" >&2; set -x;

if [ $analysis == "MULTIPLEXED" ]
then
	MSG="ANALYSIS=$analysis Program=$scriptfile Invalid pipeline program for this type of analysis. This program is for the NON-MULTIPLEXED case only"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n#############      checking input type ###############\n\n" >&2; set -x;

if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
	input_type="WGS"
	pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
	pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
elif [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
	input_type="WES"
	pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
	pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
else
	MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
	echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n#############   split/not split bam by region and then  realign###############\n\n" >&2; set -x;

if [ $splitRealign != "1" -a $splitRealign != "0" -a $splitRealign != "YES" -a $splitRealign != "NO" ]
then
	MSG="Invalid value for SPLIT2REALIGN=$splitRealign"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ $splitRealign == "1" ]
then
	splitRealign="YES"
fi
if [ $splitRealign == "0" ]
then
	splitRealign="NO"
fi  

if [ $splitRealign == "YES" -a `expr ${#packageByChr}` -lt 1 ]
then
	packageByChr="YES"
elif [ $splitRealign == "YES" -a $packageByChr == "1" ]
then
	packageByChr="YES"
elif [ $splitRealign == "YES" -a $packageByChr == "0" ]
then
	packageByChr="NO"
elif [ $splitRealign == "YES" -a $packageByChr != "YES" -a $packageByChr != "NO" ]
then
	MSG="Invalid value for SPLIT2REALIGNBYCHR=$packageByChr"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n#############      checking cleanup options   ###############\n\n" >&2; set -x;


if [ $cleanupflag != "1" -a $cleanupflag != "0" -a $cleanupflag != "YES" -a $cleanupflag != "NO" ]
then
	MSG="Invalid value for REMOVETEMPFILES=$cleanupflag"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ $cleanupflag == "1" ]
then
	cleanupflag="YES"
fi
if [ $cleanupflag == "0" ]
then
	cleanupflag="NO"
fi


set +x; echo -e "\n\n#############skip/include variant calling module ###############\n\n" >&2; set -x;

if [ $skipvcall != "1" -a $skipvcall != "0" -a $skipvcall != "YES" -a $skipvcall != "NO" ]
then
	MSG="Invalid value for SKIPVCALL=$skipvcall"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ $skipvcall == "1" ]
then
	skipvcall="YES"
fi
if [ $skipvcall == "0" ]
then
	skipvcall="NO"
fi

set +x; echo -e "\n\n#############skip/include QC step for cleanedBAM with verifyBamID  o###############\n\n" >&2; set -x;

if [ $runverify != "1" -a $runverify != "0" -a $runverify != "YES" -a $runverify != "NO" ]
then
	MSG="Invalid value for RUNVERIFYBAM=$runverify"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ $runverify == "1" ]
then
	runverify="YES"
fi
if [ $runverify == "0" ]
then
	runverify="NO"
fi   

if [ $run_cmd == "aprun" ]
then        
	run_cmd="aprun -n "
else        
	run_cmd="mpirun -np " 
fi


set +x; echo -e "\n\n#############version of java for gatk ###############\n\n" >&2; set -x;

if [ -z $javadir ]
then
	MSG="Value for JAVADIR must be specified in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n###########checking directories  and  tools ###############\n\n" >&2; set -x;


if [ ! -d $outputdir ]
then
	MSG="$outputdir ROOT directory for this run of the pipeline not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -d $picardir ]
then
	MSG="$picardir picard directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ ! -d $samdir ]
then
	MSG="$samdir samtools directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ ! -d $gatk ]
then
	MSG="$gatk GATK directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

set +x; echo -e "\n\n###########checking callsets needed for analyses ###############\n\n" >&2; set -x;

if [ ! -d $refdir ]
then
	MSG="$refdir reference genome directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ ! -s $refdir/$ref ]
then
	MSG="$ref reference genome not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -s $refdir/$dbSNP ]
then
	MSG="$dbSNP DBSNP file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi   
if [ ! -d $refdir/$indeldir ]
then
	MSG="$indeldir indel directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	exit 1;
fi
if [ ! -d $refdir/$snpdir ]
then
	MSG="$snpdir SNPs directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	exit 1;
fi

if [ ! -s $indelfile ]
then
	MSG="$indelfile INDELS file not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi   

if [ $runverify == "YES" -a ! -s $omnisites ]
then
	MSG="OMNISITES=$omnisites file not found. A file must be specified if this parameter has been specified too RUNVERIFYBAM=$runverify"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ `expr ${#targetfile}` -gt 1 -a $input_type == "WES" -a ! -s $targetfile ]
then
	MSG="TARGETFILE=$targetfile file not found. A file must be specified if this parameter has been specified too INPUTTYPE==WES"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi 


set +x; echo -e "\n\n" >&2;   
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "###########################      params ok. Now creating/resetting logs and folders        #########" >&2
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;  
  

TopOutputLogs=$outputdir/logs
RealignOutputLogs=$outputdir/logs/realign
VcallOutputLogs=$outputdir/logs/variant

if [ ! -d $TopOutputLogs ]
then
	MSG="$TopOutputLogs directory not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi
if [ ! -d $RealignOutputLogs ]
then
	mkdir -p $RealignOutputLogs
fi
if [ ! -d $VcallOutputLogs ]
then
	mkdir -p $VcallOutputLogs
fi

pipeid=$( cat $TopOutputLogs/pbs.CONFIGURE )
truncate -s 0 $TopOutputLogs/pbs.REALRECAL 
truncate -s 0 $TopOutputLogs/pbs.VCALL
truncate -s 0 $TopOutputLogs/pbs.SUMMARY


set +x; echo -e "\n\n" >&2; 
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "#######   LOOP to generate regions, targets, known/knownSites                           ############" >&2
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "#######   target is the array of targeted regions for WES ONLY which will be used by all GATK commands"  >&2
echo "#######   region is the array of known variants per chr"  >&2
echo "#######   realparms is the array with indels per chr which will be used for  gatk-IndelRealigner" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x; 

if [ $splitRealign == "YES" ]
then

	set +x; echo -e "\n\n########### samples will be split by region, therefore we need to create these files ###############\n\n" >&2; set -x;
	
	echo `date`
	i=1

	for chr in $indices
	do
		cd $refdir/$snpdir
		region[$i]=$( find $PWD -type f -name "${chr}.*.vcf.gz" | sed "s/^/:knownSites:/g" | tr "\n" ":" )
		cd $refdir/$indeldir
		realparms[$i]=$( find $PWD -type f -name "${chr}.*.vcf" | tr "\n" ":" )

		if [ -d $targetdir -a $input_type == "WES" ]
		then
			set +x; echo -e "\n\n###### let's check to see if we have a bed file for this chromosme ##############\n\n" >&2; set -x;

			if [ -s $targetdir/${chr}.bed ]
			then
				target[$i]=$targetdir/${chr}.bed
			fi
		fi
		(( i++ ))
	done

	echo `date`
fi

set +x; echo -e "\n\n" >&2;
echo -e "\n####################################################################################################" >&2 
echo "####################################################################################################" >&2
echo "################ Loop to check that aligned files exist AND to  Create output folders   ############" >&2     
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


if [ -s $outputdir/SAMPLENAMES_multiplexed.list ]
then 
	TheInputFile=$outputdir/SAMPLENAMES_multiplexed.list
else
	TheInputFile=$outputdir/SAMPLENAMES.list
fi


while read SampleLine
do
	if [ `expr ${#SampleLine}` -gt 1 ]
	then

  		set +x; echo -e "\n\n###### processing next non-empty line in SAMPLENAMES_multiplexed.list ##############\n\n" >&2; set -x;

  		sample=$( echo "$SampleLine" | cut -f 1 )

  		set +x; echo -e "\n\n###### checking aligned bam for realigning-recalibrating sample: $sample      ##############\n\n" >&2; set -x;

  		alignedfile=`find $outputdir/$sample/align -name "*.wdups.sorted.bam"`

	        set +x; echo -e "\n\n###### # now check that there is only one bam file      ##############\n\n" >&2; set -x;

	        aligned_bam=$outputdir/${sample}/align/${alignedfile}
	        aligned_bam=( $aligned_bam ) # recast variable as an array and count the number of members
	        if [ ${#aligned_bam[@]} -ne 1 ]
                then
                    MSG="more than one bam file found in $outputdir/${sample}/align"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                    exit 1;
                fi	  

  		set +x; echo -e "\n\n###### create folders for the analyses     ##############\n\n" >&2; set -x;

  		if [ ! -d $outputdir/${sample}/realign ]
  		then
      			set +x; echo -e "\n\n###### creating output folders for sample=$sample    ##############\n\n" >&2; 

      			mkdir -p $outputdir/${sample}/realign/logs
      			mkdir -p $outputdir/${sample}/variant/logs

		else
			set +x; echo -e "\n\n###### resetting logs        ##############\n\n" >&2; set -x;

			rm $outputdir/${sample}/realign/logs/*
			rm $outputdir/${sample}/variant/logs/*
  		fi

	fi  # end processing non-empty lines
done  <  $TheInputFile
# end loop over samples



set +x; echo -e "\n\n" >&2; 
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "####### BLOCK to generate and package cmds for realign recalibrate and varcall   " >&2
echo "####################################################################################################" >&2
echo "####### CASE1: NO SAMPLE SPLITTING, LAUNCHER IS USED, there will be only one anisimov job for real-recalibrate and one for varcall " >&2
echo "####### CASE2: NO SAMPLE SPLITTING, LAUNCHER IS NOT USED, there will be two qsub jobs per sample " >&2
echo "####### CASE3: SAMPLE SPLITTING BY REGION, LAUNCHER IS ALWAYS USED,PACKAGE BY CHR " >&2
echo "####### CASE4: SAMPLE SPLITTING BY REGION, LAUNCHER IS ALWAYS USED,PACKAGE BY SAMPLE " >&2	
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


if [ $splitRealign == "NO" ]
then

	set +x; echo -e "\n\n###### SPLITREALIGN=NO. Samples will not be split by region for the analysis ##############\n\n" >&2; set -x;
	set +x; echo -e "\n\n###### CASE1 and CASE2  will be covered here   ##############\n\n" >&2; set -x;

	if [ $run_method == "LAUNCHER" ]
	then
		set +x; echo -e "\n\n###### LAUNCHER Prep work: Create and/or reset anisimov files      ##############\n\n" >&2; set -x;
		RealignOutputLogs=$outputdir/logs/realign
		VcallOutputLogs=$outputdir/logs/variant

		qsubRealrecalAnisimov=$RealignOutputLogs/qsub.realrecal.AnisimovJoblist
		qsubVcallAnisimov=$VcallOutputLogs/qsub.vcall.AnisimovJoblist
		realrecalAnisimov=$RealignOutputLogs/realrecal.AnisimovJoblist			
		vcallAnisimov=$VcallOutputLogs/vcall.AnisimovJoblist

		truncate -s 0 $realrecalAnisimov
		truncate -s 0 $vcallAnisimov
	fi

	set +x; echo -e "\n\n###### Only One loop over samples is needed to create commands and populate files   ##############\n\n" >&2; set -x;	
	samplecounter=0 
	while read SampleLine
	do
		if [ `expr ${#SampleLine}` -gt 1 ]
		then

			set +x; echo -e "\n\n###### processing next non-empty line in SAMPLENAMES_multiplexed.list ##############\n\n" >&2; set -x;


			set +x; echo -e "\n\n###### sample name and input file.     ##############\n\n" >&2; set -x;

			sample=$( echo "$SampleLine" | cut -f 1 )
			aligned_bam=`find $outputdir/$sample/align -name "*.wdups.sorted.bam"`   

			set +x; echo -e "\n\n###### define output directories and files for results, cmds, logs, etc.     ##############\n\n" >&2; set -x;
		      
			if [ $run_method == "LAUNCHER" ]
			then
				set +x; echo -e "\n\n###### LAUNCHER will be used. Log folders will have these paths     ##############\n\n" >&2; set -x;
				#RealignOutputLogs=$outputdir/logs/realign
				#VcallOutputLogs=$outputdir/logs/variant
				qsubRealrecal=$qsubRealrecalAnisimov
				qsubVcall=$qsubVcallAnisimov
			else
				set +x; echo -e "\n\n###### LAUNCHER will NOT be used. Log folders will have these paths  ##############\n\n" >&2; set -x;
				RealignOutputLogs=$outputdir/${sample}/realign/logs
				VcallOutputLogs=$outputdir/${sample}/variant/logs
				qsubRealrecal=$RealignOutputLogs/qsub.Realrecal.${sample}
				qsubVcall=$VcallOutputLogs/qsub.Vcall.${sample}				
			fi			

			set +x; echo -e "\n\n###### Output results will go on here  ##############\n\n" >&2; set -x;
			
			RealignOutputDir=$outputdir/${sample}/realign
			VcallOutputDir=$outputdir/${sample}/variant
			realignedOutputfile=$RealignOutputDir/${sample}.realigned.bam
			recalibratedOutputfile=$RealignOutputDir/${sample}.recalibrated.calmd.bam
			vcallOutputFile=$VcallOutputDir/${sample}.raw.g.vcf

			set +x; echo -e "\n\n###### Logs and files with command lines will go here  ##############\n\n" >&2; set -x;
			
			jobfileRealrecal=$RealignOutputLogs/Realrecal.${sample}.jobfile
			jobfileVcall=$VcallOutputLogs/Vcall.${sample}.jobfile
			realignFailedlog=$RealignOutputDir/logs/FAILED_Realrecal.${sample}              
			vcallFailedlog=$VcallOutputDir/logs/FAILED_Vcall.${sample}

			set +x; echo -e "\n\n###### reset those files  ##############\n\n" >&2; set -x;
			
			truncate -s 0 $jobfileRealrecal
			truncate -s 0 $jobfileVcall			
			truncate -s 0 $realignFailedlog
			truncate -s 0 $vcallFailedlog


			set +x; echo -e "\n\n###### which target file to use for WES data in all GATK commands     ##############\n\n" >&2; set -x;
			
			target=$targetfile
			if [ `expr ${#target}` -lt 1 -a $input_type == "WES"  ]
			then
				target="NOTARGET"         
			fi

			set +x; echo -e "\n\n###### ready to put together the command for the realignment-recalibration analysis   ##############\n\n" >&2; set -x;

			echo "$scriptdir/realrecal_sample.sh $RealignOutputDir $aligned_bam $realignedOutputfile $recalibratedOutputfile $target $runfile $realignFailedlog $email $qsubRealrecal" > $jobfileRealrecal

			if [ $skipvcall == "NO" ]
			then
				set +x; echo -e "\n\n######  ready to put together the command for the variant calling analysis    ########\n\n" >&2; set -x;
				echo "$scriptdir/vcallgatk_sample.sh $recalibratedOutputfile $vcallOutputFile $VcallOutputDir $target $runfile $vcallFailedlog $email $qsubVcall" > $jobfileVcall
			fi

			set +x; echo -e "\n\n###### populate anisimov list if LAUNCHER will be used   ##############\n\n" >&2; set -x;
			
			if [ $run_method == "LAUNCHER" ]
			then
				Realjobfilename=$( basename $jobfileRealrecal )
				echo "$RealignOutputLogs $Realjobfilename" >> $realrecalAnisimov
				
				if [ $skipvcall == "NO" ]
				then
					Vcalljobfilename=$( basename $jobfileVcall )
					echo "$VcallOutputLogs $Vcalljobfilename" >> $vcallAnisimov
				fi
			fi						

			(( samplecounter++ ))  # this counter will be linked to numnodes we need to allocate to the launcher
		      
		fi # done processing non-empty lines

		set +x; echo -e "\n\n###### Done with this sample. Next one please   ##############\n\n" >&2; set -x; 

	done <  $TheInputFile
	
	set +x; echo -e "\n\n###### END SPLITREALIGN=NO ##############\n\n"  >&2; set -x;


elif [ $splitRealign == "YES" -a $packageByChr == "YES" ]
then

	set +x; echo -e "\n\n###### CASE3 SPLITREALIGN=YES. PACKAGEBYCHR=YES. THERE WILL BE ONE LAUNCHER PER REGION ##############\n\n" >&2; set -x;

	####stuff goes here

	set +x; echo -e "\n\n###### END SPLITREALIGN=YES. PACKAGEBYCHR=YES ##############\n\n" >&2; set -x;
	
elif [ $splitRealign == "YES" -a $packageByChr == "NO" ]
then
	set +x; echo -e "\n\n###### CASE4 SPLITREALIGN=YES. PACKAGEBYCHR=NO. THERE WILL BE ONE LAUNCHER PER SAMPLE##############\n\n" >&2; set -x;

	####stuff goes here
	
	set +x; echo -e "\n\n###### END SPLITREALIGN=YES. PACKAGEBYCHR=NO. ##############\n\n" >&2; set -x;
	
else
	MSG="THIS CASE IS NOT YET COVERED BY THIS SCRIPT. SPLIT2REALIGN=$splitRealign SPLIT2REALIGNBYCHR=$packageByChr RUNMETHOD=$run_method. EXITING NOW"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;

fi
  
# at the end of this set of nested loops, the variables chromosomecounter and samplecounter
# reflect, respectively, number_of_chromosomes+1 and number_of_samples+1,
# which is exactly the number of nodes required for anisimov launcher 



set +x; echo -e "\n\n" >&2; 
echo -e "\n####################################################################################################" >&2
echo -e "\n##########     END BLOCK to generate cmds and populate anisimov joblists          ###########" >&2
echo -e "\n####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;


#set +x; echo -e "\n ### update autodocumentation script ### \n"; set -x;
#echo -e "# @begin RealignRecalibrate_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
#echo -e "   # @in sample_chr @as aligned_bam_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
#RealignedBAMTemplate="{SampleName}/realign/{chromosome}.realrecal.${SampleName}.output.bam"
#echo -e "   # @out realrecal  @as  realigned_bam_per_chromosome @URI ${RealignedBAMTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
#echo -e "# @end RealignRecalibrate_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh

#echo -e "# @begin VariantCalling_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
#echo -e "   # @in  realrecal  @as  realigned_bam_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh
#VariantTemplate=${RealignedBAMTemplate}.raw.all.vcf
#echo -e "   # @out vcf @as output_variants @URI ${VariantTemplate}" >> $outputdir/WorkflowAutodocumentationScript.sh
#echo -e "# @end VariantCalling_per_chromosome" >> $outputdir/WorkflowAutodocumentationScript.sh

#echo -e "# @out vcf @as output_variants @URI sample_name/variant/chr_name.vcf" >> $outputdir/WorkflowAutodocumentationScript.sh
#WorkflowName=`basename $outputdir`
#echo -e "# @end $WorkflowName" >> $outputdir/WorkflowAutodocumentationScript.sh



set +x; echo -e "\n\n" >&2; 
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "####### BLOCK TO GENERATE AND SCHEDULE QSUBS   " >&2
echo "####################################################################################################" >&2
echo "####### CASE1: NO SAMPLE SPLITTING, LAUNCHER IS USED, there will be only one anisimov job for real-recalibrate and one for varcall " >&2
echo "####### CASE2: NO SAMPLE SPLITTING, LAUNCHER IS NOT USED, there will be two qsub jobs per sample " >&2
echo "####### CASE3: SAMPLE SPLITTING BY REGION, LAUNCHER IS ALWAYS USED,PACKAGE BY CHR " >&2
echo "####### CASE4: SAMPLE SPLITTING BY REGION, LAUNCHER IS ALWAYS USED,PACKAGE BY SAMPLE " >&2	
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

echo -e "\n\n" >&2; set -x;

if [ $splitRealign == "NO" -a $run_method == "LAUNCHER" ]
then

	set +x; echo -e "\n\n###### CASE1 splitRealign=NO. run_methods IS LAUNCHER ##############\n\n" >&2; set -x;

	set +x; echo -e "\n # run_method is LAUNCHER. scheduling the real-recalibrate Launcher\n" >&2; set -x

	let "numnodes = $samplecounter + 1"  # we need an additional node for the launcher

	RealignOutputLogs=$outputdir/logs/realign
	VcallOutputLogs=$outputdir/logs/variant
	
	qsubRealrecalAnisimov=$RealignOutputLogs/qsub.realrecal.AnisimovJoblist
	qsubVcallAnisimov=$VcallOutputLogs/qsub.vcall.AnisimovJoblist
	realrecalAnisimov=$RealignOutputLogs/realrecal.AnisimovJoblist			
	vcallAnisimov=$VcallOutputLogs/vcall.AnisimovJoblist


	echo "#!/bin/bash" > $qsubRealrecalAnisimov
	echo "#PBS -A $pbsprj" >> $qsubRealrecalAnisimov
	echo "#PBS -N ${pipeid}_realrecal_Anisimov" >> $qsubRealrecalAnisimov
	echo "#PBS -l walltime=$pbscpu" >> $qsubRealrecalAnisimov
	echo "#PBS -l nodes=$numnodes:ppn=$thr" >> $qsubRealrecalAnisimov
	echo "#PBS -o $RealignOutputLogs/log.realrecal.AnisimovJoblist.ou" >> $qsubRealrecalAnisimov
	echo "#PBS -e $RealignOutputLogs/log.realrecal.AnisimovJoblist.in" >> $qsubRealrecalAnisimov
	echo "#PBS -q $pbsqueue" >> $qsubRealrecalAnisimov
	echo "#PBS -m ae" >> $qsubRealrecalAnisimov
	echo "#PBS -M $email" >> $qsubRealrecalAnisimov

	echo "$run_cmd $numnodes -env OMP_NUM_THREADS=$thr $launcherdir/scheduler.x $realrecalAnisimov $bash_cmd > ${realrecalAnisimov}.log" >> $qsubRealrecalAnisimov

	echo "exitcode=\$?" >> $qsubRealrecalAnisimov
	echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubRealrecalAnisimov
	echo "   echo -e \"\n\n RealrecalAnisimov failed with exit code = \$exitcode \n logfile=$RealignOutputLogs/log.realrecal.AnisimovJoblist.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubRealrecalAnisimov
	echo "   exit 1" >> $qsubRealrecalAnisimov
	echo "fi" >> $qsubRealrecalAnisimov

	RealrecalAnisimovJoblistId=`qsub $qsubRealrecalAnisimov`
	echo $RealrecalAnisimovJoblistId >> $TopOutputLogs/pbs.REALRECAL 

	set +x; echo -e "\n # run_method is LAUNCHER. scheduling the vcall Launcher\n" >&2; set -x
	
	if [ $skipvcall == "NO" ]
	then

		echo "#!/bin/bash" > $qsubVcallAnisimov
		echo "#PBS -A $pbsprj" >> $qsubVcallAnisimov
		echo "#PBS -N ${pipeid}_Vcall_Anisimov" >> $qsubVcallAnisimov
		echo "#PBS -l walltime=$pbscpu" >> $qsubVcallAnisimov
		echo "#PBS -l nodes=$numnodes:ppn=$thr" >> $qsubVcallAnisimov
		echo "#PBS -o $VcallOutputLogs/log.Vcall.AnisimovJoblist.ou" >> $qsubVcallAnisimov
		echo "#PBS -e $VcallOutputLogs/log.Vcall.AnisimovJoblist.in" >> $qsubVcallAnisimov
		echo "#PBS -q $pbsqueue" >> $qsubVcallAnisimov
		echo "#PBS -m ae" >> $qsubVcallAnisimov
		echo "#PBS -M $email" >> $qsubVcallAnisimov
		echo "#PBS -W depend=afterok:$RealrecalAnisimovJoblistId" >> $qsubVcallAnisimov
		
		echo "$run_cmd $numnodes -env OMP_NUM_THREADS=$thr $launcherdir/scheduler.x $vcallAnisimov $bash_cmd > ${vcallAnisimov}.log" >> $qsubVcallAnisimov

		echo "exitcode=\$?" >> $qsubVcallAnisimov
		echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubVcallAnisimov
		echo "   echo -e \"\n\n RealrecalAnisimov failed with exit code = \$exitcode \n logfile=$VcallOutputLogs/log.Vcall.AnisimovJoblist.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubVcallAnisimov
		echo "   exit 1" >> $qsubVcallAnisimov
		echo "fi" >> $qsubVcallAnisimov

		VcallAnisimovJoblistId=`qsub $qsubVcallAnisimov`
		echo $VcallAnisimovJoblistId >> $TopOutputLogs/pbs.VCALL 

	
	fi  #end if skipvcall
	
	set +x; echo -e "\n\n###### END CASE1 splitRealign=NO. run_methods IS LAUNCHER ##############\n\n" >&2; set -x;	

elif [ $splitRealign == "NO" -a $run_method != "LAUNCHER" ]
then

	set +x; echo -e "\n\n###### CASE2 splitRealign=NO. run_methods IS NOT LAUNCHER ##############\n\n" >&2; set -x;

	set +x; echo -e "\n\n###### One loop over samples is needed to launch the jobs   ##############\n\n" >&2; set -x;	

	while read SampleLine
	do
		if [ `expr ${#SampleLine}` -gt 1 ]
		then

			RealignOutputLogs=$outputdir/${sample}/realign/logs
			VcallOutputLogs=$outputdir/${sample}/variant/logs
			qsubRealrecal=$RealignOutputLogs/qsub.Realrecal.${sample}
			qsubVcall=$VcallOutputLogs/qsub.Vcall.${sample}				

			jobfileRealrecal=$RealignOutputLogs/Realrecal.${sample}.jobfile
			jobfileVcall=$VcallOutputLogs/Vcall.${sample}.jobfile
			realignFailedlog=$RealignOutputDir/logs/FAILED_Realrecal.${sample}              
			vcallFailedlog=$VcallOutputDir/logs/FAILED_Vcall.${sample}


			echo "#!/bin/bash" > $qsubRealrecal
			echo "#PBS -A $pbsprj" >> $qsubRealrecal
			echo "#PBS -N ${pipeid}_realrecal_${sample}" >> $qsubRealrecal
			echo "#PBS -l walltime=$pbscpu" >> $qsubRealrecal
			echo "#PBS -l nodes=1:ppn=$thr" >> $qsubRealrecal
			echo "#PBS -o $RealignOutputLogs/log.realrecal.${sample}.ou" >> $qsubRealrecal
			echo "#PBS -e $RealignOutputLogs/log.realrecal.${sample}.in" >> $qsubRealrecal
			echo "#PBS -q $pbsqueue" >> $qsubRealrecal
			echo "#PBS -m ae" >> $qsubRealrecal
			echo "#PBS -M $email" >> $qsubRealrecal

			cat $jobfileRealrecal >> $qsubRealrecal

			echo "exitcode=\$?" >> $qsubRealrecal
			echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubRealrecal
			echo "   echo -e \"\n\n Realrecal $sample failed with exit code = \$exitcode \n logfile=$TopOutputLogs/log.realrecal.AnisimovJoblist.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubRealrecal
			echo "   exit 1" >> $qsubRealrecal
			echo "fi" >> $qsubRealrecal

			RealrecalJobId=`qsub $qsubRealrecal`
			echo $RealrecalJobId >> $TopOutputLogs/pbs.REALRECAL 

			set +x; echo -e "\n # run_method is LAUNCHER. scheduling the vcall Launcher\n" >&2; set -x

			if [ $skipvcall == "NO" ]
			then

				echo "#!/bin/bash" > $qsubVcall
				echo "#PBS -A $pbsprj" >> $qsubVcall
				echo "#PBS -N ${pipeid}_Vcall_${sample}" >> $qsubVcall
				echo "#PBS -l walltime=$pbscpu" >> $qsubVcall
				echo "#PBS -l nodes=1:ppn=$thr" >> $qsubVcall
				echo "#PBS -o $VcallOutputLogs/log.Vcall.${sample}.ou" >> $qsubVcall
				echo "#PBS -e $VcallOutputLogs/log.Vcall.${sample}.in" >> $qsubVcall
				echo "#PBS -q $pbsqueue" >> $qsubVcall
				echo "#PBS -m ae" >> $qsubVcall
				echo "#PBS -M $email" >> $qsubVcall
				echo "#PBS -W depend=afterok:$RealrecalJobId" >> $qsubVcall
				
				cat  $jobfileVcall >> $qsubVcall
				
				echo "exitcode=\$?" >> $qsubVcall
				echo -e "if [ \$exitcode -ne 0 ]\nthen " >> $qsubVcall
				echo "   echo -e \"\n\n Vcall ${sample} failed with exit code = \$exitcode \n logfile=$VcallOutputLogs/log.Vcall.${sample}.in\n\" | mail -s \"[Task #${reportticket}]\" \"$redmine,$email\"" >> $qsubVcall
				echo "   exit 1" >> $qsubVcall
				echo "fi" >> $qsubVcall

				VcallJobId=`qsub $qsubVcall`
				echo $VcallJobId >> $TopOutputLogs/pbs.VCALL 


			fi  #end if skipvcall

		fi # done processing non-empty lines

		set +x; echo -e "\n\n###### bottom of the loop over samples. Next one please "  >&2; set -x; 

	done <  $TheInputFile

	set +x; echo -e "\n\n###### END CASE2 splitRealign=NO. run_methods IS NOT LAUNCHER ##############\n\n" >&2; set -x;
	
elif [ $splitRealign == "YES" -a $packageByChr == "YES" ]
then

	set +x; echo -e "\n\n###### CASE3 SPLITREALIGN=YES. PACKAGEBYCHR=YES. THERE WILL BE ONE LAUNCHER PER REGION ##############\n\n" >&2; set -x;

	####stuff goes here

	set +x; echo -e "\n\n###### END SPLITREALIGN=YES. PACKAGEBYCHR=YES ##############\n\n" >&2; set -x;
	
elif [ $splitRealign == "YES" -a $packageByChr == "NO" ]
then
	set +x; echo -e "\n\n###### CASE4 SPLITREALIGN=YES. PACKAGEBYCHR=NO. THERE WILL BE ONE LAUNCHER PER SAMPLE##############\n\n" >&2; set -x;

	####stuff goes here
	
	set +x; echo -e "\n\n###### END SPLITREALIGN=YES. PACKAGEBYCHR=NO. ##############\n\n" >&2; set -x;

fi

set +x; echo -e "\n\n" >&2; 
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "####### BLOCK to wrap up and produce summary table                      ############################" >&2
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

if [ $skipvcall == "NO" ]
then
      summarydependids=$( cat $TopOutputLogs/pbs.VCALL  | tr "\n" ":" | sed "s/^://" | sed "s/:$//" )     #for job dependency argument
else
      summarydependids=$( cat $TopOutputLogs/pbs.REALRECAL  | tr "\n" ":" | sed "s/^://" | sed "s/:$//" ) #for job dependency argument
fi

lastjobid=""
cleanjobid=""

if [ $cleanupflag == "YES" ]
then
       qsub_cleanup=$TopOutputLogs/qsub.cleanup
       echo "#PBS -A $pbsprj" >> $qsub_cleanup
       echo "#PBS -N ${pipeid}_cleanup" >> $qsub_cleanup
       echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
       echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
       echo "#PBS -o $TopOutputLogs/log.cleanup.ou" >> $qsub_cleanup
       echo "#PBS -e $TopOutputLogs/log.cleanup.in" >> $qsub_cleanup
       echo "#PBS -q $pbsqueue" >> $qsub_cleanup
       echo "#PBS -m a" >> $qsub_cleanup
       echo "#PBS -M $email" >> $qsub_cleanup
       echo "#PBS -W depend=afterok:$summarydependids" >> $qsub_cleanup
       echo "$scriptdir/cleanup.sh $outputdir $analysis $TopOutputLogs/log.cleanup.in $TopOutputLogs/log.cleanup.ou $email $TopOutputLogs/qsub.cleanup" >> $qsub_cleanup
       #`chmod a+r $qsub_cleanup`
       cleanjobid=`qsub $qsub_cleanup`
       echo $cleanjobid >> $TopOutputLogs/pbs.CLEANUP
fi

`sleep 10s`

set +x; echo -e "\n\n###### Summary report when all jobs finished successfully ##############\n\n" >&2; set -x;

qsub_summary=$TopOutputLogs/qsub.summary.allok
echo "#PBS -A $pbsprj" >> $qsub_summary
echo "#PBS -N ${pipeid}_summaryok" >> $qsub_summary
echo "#PBS -l walltime=01:00:00" >> $qsub_summary # 1 hour should be more than enough
echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
echo "#PBS -o $TopOutputLogs/log.summary.ou" >> $qsub_summary
echo "#PBS -e $TopOutputLogs/log.summary.in" >> $qsub_summary
echo "#PBS -q $pbsqueue" >> $qsub_summary
echo "#PBS -m a" >> $qsub_summary
echo "#PBS -M $email" >> $qsub_summary

if [ $cleanupflag == "YES" ]
then
       echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub_summary
else
       echo "#PBS -W depend=afterok:$summarydependids" >> $qsub_summary
fi
echo "$scriptdir/summary.sh $runfile $email exitok $reportticket"  >> $qsub_summary
#`chmod a+r $qsub_summary`
lastjobid=`qsub $qsub_summary`
echo $lastjobid >> $TopOutputLogs/pbs.SUMMARY


set +x; echo -e "\n\n###### Summary report when at least one job failed ##############\n\n" >&2; set -x;

qsub_summary=$TopOutputLogs/qsub.summary.afterany
echo "#PBS -A $pbsprj" >> $qsub_summary
echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub_summary
echo "#PBS -l walltime=01:00:00" >> $qsub_summary # 1 hour should be more than enough
echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
echo "#PBS -o $TopOutputLogs/log.summary.afterany.ou" >> $qsub_summary
echo "#PBS -e $TopOutputLogs/log.summary.afterany.in" >> $qsub_summary
echo "#PBS -q $pbsqueue" >> $qsub_summary
echo "#PBS -m a" >> $qsub_summary
echo "#PBS -M $email" >> $qsub_summary

if [ $cleanupflag == "YES" ]
then
       echo "#PBS -W depend=afternotok:$cleanjobid" >> $qsub_summary
else
       echo "#PBS -W depend=afternotok:$summarydependids" >> $qsub_summary
fi
echo "$scriptdir/summary.sh $runfile $email exitnotok $reportticket"  >> $qsub_summary
#`chmod a+r $qsub_summary`
badjobid=`qsub $qsub_summary`
echo $badjobid >> $TopOutputLogs/pbs.SUMMARY


set +x; echo -e "\n\n" >&2; 
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo "####### DONE. EXITING NOW                                               ############################" >&2
echo "####################################################################################################" >&2
echo "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

