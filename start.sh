#!/bin/bash

################################################################################################ 
# Program to calculate raw variants from human samples of WES short reads
# In order to run this pipeline please type at the command line
# start.sh <runfile>
################################################################################################
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 1 ]
then
        MSG="Parameter mismatch.\nRerun like this: $0 <runfile>\n"
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
        exit 1;
fi


echo -e "\n\n########################################################################################"
echo -e     "#############                BEGIN VARIANT CALLING WORKFLOW              ###############"
echo -e     "########################################################################################\n\n"

umask 0027
set -x
echo `date`	
scriptfile=$0
runfile=$1
if [ !  -s $runfile ]
then
   MSG="program=$0 stopped at line=$LINENO. $runfile configuration file not found."
   exit 1;
fi

echo -e "\n\n########################################################################################"
echo -e "#############                CHECKING PARAMETERS                         ###############"
echo -e "########################################################################################\n\n"

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )  
email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
inputformat=$( cat $runfile | grep -w INPUTFORMAT | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )        
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
dup_cutoff=$( cat $runfile | grep -w  DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w  MAP_CUTOFF | cut -d '=' -f2 )
paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
aligner_tool=$( cat $runfile | grep -w ALIGNERTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
aligner_mod=$( cat $runfile | grep -w ALIGNERMODULE | cut -d '=' -f2 )
samtools_mod=$( cat $runfile | grep -w SAMTOOLSMODULE | cut -d '=' -f2 )
vcftools_mod=$( cat $runfile | grep -w VCFTOOLSMODULE | cut -d '=' -f2 )
samblaster_mod=$( cat $runfile | grep -w SAMBLASTERMODULE | cut -d '=' -f2 )
sorttool_mod=$( cat $runfile | grep -w SORTMODULE | cut -d '=' -f2 )
markduplicates=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
gatk_mod=$( cat $runfile | grep -w GATKMODULE | cut -d '=' -f2 )        
gatk_dir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
picard_mod=$( cat $runfile | grep -w PICARDMODULE | cut -d '=' -f2 )
java_mod=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
indices=$( cat $runfile | grep -w CHRNAMES | cut -d '=' -f2 | tr ':' ' ' )
tabix_mod=$( cat $runfile | grep -w TABIXMODULE | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
queue=$( cat $runfile | grep -w PBSQUEUE | cut -d '=' -f2 )
analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

if [ `expr ${#tmpdir}` -lt 1  ]
then
	MSG="Invalid value specified for TMPDIR in the configuration file."
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -d  $refdir  ]
then
	MSG="Invalid value specified for REFGENOMEDIR=$refdir in the configuration file."
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -s  $refdir/$refgenome  ]
then
	MSG="Invalid value specified for REFGENOME=$refgenome in the configuration file."
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -s  $refdir/$dbSNP  ]
then
	MSG="Invalid value specified for DBSNP=$dbSNP in the configuration file."
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ ! -d  $gatk_dir  ]
then
	MSG="Invalid value specified for GATKDIR=$gatk_dir in the configuration file."
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi


if [ $inputformat != "FASTQ"  ]
then
    MSG="Incorrect value for INPUTFORMAT=$inputformat in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ $aligner_tool != "BWAMEM"  -a $aligner_tool != "BWA_MEM"  ]
then
    MSG="Incorrect value for ALIGNERTOOL=$aligner_tool in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $scriptdir ]
then
   MSG="SCRIPTDIR=$scriptdir directory not found"
   echo -e "Program $0 stopped at line=$LINENO.\n\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi
if [ -z $email ]
then
   MSG="Invalid value for parameter PBSEMAIL=$email in the configuration file"
   echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi

if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ] 
then
	MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
	echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

if [ `expr ${#dup_cutoff}` -lt 1 -o `expr ${#map_cutoff}` -lt 1 ]
then
   MSG="Invalid value for MAP_CUTOFF=$map_cutoff or for DUP_CUTOFF=$dup_cutoff  in the configuration file"
   echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi

if [ `expr ${#indices}` -lt 1 ]
then
   MSG="Invalid value for CHRNAMES in the configuration file"
   echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi


if [ $markduplicates != "NOVOSORT" -a $markduplicates != "SAMBLASTER" -a $markduplicates != "PICARD" ]
then
    MSG="Invalid value for parameter MARKDUPLICATESTOOL=$markduplicates  in the configuration file."
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -s $sampleinfo ]
then
    MSG="SAMPLEINFORMATION=$sampleinfo  file not found."
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ $numsamples -lt 1 ]
then
    MSG="SAMPLEINFORMATION=$sampleinfo  file is empty."
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;	
fi

if [ $paired -ne 1 -a $paired != "YES" ]
then
    MSG="Invalid value for parameter PAIRED=$paired in configuration file "
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

echo -e "\n\n########################################################################################"
echo -e "###########                      checking PBS params                      ##################"
echo -e "########################################################################################\n\n"

if [ `expr ${#thr}` -lt 1 ]
then
    thr=$PBS_NUM_PPN
fi

if [ `expr ${#nodes}` -lt 1 ]
then
    nodes=1
fi

if [ `expr ${#queue}` -lt 1 ]
then
    queue="default"
fi


echo -e "\n\n########################################################################################"
echo -e "###########                      checking tools                       ##################"
echo -e "########################################################################################\n\n"

`module load $aligner_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for ALIGNERMODULE=$aligner_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
`module unload $aligner_mod`

`module load $samblaster_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for SAMBLASTERMODULE=$samtools_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
`module unload $samblaster_mod`

`module load $samtools_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for SAMTOOLSMODULE=$samtools_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
`module unload $samtools_mod`

`module load $vcftools_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for VCFTOOLSMODULE=$vcftools_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
`module unload $samtools_mod`

`module load $sorttool_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for SORTMODULE=$sorttool_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
`module unload $sorttool_mod`

`module load $tabix_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for TABIXMODULE=$tabix_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
`module unload $tabix_mod`

`module load $gatk_mod`
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for GATKMODULE=$gatk_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi

`module unload $gatk_mod`
`module unload java`
`module load $java_mod`
`module load $picard_mod`

exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for PICARDMODULE=$picard_mod in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi



echo -e "\n\n########################################################################################"
echo -e "#############  Everything seems ok. Now setup/configure output folders and files   #########"
echo -e "########################################################################################\n\n"

if [ ! -d $outputdir/logs  ]
then
        # the output directory does not exist. create it
        mkdir -p $outputdir/logs
fi

if [ ! -d $outputdir/$deliverydir/docs  ]
then
        # the delivery directory does not exist. create it
	mkdir -p $outputdir/$deliverydir/docs
fi

`cp $runfile    $outputdir/$deliverydir/docs/runfile.txt`
`cp $sampleinfo $outputdir/$deliverydir/docs/sampleinfo.txt`
truncate -s 0   $outputdir/$deliverydir/docs/Summary.Report
truncate -s 0   $outputdir/$deliverydir/docs/QC_test_results.txt 

runfile=$outputdir/$deliverydir/docs/runfile.txt
TopOutputLogs=$outputdir/logs

truncate -s 0 $TopOutputLogs/pbs.ALIGN
truncate -s 0 $TopOutputLogs/pbs.MERGE

generic_qsub_header=$TopOutputLogs/qsubGenericHeader
truncate -s 0 $generic_qsub_header
echo "#!/bin/bash" > $generic_qsub_header
echo "#PBS -q $queue" >> $generic_qsub_header
echo "#PBS -m ae" >> $generic_qsub_header
echo "#PBS -M $email" >> $generic_qsub_header
echo "#PBS -l nodes=$nodes:ppn=$thr" >> $generic_qsub_header


echo -e "##### let's check that it worked and that the file was created                     ####"

if [ ! -s $generic_qsub_header ]
then 
    MSG="$generic_qsub_header is empty"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

echo -e "\n\n########################################################################################"
echo -e "################### Documenting progress on redmine with this message ##################"
echo -e "########################################################################################"
echo -e "##### the first part of the Report also needs to be stored in Summary.Report      ######"
echo -e "########################################################################################\n\n"



MSG="Variant calling workflow  started by username:$USER at: "$( echo `date` )
LOGS="Documentation about this run such as config files and results of QC tests will be placed in this folder:\n\n$outputdir/$deliverydir/docs/ \n\n"
echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
echo -e "$MSG\n\nDetails:\n\n$LOGS" >> $outputdir/$deliverydir/docs/Summary.Report


echo -e "\n\n########################################################################################"
echo -e "########################################################################################"
echo -e "########################################################################################"
echo -e "#####                               MAIN LOOP STARTS HERE                      #########"
echo -e "########################################################################################"
echo -e "#####  Trimming has been performed already                                     #########"
echo -e "#####  Alignment-dedup analysis: one qsub per sample                           #########"
echo -e "#####  Realignment-recalibration-variantCalling: 25 qsubs per sample, one per chr   ####"
echo -e "########################################################################################\n\n"

while read sampleLine
do
    if [ `expr ${#sampleLine}` -lt 1 ]
    then
	echo -e "\n\n########################################################################################"
	echo -e "##############                 skipping empty line        ##############################"
	echo -e "########################################################################################\n\n"
    else
	echo -e "\n\n########################################################################################"
	echo -e "#####         Processing next line $sampleLine                                ##########"
	echo -e "##### col1=sample_name col2=read1 col3=read2  including full paths            ##########"
	echo -e "##### sample_name will be used for directory namas and in RG line of BAM files##########"
	echo -e "########################################################################################\n\n"

	sample=$( echo "$sampleLine" | cut -d ' ' -f 1 )  
	FQ_R1=$( echo "$sampleLine" | cut -d ' '  -f 2 )
	FQ_R2=$( echo "$sampleLine" | cut -d ' ' -f 3 )

	if [ `expr ${#sample}` -lt 1 ]
	then
	     MSG="unable to parse line $sampleLine"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	     exit 1
	fi

	if [ `expr ${#FQ_R1}` -lt 1 ]
	then
	     MSG="unable to parse line $sampleLine"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	     exit 1
	elif [ ! -s $FQ_R1 ]
	then
	     MSG="$FQ_R1 read1 file not found"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                                          
	     exit 1                
	fi

	if [ `expr ${#FQ_R2}` -lt 1 ]
	then
	     MSG="unable to parse line $sampleLine"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	     exit 1
	elif [ ! -s $FQ_R2 ]
	then
	     MSG="$FQ_R2 read2 file not found"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	     exit 1                
	fi

	echo -e "\n\n########################################################################################"
	echo -e "###   Everything seems in order. Now creating folders where results will go  ###########"
	echo -e "########################################################################################\n\n"

	if [ -d $outputdir/$sample ]
	then
	     ### $outputdir/$sample already exists. Resetting it now. 
	     ### Maybe not. We already run trimming and we want to keep those results
	     ### rm -R $outputdir/$sample
	     mkdir -p $outputdir/$sample/align
	     mkdir -p $outputdir/$sample/realign
	     mkdir -p $outputdir/$sample/variant
	     mkdir -p $outputdir/$deliverydir/$sample
	else 
	     mkdir -p $outputdir/$sample/align
	     mkdir -p $outputdir/$sample/realign
	     mkdir -p $outputdir/$sample/variant
	     mkdir -p $outputdir/$deliverydir/$sample	     
	fi

	echo -e "\n\n########################################################################################"                
	echo -e "####   Launching Alignment script for SAMPLE $sample with R1=$FQ_R1 R2=$FQ_R2     ##########"
	echo -e "########################################################################################\n\n"

	qsub1=$TopOutputLogs/qsub.alignDedup.$sample
	cat $generic_qsub_header > $qsub1
	echo "#PBS -N alignDedup.$sample" >> $qsub1
	echo "#PBS -o $TopOutputLogs/log.alignDedup.$sample.ou" >> $qsub1
	echo "#PBS -e $TopOutputLogs/log.alignDedup.$sample.in" >> $qsub1
	echo "echo `date`" >> $qsub1
	echo "$scriptdir/align_dedup.sh $runfile $sample $FQ_R1 $FQ_R2 $TopOutputLogs/log.alignDedup.$sample.in $TopOutputLogs/log.alignDedup.$sample.ou $TopOutputLogs/qsub.alignDedup.$sample" >> $qsub1
	echo "echo `date`" >> $qsub1
	`chmod a+r $qsub1`               
	jobid=`qsub $qsub1` 
        `qhold -h u $jobid`
	echo $jobid >> $TopOutputLogs/pbs.ALIGN
	echo `date`

        if [ `expr ${#jobid}` -lt 1 ]
        then
	     MSG="unable to launch qsub align job for $sample. Exiting now"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	     exit 1        
        
        fi
	echo -e "\n\n#######  this jobid=$jobid will be used to hold execution of realign_varcall.sh     ########\n\n"

        if [ $analysis == "ALIGNMENT" -o $analysis == "ALIGN" -o $analysis == "ALIGN_ONLY" ]
	then
	    set +x; echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" >&2; set -x;
            # release all held jobs
            `qrls -h u $jobid`
        else




	   echo -e "\n\n########################################################################################"                
	   echo -e "####   Next loop2 for Launching Realign-Vcall script for SAMPLE $sample on all chr     #####"
	   echo -e "########################################################################################\n\n"

           alnjobid=$( echo $jobid | cut -d '.' -f 1 )
        
	   truncate -s 0 $TopOutputLogs/pbs.VCALL.$sample
	
           for chr in $indices
           do


		echo -e "\n\n########################################################################################"                
		echo -e "####   Realign-Vcall script for SAMPLE $sample chr=$chr                              #######"
		echo -e "########################################################################################\n\n"

		qsub1=$TopOutputLogs/qsub.realVcall.$sample.$chr
		cat $generic_qsub_header > $qsub1
		echo "#PBS -N realVcall.$sample.$chr" >> $qsub1
		echo "#PBS -o $TopOutputLogs/log.realVcall.$sample.$chr.ou" >> $qsub1
		echo "#PBS -e $TopOutputLogs/log.realVcall.$sample.$chr.in" >> $qsub1
                echo "#PBS -W depend=afterok:$alnjobid" >> $qsub1
 	        echo "echo `date`" >> $qsub1               
		echo "$scriptdir/realign_varcall_by_chr.sh $runfile $sample $chr $TopOutputLogs/log.realVcall.$sample.$chr.in $TopOutputLogs/log.realVcall.$sample.$chr.ou $TopOutputLogs/qsub.realVcall.$sample.$chr" >> $qsub1
	        echo "echo `date`" >> $qsub1
		`chmod a+r $qsub1`               
		realjobid=`qsub $qsub1` 
		echo $realjobid >> $TopOutputLogs/pbs.VCALL.$sample
		echo `date`
		
		if [ `expr ${#realjobid}` -lt 1 ]
		then
		     MSG="unable to launch qsub realign job for $sample. Exiting now"
		     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
		     exit 1        

		fi

           done
        
	   echo -e "\n\n########################################################################################"                
	   echo -e "####   Out of loop2. Now launching merge_vcfs script for SAMPLE $sample       ##########"
	   echo -e "########################################################################################\n\n"

           vcalljobids=$( cat $TopOutputLogs/pbs.VCALL.$sample | sed "s/\.[a-z]*//g" | tr "\n" ":" )

	   echo -e "\n\n### this list of jobids=[$vcalljobids] will be used to hold execution of merge_vcfs.sh #####\n\n"

	   qsub1=$TopOutputLogs/qsub.merge.$sample
	   cat $generic_qsub_header > $qsub1
	   echo "#PBS -N merge.$sample" >> $qsub1
	   echo "#PBS -o $TopOutputLogs/log.merge.$sample.ou" >> $qsub1
	   echo "#PBS -e $TopOutputLogs/log.merge.$sample.in" >> $qsub1
           echo "#PBS -W depend=afterok:$vcalljobids" >> $qsub1
	   echo "echo `date`" >> $qsub1        
	   echo "$scriptdir/merge_vcf_and_bam.sh $runfile $sample $TopOutputLogs/log.mergeVcf.$sample.in $TopOutputLogs/log.merge.$sample.ou $TopOutputLogs/qsub.merge.$sample" >> $qsub1
	   echo "echo `date`" >> $qsub1
	   `chmod a+r $qsub1`               
	   mergejobid=`qsub $qsub1` 
	   echo $mergejobid >> $TopOutputLogs/pbs.MERGE
	   echo `date`
	
           if [ `expr ${#mergejobid}` -lt 1 ]
              then
	        MSG="unable to launch qsub merge job for $sample. Exiting now"
	        echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	        exit 1        
           fi	
        fi # close the if statement checking whether the workflow end with alignment or not
        
   fi  # end non-empty line
done <  $sampleinfo

echo -e "\n\n########################################################################################"
echo -e "########################################################################################"
echo -e "#################           MAIN LOOP ENDS HERE                  #######################"
echo -e "########################################################################################"
echo -e "########################################################################################"
echo -e "#################     Now, we need to generate summary           #######################"
echo -e "########################################################################################"
echo -e "########################################################################################\n\n"

alljobids=$( cat $TopOutputLogs/pbs.MERGE | sed "s/\.[a-z]*//g" | tr "\n" ":" )

echo -e "\n\n### this list of jobids=[$alljobids] will be used to hold execution of summary.sh #####\n\n"

qsub2=$TopOutputLogs/qsub.summary
cat $generic_qsub_header > $qsub2
echo "#PBS -N Summary_vcall" >> $qsub2
echo "#PBS -o $TopOutputLogs/log.summary.ou" >> $qsub2
echo "#PBS -e $TopOutputLogs/log.summary.in" >> $qsub2
echo "#PBS -W depend=afterok:$alljobids " >> $qsub2
echo "$scriptdir/summary.sh $runfile $TopOutputLogs/log.summary.in $TopOutputLogs/log.summary.ou $TopOutputLogs/qsub.summary" >> $qsub2
`chmod a+r $qsub2`
lastjobid=`qsub $qsub2`
echo $lastjobid >> $TopOutputLogs/pbs.SUMMARY
echo `date`     


if [ `expr ${#lastjobid}` -lt 1 ]
then
     MSG="unable to launch qsub summary job. Exiting now"
     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
     exit 1        

fi
        
echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
