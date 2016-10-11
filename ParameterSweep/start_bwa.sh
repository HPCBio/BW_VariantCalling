#!/bin/bash

################################################################################################ 
# Program to evaluate the effect of changing the parameters in alignment of  samples of WES short reads
# In order to run this pipeline please type at the command line
# start_bwa.sh <runfile>
################################################################################################
set -x
redmine=hpcbio-redmine@igb.illinois.edu

echo -e "\n\n########################################################################################" >&2
echo -e "#############                CHECKING PARAMETERS                         ###############" >&2
echo -e "########################################################################################\n\n" >&2
set -x

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )  
scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )        
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
bwamemdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
queue=$( cat $runfile | grep -w PBSQUEUE | cut -d '=' -f2 )
pbswalltime=$( cat $runfile | grep -w PBSWALLTIME | cut -d '=' -f2 )

if [ -z $email ]
then
   MSG="Invalid value for parameter PBSEMAIL=$email in the runfile"
   echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi

set +x 
echo -e "\n\n########################################################################################" >&2
echo -e "###########                      checking PBS params                      ##################" >&2
echo -e "########################################################################################\n\n" >&2
set -x

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

echo -e "\n\n########################################################################################" >&2
echo -e "###########                      checking tools                       ##################" >&2
echo -e "########################################################################################\n\n" >&2
set -x

########################## Insert commands to check the full paths of tools :)

ash  $samtoolsdir/samtools 2>/dev/null || { echo >&2 "I require sambtools but it's not installed.  Aborting."; exit 1; }

if [ ! -d  $bwamemdir  ]
then
        MSG="Invalid value specified for BWAMEMDIR=$bwamemdir in the runfile."
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        exit 1;
fi

set +x
echo -e "\n\n########################################################################################" >&2
echo -e "#############  Everything seems ok. Now setup/configure output folders and files   #########" >&2
echo -e "########################################################################################\n\n" >&2
set -x

if [ ! -d $outputdir ]; then
	mkdir $outputdir
fi

setfacl -Rm   g::rwx $outputdir  #gives the group rwx permission, and to subdirectories
setfacl -Rm d:g::rwx $outputdir  #passes the permissions to newly created files/folders

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
truncate -s 0   $outputdir/$deliverydir/docs/Summary.Report

runfile=$outputdir/$deliverydir/docs/runfile.txt
TopOutputLogs=$outputdir/logs

truncate -s 0 $TopOutputLogs/pbs.ALIGN_sweep
truncate -s 0 $TopOutputLogs/pbs.summary_dependencies

generic_qsub_header=$TopOutputLogs/qsubGenericHeader
truncate -s 0 $generic_qsub_header
echo "#!/bin/bash" > $generic_qsub_header
echo "#PBS -q $queue" >> $generic_qsub_header
echo "#PBS -m ae" >> $generic_qsub_header
echo "#PBS -M $email" >> $generic_qsub_header
echo "#PBS -l nodes=$nodes:ppn=$thr" >> $generic_qsub_header
echo "#PBS -l walltime=${pbswalltime}" >> $generic_qsub_header

set +x
echo -e "##### let's check that it worked and that the file was created                     ####" >&2
set -x

if [ ! -s $generic_qsub_header ]
then 
    MSG="$generic_qsub_header is empty"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

set +x
echo -e "\n\n########################################################################################" >&2
echo -e "################### Documenting progress on redmine with this message ##################" >&2
echo -e "########################################################################################" >&2
echo -e "##### the first part of the Report also needs to be stored in Summary.Report      ######" >&2
echo -e "########################################################################################\n\n" >&2
set -x


MSG="Sweeping parameters for bwas started by username:$USER at: "$( echo `date` )
LOGS="Documentation about this run such as config files and results of QC tests will be placed in this folder:\n\n$outputdir/$deliverydir/docs/ \n\n"
echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
echo -e "$MSG\n\nDetails:\n\n$LOGS" >> $outputdir/$deliverydir/docs/Summary.Report

set +x
echo -e "\n\n########################################################################################" >&2
echo -e "#####                               MAIN LOOP STARTS HERE                      #########" >&2
echo -e "########################################################################################" >&2
echo -e "#####  Alignment-dedup analysis: one qsub per sample                           #########" >&2
echo -e "########################################################################################\n\n" >&2
set -x

while read sampleLine
do
    if [ `expr ${#sampleLine}` -lt 1 ]
    then
	set +x 
	echo -e "\n\n########################################################################################" >&2
	echo -e "##############                 skipping empty line        ##############################" >&2
	echo -e "########################################################################################\n\n" >&2
    else
	echo -e "\n\n########################################################################################" >&2
	echo -e "#####         Processing next line $sampleLine                                ##########" >&2
	echo -e "##### col1=sample_name col2=read1 col3=read2  including full paths            ##########" >&2
	echo -e "##### sample_name will be used for directory namas and in RG line of BAM files##########" >&2
	echo -e "########################################################################################\n\n" >&2
	set -x

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
	
	set +x
	echo -e "\n\n########################################################################################" >&2
	echo -e "###   Everything seems in order. Now creating folders where results will go  ###########" >&2
	echo -e "########################################################################################\n\n" >&2
	set -x

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
	
	set +x
	echo -e "\n\n########################################################################################" >&2               
	echo -e "####   Launching Alignment script for SAMPLE $sample with R1=$FQ_R1 R2=$FQ_R2     ##########" >&2
	echo -e "########################################################################################\n\n" >&2
	set -x

	qsub1=$TopOutputLogs/qsub.alignSweep.$sample
	cat $generic_qsub_header > $qsub1
	echo "#PBS -N alignSweep.$sample" >> $qsub1
	echo "#PBS -o $TopOutputLogs/log.alignSweep.$sample.ou" >> $qsub1
	echo "#PBS -e $TopOutputLogs/log.alignSweep.$sample.in" >> $qsub1
	echo "$scriptdir/ParameterSweep/bwa_sweep.sh $runfile $sample $FQ_R1 $FQ_R2 " >> $qsub1
	`chmod a+r $qsub1`               
	alignjobid=`qsub $qsub1` 
       	`qhold -h u $alignjobid`
	echo $alignjobid >> $TopOutputLogs/pbs.ALIGN_sweep
	echo $alignjobid >> $TopOutputLogs/pbs.summary_dependencies
	echo `date`

       	if [ `expr ${#alignjobid}` -lt 1 ]
        then
	     MSG="unable to launch qsub align job for $sample. Exiting now"
	     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
	     exit 1        
	        
        fi
	
 	set +x
	echo -e "\n\n#######  this jobid=$alignjobid will be used to hold execution of alignmentplotting.R     ########\n\n" >&2
	set -x	

   set +x
	   echo -e "\n\n########################################################################################" >&2         
	   echo -e "####   Next loop2 for Launching alignmentplotting.R script for SAMPLE $sample      #####" >&2
	   echo -e "########################################################################################\n\n" >&2
	   set -x
           alnjobid=$( echo $alignjobid | cut -d '.' -f 1 )
        
	   truncate -s 0 $TopOutputLogs/pbs.plot.$sample
	
		echo -e "\n\n########################################################################################" >&2 
		echo -e "####   Plotting script for SAMPLE $sample                               #######" >&2
		echo -e "########################################################################################\n\n" >&2
		set -x

		qsub1=$TopOutputLogs/qsub.alignPlot.$sample
		cat $generic_qsub_header > $qsub1
		echo "#PBS -N alignPlot.$sample" >> $qsub1
		echo "#PBS -o $TopOutputLogs/log.alignPlot.$sample.ou" >> $qsub1
		echo "#PBS -e $TopOutputLogs/log.alignPlot.$sample.in" >> $qsub1
                echo "#PBS -W depend=afterok:$alnjobid" >> $qsub1
################################################# azza: here should only realign/recalibrate!
		echo "$scriptdir/ParameterSweep/realign_varcall_by_chr.sh $runfile $sample $chr $TopOutputLogs/log.alignPlot.$sample.$chr.in $TopOutputLogs/log.alignPlot.$sample.$chr.ou $TopOutputLogs/qsub.alignPlot.$sample.$chr" >> $qsub1
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
           set +x  
	   echo -e "\n\n########################################################################################" >&2            
	   echo -e "####   Out of loop2. Now launching merge_vcf script for SAMPLE $sample       ##########" >&2
	   echo -e "########################################################################################\n\n" >&2
	   set -x 

           vcalljobids=$( cat $TopOutputLogs/pbs.VCALL.$sample | sed "s/\.[a-z]*//g" | tr "\n" ":" )

	   echo -e "\n\n### this list of jobids=[$vcalljobids] will be used to hold execution of merge_vcf.sh #####\n\n"

	   qsub1=$TopOutputLogs/qsub.merge.$sample
	   cat $generic_qsub_header > $qsub1
	   echo "#PBS -N merge.$sample" >> $qsub1
	   echo "#PBS -o $TopOutputLogs/log.merge.$sample.ou" >> $qsub1
	   echo "#PBS -e $TopOutputLogs/log.merge.$sample.in" >> $qsub1
           echo "#PBS -W depend=afterok:$vcalljobids" >> $qsub1
################################################### azza: here should only be merge_bams of each sample
	   echo "$scriptdir/merge_vcf.sh $runfile $sample $TopOutputLogs/log.mergeVcf.$sample.in $TopOutputLogs/log.merge.$sample.ou $TopOutputLogs/qsub.merge.$sample" >> $qsub1
	   `chmod a+r $qsub1`               
	   mergejobid=`qsub $qsub1` 
	   echo $mergejobid >> $TopOutputLogs/pbs.summary_dependencies
	   echo `date`
	
           if [ `expr ${#mergejobid}` -lt 1 ]
              then
	        MSG="unable to launch qsub merge job for $sample. Exiting now"
	        echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email" 
	        exit 1        
           fi
#azza: here should be calling variants using HC

        fi # close the if statement checking whether the workflow end with alignment or not
        # release all held jobs
        `qrls -h u $alignjobid`
   fi  # end non-empty line
done <  $sampleinfo	

############################################################################################################################# azza's GenotypeGVCF

	   set +x  
           echo -e "\n\n########################################################################################" >&2
           echo -e "####    Now launching joint_genotyping script for all SAMPLEs: each 200 together        ##########" >&2
           echo -e "########################################################################################\n\n" >&2
           set -x 

           mergedjobsids=$( cat $TopOutputLogs/pbs.summary_dependencies | sed "s/\.[a-z]*//g" | tr "\n" ":" )

           echo -e "\n\n### this list of jobids=[$mergedjobsids] will be used to hold execution of joint_vcfs.sh #####\n\n"

           qsub1=$TopOutputLogs/qsub.jointcall
           cat $generic_qsub_header > $qsub1
           echo "#PBS -N JointCalling" >> $qsub1
           echo "#PBS -o $TopOutputLogs/log.jointcall.ou" >> $qsub1
           echo "#PBS -e $TopOutputLogs/log.jointcall.in" >> $qsub1
           echo "#PBS -W depend=afterok:$mergedjobsids" >> $qsub1
           echo "$scriptdir/joint_vcf.sh $runfile $TopOutputLogs/log.jointcall.in $TopOutputLogs/log.jointcall.ou $TopOutputLogs/qsub.jointcall" >> $qsub1
           `chmod a+r $qsub1`
           jointcalljobid=`qsub $qsub1`
           echo $jointcalljobid >> $TopOutputLogs/pbs.summary_dependencies
           echo `date`

           if [ `expr ${#jointcalljobid}` -lt 1 ]
              then
                MSG="unable to launch qsub jointVCFcall job for 200 samples. Exiting now"
                echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1
           fi

################################################################################################################################### end azza's block

set +x
echo -e "\n\n########################################################################################" >&2
echo -e "########################################################################################" >&2
echo -e "#################           MAIN LOOP ENDS HERE                  #######################" >&2
echo -e "########################################################################################" >&2
echo -e "########################################################################################" >&2
echo -e "#################     Now, we need to generate summary           #######################" >&2
echo -e "########################################################################################" >&2
echo -e "########################################################################################\n\n" >&2
set -x

alljobids=$( cat $TopOutputLogs/pbs.summary_dependencies | sed "s/\.[a-z]*//g" | tr "\n" ":" )

set +x
echo -e "\n\n### this list of jobids=[$alljobids] will be used to hold execution of summary.sh #####\n\n" >&2
set -x

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

 # release all held jobs
 `qrls -h u $alignjobid`

set +x        
echo -e "\n\n########################################################################################" >&2
echo -e "##############                 EXITING NOW                            ##################" >&2	
echo -e "########################################################################################\n\n" >&2
set -x
