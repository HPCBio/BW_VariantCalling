#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# -ne 10 ]
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
else
    set -x
    echo `date`
    scriptfile=$0
    runfile=$1
    sample=$2
    sLB=$3
    lanefiles=$4
    outfile=$5
    rootdir=$6
    elog=$7
    olog=$8
    email=$9
    qsubfile=${10}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
    
   echo "####################################################################################################"
   echo "#####################################################################################################"
   echo "#####################################                       ########################################"
   echo "##################################### PARSING RUN INFO FILE ########################################"
   echo "#####################################  SANITY CHECK         ########################################"
   echo "####################################################################################################" 
   
    memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
    javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
    chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
    dup_cutoff=$( cat $runfile | grep -w DUP_CUTOFF | cut -d '=' -f2 )
    map_cutoff=$( cat $runfile | grep -w MAP_CUTOFF | cut -d '=' -f2 )
    chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
    indices=$( echo $chrindex | sed 's/:/ /g' )
    thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
    sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
    sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
    sPU=$( cat $runfile | grep -w SAMPLEPU | cut -d '=' -f2 )
    sID=$( echo $sample )
    sSM=$( echo $sample )

    parameters="RGID=${sID} RGSM=${sSM} RGLB=${sLB} RGPL=${sPL} RGPU=${sPU} RGCN=${sCN}"
    RGline=$( echo -e "@RG\tID:${sID}\tSM:${sSM}\tLB:${sLB}\tPL:${sPL}\tPU:${sPU}\tCN:${sCN}" ) 
    RealignOutput=$rootdir/$sample/realign


    if [ `expr ${#sCN}` -lt 2 -o `expr ${#sPL}` -lt 2 -o `expr ${#sPU}` -lt 2 -o `expr ${#sID}` -lt 2 -o `expr ${#sSM}` -lt 2 ]
    then
       MSG="some required fields to construct RG line are missing"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    if [ ! -d $rootdir ]
    then
       MSG="$rootdir output directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $RealignOutput ]
    then
       echo -e "creating $RealignOutput"
       mkdir -p $RealignOutput
    fi
    if [ ! -d $picardir ]
    then
       MSG="$picardir picard directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="$picardir samtools directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    if [ -z $javadir ]
    then
	MSG="A value must be specified for JAVADIR in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit 1;
    #else
        #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    #        `module load $javamodule`
    fi          

    echo "####################################################################################################"
    echo "#####################################################################################################"
    echo "           creating text files with QC results for EACH LANE"
    echo "####################################################################################################"
    echo "#####################################################################################################"
   
    truncate -s 0 $rootdir/QC_Results.${sample}.txt

    good2mergeLanes=""
    bad2mergeLanes=""

   echo "####################################################################################################"
   echo "#####################################################################################################"
   echo "            MAIN LOOP BEGINS HERE "
   echo "####################################################################################################"
   echo "#####################################################################################################"

    lanefiles=$( echo $lanefiles | sed 's/::/ /g' | sed 's/:/ /g' )
    for lane in $lanefiles
    do
        if [ `expr ${#lane}` -gt 2 ]
        then
	   echo "#####################################################################################################"
           echo "                   now cheching that we have data to work with"
           echo "####################################################################################################"
          
           inputdir=$rootdir/${sample}/${lane}/realign
           if [ ! -d $inputdir ]
           then
	       MSG="$inputdir directory not found for sample: $sample lane: $lane"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit 1;
           fi
           `echo date`
           cd $inputdir

	   echo "#####################################################################################################"
           echo "              now creating temp files per lane"
	   echo "#####################################################################################################"
           
           mergedLane=${RealignOutput}/${lane}.merged.bam
           #bamsList=${lane}.bams.list
           prestats=${RealignOutput}/${lane}.merged.bam.flagstat
           mergedLaneMappedPaired=${RealignOutput}/${lane}.merged.mappedPaired.bam
           poststats=${RealignOutput}/${lane}.merged.mappedPaired.bam.flasgstat

           #truncate -s 0 $bamsList
           bamsList=`find ./ -name "*.${lane}.output.bam"`
           bamsListPlain=$( echo $bamsList | sed "s/\.\// /g" | tr "\n" " " )
           bamsList=$( echo $bamsList | sed "s/\.\// INPUT=/g" | tr "\n" " " )
           if [ `expr ${#bamsList}` -lt 2 ]
           then
 	        echo "#####################################################################################################"
                echo "           empty list. nothing to merge. Stopping analysis for sample=$sample"
 	        echo "#####################################################################################################"
               
	        MSG="In dir:$inputdir no bam files found for lane: $lane. Nothing to merge"
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	        #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	        exit 1;
           fi
           `echo date`

 	   echo "#####################################################################################################"
           echo "         FILE1: bam file with all reads for lane=$lane sample=$sample"
           echo "         merging all chr bam files with picard and producing stats with samtools flagstat"
 	   echo "#####################################################################################################"


           if [ ! -s $mergedLane ]
           then
                #$javadir/java -Xmx8g -Xms1024m -jar $picardir/MergeSamFiles.jar $bamsList OUTPUT=$mergedLane USE_THREADING=true 
                $novodir/novosort --index --threads $thr --tmpdir $inputdir -m 16g --kt --compression 1 -o $mergedLane $bamsListPlain
                exitcode=$?
               `echo date`
               if [ $exitcode -ne 0 ]
               then
		   MSG="novosort command failed exitcode=$exitcode while trying to merge $bamsList"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		   exit $exitcode;
               fi
           fi

           if [ ! -s $mergedLane ]
           then
	       MSG="novosort command produced an empty file while trying to merge $bamsList"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi

           if [ ! -s $prestats ]
           then
               $samdir/samtools flagstat $mergedLane > $prestats
               exitcode=$?
               `echo date`
               if [ $exitcode -ne 0 ]
               then
		   MSG="samtools flagstat command failed exitcode=$exitcode with $mergedLane"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		   exit $exitcode;
               fi
	   fi
           if [ ! -s $prestats ]
           then
	       MSG="samtools flagstat command produced an empty file with $mergedLane"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi


 	   echo "#####################################################################################################"
           echo "         FILE2: bam file with ONLY properly mapped reads for lane=$lane sample=$sample"
           echo "         and producing stats with samtools flagstat"
 	   echo "#####################################################################################################"

           if [ ! -s $mergedLaneMappedPaired ]
           then
              $samdir/samtools view -bu -F 12 $mergedLane > $mergedLaneMappedPaired 
              exitcode=$?
              `echo date`
              if [ $exitcode -ne 0 ]
              then
		  MSG="samtools flagstat command failed exitcode=$exitcode with $mergedLaneMappedPaired"
		  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	          #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		  exit $exitcode;
              fi
           fi

           if [ ! -s $mergedLaneMappedPaired ]
           then
	       MSG="samtools flagstat command produced an empty file with $mergedLaneMappedPaired"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi


           if [ ! -s $poststats ]
           then
               $samdir/samtools flagstat $mergedLaneMappedPaired > $poststats 
               exitcode=$?
               `echo date`
               if [ $exitcode -ne 0 ]
               then
		   MSG="samtools flagstat command failed exitcode=$exitcode with $mergedLaneMappedPaired"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		   exit $exitcode;
               fi
	   fi

           if [ ! -s $poststats ]
           then
	       MSG="samtools flagstat command produced an empty file with $mergedLane"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi

 	   echo "#####################################################################################################"
 	   echo "#####################################################################################################"
           echo "             now extracting the stats of interest from the samtools flagstat cmd "
 	   echo "#####################################################################################################"
 	   echo "#####################################################################################################"
           

           tot_mapped=$( cat $poststats | grep "mapped (" | cut -d ' ' -f1 )
           tot_reads=$( cat $prestats | grep "in total" | cut -d ' ' -f1 )
           tot_dups=$( cat $poststats | grep "duplicates" | cut -d ' ' -f1 )

           #now testing if these variables have integer numbers

           if [ $tot_dups -eq $tot_dups 2>/dev/null ]
           then
               echo -e "ok val"
           else
	       MSG="$prestats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
	   fi
           if [ $tot_reads -eq $tot_reads 2>/dev/null ]
           then
               echo -e "ok val"
           else
	       MSG="$prestats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
	   fi

           if [ $tot_mapped -eq $tot_mapped 2>/dev/null ]
           then
               echo -e "ok val"
           else
	       MSG="$prestats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi

 	   echo "#####################################################################################################"
 	   echo "#####################################################################################################"
           echo "      now calculating percentages to be used for the filtering step"
 	   echo "#####################################################################################################"
 	   echo "#####################################################################################################"
           

           perc_dup=$(( tot_dups * 100 / tot_reads ))
           perc_mapped=$(( tot_mapped * 100 / tot_reads ))

           if [ $perc_dup -eq $perc_dup 2>/dev/null ]
           then
               echo -e "ok val"
           else
	       MSG="$stats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi

           if [ $perc_mapped -eq $perc_mapped 2>/dev/null ]
           then
               echo -e "ok val"
           else
	       MSG="$stats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
           fi

 	   echo "#####################################################################################################"
 	   echo "#####################################################################################################"
           echo "                     now applying the filters and populating the output reports  "
 	   echo "#####################################################################################################"
 	   echo "#####################################################################################################"

           
           if [ $perc_dup -lt $dup_cutoff ]
           then
               echo -e "$lane passed filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
               if [ $perc_mapped -gt $map_cutoff ]
               then
                  echo -e "$lane passed filter percent_mapped with value $perc_mapped, minimum cutoff is $map_cutoff"
                  echo -e "adding $lane to list of accepted bams to merge per sample"
	          good2mergeLanes="INPUT=$mergedLane "$good2mergeLanes
                  detail="$lane\tPASSED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	          echo -e "$detail" >> $rootdir/QC_Results.${sample}.txt
               else
                  echo -e "$lane DID NOT pass filter percent_mapped value $perc_mapped, minimum cutoff is $map_cutoff"
                  echo -e "adding $lane to list of rejected bams to merge per sample"
                  detail="$lane\tFAILED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	          echo -e "$detail"  >> $rootdir/QC_Results.${sample}.txt
               fi
           else
               echo -e "$lane DID NOT pass filter percent_duplicates value $perc_dup, minimum cutoff is $dup_cutoff"
               echo -e "adding $lane to list of rejected bams to merge per sample"
               detail="$lane\tFAILED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	       echo -e "$detail"  >> $rootdir/QC_Results.${sample}.txt
           fi
        fi # skip empty lines
        # lane processed, next one please
        `echo date`
    done

    
    `echo date`
    echo "#####################################################################################################"
    echo "#####################################################################################################"
    echo "       now we need to make sure that there are goodlanes to work with"
    echo "#####################################################################################################"
    echo "#####################################################################################################"
    
    if [ `expr ${#good2mergeLanes}` -lt 2 ]
    then
        echo "#####################################################################################################"
        echo "  no good lanes to merge. Stopping the analysis"
        echo "#####################################################################################################"
       
	MSG="no good lanes to merge"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
        exit 1;
    fi

    echo "#####################################################################################################"
    echo "#####################################################################################################"
    echo "           Now merging lanes for sample=$sample that passed BOTH filters= $good2mergeLanes"
    echo "           First the merge, then add readgroup, then sort and finaly index the bam file    "
    echo "#####################################################################################################"
    echo "#####################################################################################################"
 
    cd $RealignOutput

    presorted=${outfile}.presorted
    withRG=${outfile}.presorted_wrg
    bamsToMerge=$( echo $good2mergeLanes | sed "s/INPUT=/ /g" )

    #if [ ! -s $presorted ]
    #then
	#$javadir/java -Xmx8g -Xms1024m -jar $picardir/MergeSamFiles.jar \
        #    $good2mergeLanes \
        #    OUTPUT=$presorted \
        #    USE_THREADING=true 

	#exitcode=$?
	#`echo date`
	#if [ $exitcode -ne 0 ]
	#then
	#    MSG="picard MergeSamFiles command failed exitcode=$exitcode while trying to merge $good2mergeLanes"
	#    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#    exit $exitcode;
	#fi
    #fi
    #if [ ! -s $presorted ]
    #then
	#MSG="picard MergeSamFiles command produced an empty file while trying to merge $good2mergeLanes"
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#exit $exitcode;
    #fi

    #$javadir/java -Xmx8g -Xms1024m -jar $picardir/AddOrReplaceReadGroups.jar \
    #    INPUT=${presorted} \
    #    OUTPUT=${withRG} \
    #    TMP_DIR=$RealignOutput \
    #    SORT_ORDER=coordinate \
    #    $parameters
 
    #exitcode=$?
    #`echo date`
    #if [ $exitcode -ne 0 ]
    #then
	#MSG="picard addreadgroup command failed exitcode=$exitcode with $outfile"
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#exit $exitcode;
    #fi
    #if [ ! -s $withRG ]
    #then
	#MSG="picard addreadgroup command produced an empty file"
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#exit $exitcode;
    #fi


    #$javadir/java -Xmx8g -Xms1024m -jar $picardir/SortSam.jar \
    #    INPUT=$withRG \
    #    OUTPUT=$outfile \
    #    TMP_DIR=$RealignOutput \
    #    SORT_ORDER=coordinate \
    #    CREATE_INDEX=true 
 
    
    $novodir/novosort --index --threads $thr --tmpdir $RealignOutput --rg "${RGline}"  -m 16g --kt --compression 1 -o $outfile $bamsToMerge 
    exitcode=$?
    `echo date`
    if [ $exitcode -ne 0 ]
    then
	MSG="novosort command failed exitcode=$exitcode with $presorted"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    if [ ! -s $outfile ]
    then
	MSG="novosort command produced an empty file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi


    #$samdir/samtools index $outfile 
    #exitcode=$?
    #`echo date`
    #if [ $exitcode -ne 0 ]
    #then
	#MSG="samtools index command failed exitcode=$exitcode mergedbylane $outfile"
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#exit $exitcode;
    #fi

    # exiting now
fi
