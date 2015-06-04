#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#redmine=hpcbio-redmine@igb.illinois.edu
if [ $# -ne 11 ]
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
    RealignOutputLogs=${11}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
    javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
    dup_cutoff=$( cat $runfile | grep -w DUP_CUTOFF | cut -d '=' -f2 )
    map_cutoff=$( cat $runfile | grep -w MAP_CUTOFF | cut -d '=' -f2 )
    chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
    indices=$( echo $chrindex | sed 's/:/ /g' )
    sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
    sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
    sPU=$( cat $runfile | grep -w SAMPLEPU | cut -d '=' -f2 )
    sID=$( echo $sample )
    sSM=$( echo $sample )
    parameters="RGID=${sID} RGSM=${sSM} RGLB=${sLB} RGPL=${sPL} RGPU=${sPU} RGCN=${sCN}"
    RealignOutput=$rootdir/$sample/realign


    if [ ! -d $rootdir ]
    then
       MSG="$rootdir output directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
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
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="$picardir samtools directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    if [ -z $javadir ]
    then
	MSG="A value must be specified for JAVADIR in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    #else
        #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    #        `module load $javamodule`
    fi          


    #now merging all chr*.dedu-real-recal.bams for a lane into a single file
    #in order to calculate stats
    good2mergeLanes=""
    bad2mergeLanes=""
    lanefiles=$( echo $lanefiles | sed 's/::/ /g' | sed 's/:/ /g' )
    for lane in $lanefiles
    do
        if [ `expr ${#lane}` -gt 2 ]
        then
        #now cheching that we have data to work with
        inputdir=$rootdir/${sample}/${lane}/realign
        if [ ! -d $inputdir ]
        then
	    MSG="$inputdir directory not found for sample: $sample lane: $lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        `echo date`
        cd $inputdir

        #now creating temp files per lane
        mergedLane=${RealignOutput}/${lane}.merged.bam
        #bamsList=${lane}.bams.list
        prestats=${RealignOutput}/${lane}.merged.bam.flagstat
        mergedLaneMappedPaired=${RealignOutput}/${lane}.merged.mappedPaired.bam
        poststats=${RealignOutput}/${lane}.merged.mappedPaired.bam.flasgstat

        #truncate -s 0 $bamsList
        bamsList=`find ./ -name "*.${lane}.output.bam"`
        bamsList=$( echo $bamsList | sed "s/\.\// INPUT=/g" | tr "\n" " " )
        if [ `expr ${#bamsList}` -lt 1 ]
        then
            # empty list nothing to merge
	    MSG="In dir:$inputdir no bam files found for lane: $lane. Nothing to merge"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        `echo date`

        #now merging all chr bam files with picard
        #remove this if later on after done testing

        if [ ! -s $mergedLane ]
        then
            # we did not perform this step already
            $javadir/java -Xmx1024m -Xms1024m -jar $picardir/MergeSamFiles.jar $bamsList OUTPUT=$mergedLane USE_THREADING=true 
            exitcode=$?
            `echo date`
            if [ $exitcode -ne 0 ]
            then
		MSG="picard MergeSamFiles command failed exitcode=$exitcode while trying to merge $bamsList"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
        fi

        if [ ! -s $mergedLane ]
        then
	    MSG="picard MergeSamFiles command produced an empty file while trying to merge $bamsList"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi



        #now applying filters and deciding whether to discard or keep this lane

        echo -e "filtering begins"
        echo -e "producing properly-mapped bam file and then flagstat on both, prefilter and post filter"

        # remove this if later on after done with testing
        if [ ! -s $prestats ]
        then
            $samdir/samtools flagstat $mergedLane > $prestats
            exitcode=$?
            `echo date`
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools flagstat command failed exitcode=$exitcode with $mergedLane"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
	fi
        if [ ! -s $prestats ]
        then
	    MSG="samtools flagstat command produced an empty file with $mergedLane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

        # remove this if later on after done with testing
        if [ ! -s $mergedLaneMappedPaired ]
        then
            $samdir/samtools view -b -F 12 $mergedLane > $mergedLaneMappedPaired 
            exitcode=$?
            `echo date`
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools flagstat command failed exitcode=$exitcode with $mergedLaneMappedPaired"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
        fi

        if [ ! -s $mergedLaneMappedPaired ]
        then
	    MSG="samtools flagstat command produced an empty file with $mergedLaneMappedPaired"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

        # remove this if later on after done with testing
        if [ ! -s $poststats ]
        then
            $samdir/samtools flagstat $mergedLaneMappedPaired > $poststats 
            exitcode=$?
            `echo date`
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools flagstat command failed exitcode=$exitcode with $mergedLaneMappedPaired"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
	fi

        if [ ! -s $poststats ]
        then
	    MSG="samtools flagstat command produced an empty file with $mergedLane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi


        # now extracting the stats of interest from the flagstats

        tot_mapped=$( cat $poststats | grep "mapped (" | cut -d ' ' -f1 )
        tot_reads=$( cat $prestats | grep "in total" | cut -d ' ' -f1 )
        tot_dups=$( cat $poststats | grep "duplicates" | cut -d ' ' -f1 )

        #now testing if these variables have numbers

        if [ $tot_dups -eq $tot_dups 2>/dev/null ]
        then
           echo -e "ok val"
        else
	    MSG="$prestats samtools flagstat file parsed incorrectly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
        if [ $tot_reads -eq $tot_reads 2>/dev/null ]
        then
           echo -e "ok val"
        else
	    MSG="$prestats samtools flagstat file parsed incorrectly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi

        if [ $tot_mapped -eq $tot_mapped 2>/dev/null ]
        then
           echo -e "ok val"
        else
	    MSG="$prestats samtools flagstat file parsed incorrectly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

        # now calculating percentages to be used for the filtering step

        perc_dup=$(( tot_dups * 100 / tot_reads ))
        perc_mapped=$(( tot_mapped * 100 / tot_reads ))

        if [ $perc_dup -eq $perc_dup 2>/dev/null ]
        then
           echo -e "ok val"
        else
	    MSG="$stats samtools flagstat file parsed incorrectly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

        if [ $perc_mapped -eq $perc_mapped 2>/dev/null ]
        then
           echo -e "ok val"
        else
	    MSG="$stats samtools flagstat file parsed incorrectly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

        echo -e "now applying the filters"

        if [ $perc_dup -lt $dup_cutoff ]
        then
           echo -e "$lane passed filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
           if [ $perc_mapped -gt $map_cutoff ]
           then
              echo -e "$lane passed filter percent_mapped with value $perc_mapped, minimum cutoff is $map_cutoff"
              echo -e "adding $lane to list of accepted bams to merge per sample"
	      good2mergeLanes="INPUT=$mergedLane "$good2mergeLanes
	      #good2mergeLanes="INPUT=$mergedLaneMappedPaired  "$good2mergeLanes

           else
              echo -e "$lane DID NOT pass filter percent_mapped value $perc_mapped, minimum cutoff is $map_cutoff"
              echo -e "adding $lane to list of rejected bams to merge per sample"
	      bad2mergeLanes="INPUT=$mergedLane "$bad2mergeLanes
           fi
        else
            echo -e "$lane DID NOT pass filter percent_duplicates value $perc_dup, minimum cutoff is $dup_cutoff"
            echo -e "adding $lane to list of rejected bams to merge per sample"
	    bad2mergeLanes="INPUT=$mergedLane "$bad2mergeLanes
        fi
        fi # skip empty lines
        # lane processed, next one please
        `echo date`
    done
    echo `date`
    echo -e "now we need to make sure that there are goodlanes to work with"
    if [ `expr ${#good2mergeLanes}` -lt 1 ]
    then
        echo -e "no good lanes to merge"
	MSG="no good lanes to merge"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
        exit 1;
    fi


    echo -e "now merging lanes that passed BOTH filters= $good2mergeLanes"
 
    cd $RealignOutput
    presorted=$RealignOutput/presorted.$outfile
    sorted_norg=$RealignOutput/sorted.norg.$outfile

    #remove this if after we are done with testing
    if [ ! -s $presorted ]
    then
	$javadir/java -Xmx1024m -Xms1024m -jar $picardir/MergeSamFiles.jar \
            $good2mergeLanes \
            OUTPUT=$presorted \
            USE_THREADING=true 

	exitcode=$?
	`echo date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="picard MergeSamFiles command failed exitcode=$exitcode while trying to merge $good2mergeLanes"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
    fi
    if [ ! -s $presorted ]
    then
	MSG="picard MergeSamFiles command produced an empty file while trying to merge $good2mergeLanes"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi

    $javadir/java -Xmx1024m -Xms1024m -jar $picardir/SortSam.jar \
        INPUT=$presorted \
        OUTPUT=${sorted_norg}\
        TMP_DIR=$RealignOutput \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true 
 
    exitcode=$?
    `echo date`
    if [ $exitcode -ne 0 ]
    then
	MSG="picard SortSam command failed exitcode=$exitcode with $presorted"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    if [ ! -s ${sorted_norg} ]
    then
	MSG="picard SortSam command produced an empty file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi

    $javadir/java -Xmx1024m -Xms1024m -jar $picardir/AddOrReplaceReadGroups.jar \
        INPUT=${sorted_norg} \
        OUTPUT=$RealignOutput/$outfile \
        TMP_DIR=$RealignOutput \
        SORT_ORDER=coordinate \
        $parameters
 
    exitcode=$?
    `echo date`
    if [ $exitcode -ne 0 ]
    then
	MSG="picard addreadgroup command failed exitcode=$exitcode with $outfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    if [ ! -s $outfile ]
    then
	MSG="picard addreadgroup command produced an empty file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi

    $samdir/samtools index $outfile 
    exitcode=$?
    `echo date`
    if [ $exitcode -ne 0 ]
    then
	MSG="samtools index command failed exitcode=$exitcode mergedbylane $outfile"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi

    echo -e "creating text files one with the list of good and one with the list of discarded"
    echo -e "bams per lane used for merging and crating a single bam by sample"
    
    truncate -s o $RealignOutput/goodLanes4MergingBySample.list
    truncate -s o $RealignOutput/badLanes4MergingBySample.list

    echo $good2mergeLanes > $RealignOutput/goodLanes4MergingBySample.list
    echo $bad2mergeLanes  > $RealignOutput/badLanes4MergingBySample.list
    `echo date`
    # exiting now
fi
