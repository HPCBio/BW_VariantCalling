#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# -ne 9 ]
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

    set -x
    echo `date`
    umask 0027
    scriptfile=$0
    runfile=$1
    sample=$2
    RealignOutput=$3
    outfile=$4
    rootdir=$5
    elog=$6
    olog=$7
    email=$8
    qsubfile=$9
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

    memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
    verifydir=$( cat $runfile | grep -w VERIFYBAMDIR | cut -d '=' -f2 )    
    freemix_cutoff=$( cat $runfile | grep -w FREEMIX_CUTOFF | cut -d '=' -f2 )
    refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
    sites=$refdir/omni.sites.vcf.gz
    deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
    threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
    thr=`expr $threads "-" 1`
 
    if [ ! -d $rootdir ]
    then
       MSG="$rootdir output root directory not found. verifySample failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ `expr ${#deliveryfolder}` -lt 2 ]
    then
        delivery=$rootdir/delivery
    else
        delivery=$rootdir/$deliveryfolder
    fi
    if [ ! -d $delivery/Cleaned_BAMS ]
    then
        mkdir -p $delivery/Cleaned_BAMS
        #mkdir -p $delivery/QC_results
    fi

    if [ ! -d $RealignOutput ]
    then
       MSG="$RealignOutput  realign directory not found. verifySample failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -s $sites ]
    then
       MSG="$sites genotyping info file not found. verifySample failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    set +x; echo -e "\n\n#############    resetting files and logs  ###############\n\n" >&2; set -x;
    
    goodSamples=$rootdir/VERIFIED_ok_SAMPLES.list
    badSamples=$rootdir/VERIFIED_notok_SAMPLES.list
    mergedfile=${sample}.realignedSample.calmd.bam
    outputfile=$delivery/Cleaned_BAMS/${sample}.Improved.realignedSample.bam
    qc_result=$rootdir/QC_Results.txt

    if [ ! -s $goodSamples ]
    then
        truncate -s 0 $goodSamples
    fi
    if [ ! -s $badSamples ]
    then
        truncate -s 0 $badSamples
    fi
    if [ ! -s $rootdir/QC_Results.txt ]
    then
        truncate -s 0 $rootdir/QC_Results.txt
    fi
    #if [ $freemix_cutoff -eq $freemix_cutoff 2>/dev/null ]
    #then
    #    echo -e "ok val, it is numeric"
    #else
	#MSG="FREEMIX_CUTFF=$freemix_cutff invalid value. verifySample failed"
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#exit $exitcode;
    #fi


    set +x; echo -e "\n\n#############   merging all chr for calmd.bam corresponding to sample=$sample  ###############\n\n" >&2; set -x;
    
    `echo date`

    cd $RealignOutput

    chrbamfiles=`find ./ -name "*${sample}*.calmd.bam"`
    #bamsList=$( echo $chrbamfiles | sed "s/\.\// INPUT=/g" | tr "\n" " " )
    bamsList=$( echo $chrbamfiles | sed "s/\.\// /g" | tr "\n" " " )
    if [ `expr ${#bamsList}` -lt 2 ]
    then
            # empty list nothing to merge
	MSG="$RealignOutput no bam files to merge for sample: $sample. verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit 1;
    fi
    `echo date`

    if [ ! -s $mergedfile ]
    then
            # we did not perform this step already
        #$javadir/java -Xmx1024m -Xms1024m -jar $picardir/MergeSamFiles.jar $bamsList OUTPUT=$mergedfile USE_THREADING=true 
        $novodir/novosort --index --tmpdir $RealignOutput --threads $thr -m 16g --kt -o $mergedfile $bamsList
        exitcode=$?
        `echo date`
        if [ $exitcode -ne 0 ]
        then
	    MSG="novosort command failed exitcode=$exitcode verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
    fi

    set +x; echo -e "\n\n#############  checking that file was created  ###############\n\n" >&2; set -x;

    if [ ! -s $mergedfile ]
    then
	MSG="novosort command produced an empty file while trying to merge $bamsList verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
    fi


    set +x; echo -e "generating an index file for   $mergedfile because verifyBamID needs it" >&2; set -x;


    #$samdir/samtools index $mergedfile 
    #`echo date`
    if [ ! -s $mergedfile.bai ]
    then
	MSG="novosort did not produced index file for $mergedfile  verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
    fi

    set +x; echo -e "\n\n" >&2; 
    echo -e "##########################################################################" >&2
    echo -e "##########################################################################" >&2
    echo -e "#############  populating delivery folder                #################" >&2
    echo -e "#############  with improved bam folr sample=$sample     #################" >&2
    echo -e "##########################################################################" >&2
    echo -e "\n\n" >&2; set -x;

    `echo date`
    cp $mergedfile $outputfile
    cp $mergedfile.bai $outputfile.bai

    if [ ! -s $outputfile ]
    then
	MSG="WARNING=$outputfile output file for sample=$sample was not copied properly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    fi

    set +x; echo -e "\n\n" >&2; 
    echo -e "##########################################################################" >&2
    echo -e "##########################################################################" >&2
    echo -e "#############  running verifyBamID for sample=$sample    #################" >&2
    echo -e "##########################################################################" >&2
    echo -e "##########################################################################" >&2
    echo -e "\n\n" >&2; set -x;
    
    `echo date`


    $verifydir/verifyBamID --vcf $sites --bam $mergedfile --out $outfile --verbose --ignoreRB --chip-none 
    exitcode=$?
    `echo date`
    if [ $exitcode -ne 0 ]
    then
	MSG="verifyBamID command failed exitcode=$exitcode verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    if [ ! -s $outfile* ]
    then
	MSG="verifyBamID did not produced output files  exitcode=$exitcode verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    set +x; echo -e "\n\n" >&2; 
    echo -e "##########################################################################" >&2
    echo -e "#############  filtering for sample=$sample              #################" >&2
    echo -e "               parsing verifybam output file" >&2
    echo -e "##########################################################################" >&2
    echo -e "\n\n" >&2; set -x;
    
    
    verifyfile=`find ./ -name "*.selfSM"`

    freemix=$( cat $verifyfile | sed '2q;d' | cut -f 7 ) 

    set +x; echo -e "\n\n" >&2; 
    echo -e "##########################################################################" >&2
    echo -e "checking that we actually grabbed something that is a number" >&2
    echo -e "using a weird trick to do the comparison in bash between noninteger variables" >&2
    echo -e "##########################################################################" >&2
    echo -e "\n\n" >&2; set -x;    
 
    if [ $(echo "$freemix < $freemix_cutoff"|bc) -eq 1 ]
    then
        set +x; echo -e "\n\n" >&2;    
	echo -e "##########################################################################" >&2
        echo -e "sample-$sample PASSED verifyBamID filter\nadd to list of good-samples" >&2
	echo -e "##########################################################################" >&2
        echo -e "\n\n" >&2; set -x;	
	
        filterflag="PASSED"	
        detail=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
        echo -e "$detail" >> $goodSamples
        echo -e "$detail" >> $qc_result
    else
    
        set +x; echo -e "\n\n" >&2;  
	echo -e "##########################################################################" >&2
        echo -e "sample-$sample FAILED verifyBamID filter\nadd to list of bad-samples" >&2
	echo -e "##########################################################################" >&2
        echo -e "\n\n" >&2; set -x;	

        filterflag="FAILED"
        detail=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
        echo -e "$detail" >> $badSamples
        echo -e "$detail" >> $qc_result
    fi
    echo `date`

    if [ $filterflag == "FAILED" ]
    then

        set +x; echo -e "\n\n" >&2;  
	echo -e "##########################################################################" >&2
        echo -e "############# $sample FAILED the verifyBamId filter" >&2
	echo -e "##########################################################################" >&2
        echo -e "\n\n" >&2; set -x;	
        

	MSG="$sample DID NOT PASS the VerifyBamID filter. Stopping analysis now"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	#exit 1;
    else
        set +x; echo -e "\n\n" >&2;  
	echo -e "##########################################################################" >&2
        echo -e "############# $sample PASSED the verifyBamId filter" >&2
	echo -e "##########################################################################" >&2
        echo -e "\n\n" >&2; set -x;	
        
        exit 0;
    fi
