#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#redmine=hpcbio-redmine@igb.illinois.edu
if [ $# -ne 9 ]
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

    set -x
    echo `date`
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

    echo -e "##########################################################################"
    echo -e "################  parse runfile and sanity check #########################"
    echo -e "##########################################################################"

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
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ `expr ${#deliveryfolder}` -lt 2 ]
    then
        delivery=$rootdir/delivery
    else
        delivery=$rootdir/$deliveryfolder
    fi
    if [ ! -d $delivery ]
    then
        mkdir -p $delivery/Cleaned_BAMS
        mkdir -p $delivery/QC_results
    fi

    if [ ! -d $RealignOutput ]
    then
       MSG="$RealignOutput  realign directory not found. verifySample failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -s $sites ]
    then
       MSG="$sites genotyping info file not found. verifySample failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    goodSamples=$rootdir/VERIFIED_ok_SAMPLES.list
    badSamples=$rootdir/VERIFIED_notok_SAMPLES.list
    mergedfile=${sample}.realignedSample.calmd.bam
    outputfile=$delivery/Cleaned_BAMS/${sample}.Improved.realignedSample.bam
    qc_result=$rootdir/QC_Results.${sample}.txt
    qc_result2=$RealignOutput/QC_Results.${sample}.txt
    qc_result3=$delivery/QC_results/QC_Results.${sample}.txt
    if [ ! -s $goodSamples ]
    then
        truncate -s 0 $goodSamples
    fi
    if [ ! -s $badSamples ]
    then
        truncate -s 0 $badSamples
    fi
    if [ ! -s $rootdir/QC_Results.${sample}.txt ]
    then
        truncate -s 0 $rootdir/QC_Results.${sample}.txt
    fi
    #if [ $freemix_cutoff -eq $freemix_cutoff 2>/dev/null ]
    #then
    #    echo -e "ok val, it is numeric"
    #else
	#MSG="FREEMIX_CUTFF=$freemix_cutff invalid value. verifySample failed"
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	#exit $exitcode;
    #fi

    echo -e "##########################################################################"
    echo -e "##########################################################################"
    echo -e "#############  merging all chr files for sample=$sample  #################"
    echo -e "##########################################################################"
    echo -e "##########################################################################"
    `echo date`

    cd $RealignOutput

    chrbamfiles=`find ./ -name "*${sample}.*.calmd.bam"`
    #bamsList=$( echo $chrbamfiles | sed "s/\.\// INPUT=/g" | tr "\n" " " )
    bamsList=$( echo $chrbamfiles | sed "s/\.\// /g" | tr "\n" " " )
    if [ `expr ${#bamsList}` -lt 2 ]
    then
            # empty list nothing to merge
	MSG="$RealignOutput no bam files to merge for sample: $sample. verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
    fi

    echo -e "making sure  $mergedfile was created"

    if [ ! -s $mergedfile ]
    then
	MSG="novosort command produced an empty file while trying to merge $bamsList verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
    fi

    echo -e "generating an index file for   $mergedfile because verifyBamID needs it"


    #$samdir/samtools index $mergedfile 
    #`echo date`
    if [ ! -s $mergedfile.bai ]
    then
	MSG="novosort did not produced index file for $mergedfile  verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
    fi

    echo -e "##########################################################################"
    echo -e "##########################################################################"
    echo -e "#############  populating delivery folder                #################"
    echo -e "#############  with improved bam folr sample=$sample     #################"
    echo -e "##########################################################################"
    `echo date`
    cp $mergedfile $outputfile

    if [ ! -s $outputfile ]
    then
	MSG="WARNING=$outputfile output file for sample=$sample was not copied properly"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
    fi

    echo -e "##########################################################################"
    echo -e "##########################################################################"
    echo -e "#############  running verifyBamID for sample=$sample    #################"
    echo -e "##########################################################################"
    echo -e "##########################################################################"

    `echo date`


    $verifydir/verifyBamID --vcf $sites --bam $mergedfile --out $outfile --verbose --ignoreRB --chip-none 
    exitcode=$?
    `echo date`
    if [ $exitcode -ne 0 ]
    then
	MSG="verifyBamID command failed exitcode=$exitcode verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    if [ ! -s $outfile* ]
    then
	MSG="verifyBamID did not produced output files  exitcode=$exitcode verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    echo -e "##########################################################################"
    echo -e "#############  filtering for sample=$sample              #################"
    echo -e "               parsing verifybam output file"
    echo -e "##########################################################################"
   
    verifyfile=`find ./ -name "*.selfSM"`

    freemix=$( cat $verifyfile | sed '2q;d' | cut -f 7 ) 

    echo -e "##########################################################################"
    echo -e "checking that we actually grabbed something that is a number"
    echo -e "using a weird trick to do the comparison in bash between noninteger variables"
    echo -e "##########################################################################"
    
 
    if [ $(echo "$freemix < $freemix_cutoff"|bc) -eq 1 ]
    then
	echo -e "##########################################################################"
        echo -e "sample-$sample passed verifyBamID filter\nadd to list of good-samples"
	echo -e "##########################################################################"
        filterflag="PASSED"	
        detail=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
        echo -e "$detail" >> $goodSamples
        sed -i "1i $detail"  $qc_result
    else
	echo -e "############################################################################"
        echo -e "sample-$sample DID NOT passed verifyBamID filter\nadd to list of bad-samples"
	echo -e "############################################################################"
        filterflag="FAILED"
        detail=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
        echo -e "$detail" >> $badSamples
        sed -i "1i $detail"  $qc_result
    fi
    echo `date`
    cat $qc_result > $qc_result2
    cat $qc_result > $qc_result3

    if [ $filterflag == "FAILED" ]
    then
	echo -e "##########################################################################"
	echo -e "#############  $sample FAILED the verifyBamId filter     #################"
	echo -e "#############  analysis is STOPPED now                  #################"
	echo -e "##########################################################################"
	MSG="$sample DID NOT PASS the VerifyBamID filter. Stopping analysis now"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    else
	echo -e "##########################################################################"
	echo -e "#############  $sample PASSED the verifyBamId filter     #################"
	echo -e "#############  analysis continues with vcallgatk.sh      #################"
	echo -e "##########################################################################"
        exit 0;
    fi
