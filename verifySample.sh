#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#redmine=hpcbio-redmine@igb.illinois.edu
if [ $# -ne 11 ]
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
    realignOutput=$3
    filterflag=$4
    outfile=$5
    rootdir=$6
    elog=$7
    olog=$8
    email=$9
    qsubfile=${10}
    RealignOutputLogs=${11}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    echo -e "##########################################################################"
    echo -e "################  parse runfile and sanity check #########################"
    echo -e "##########################################################################"

    memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
    javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
    verifydir=$( cat $runfile | grep -w VERIFYBAMDIR | cut -d '=' -f2 )    
    freemix_cutoff=$( cat $runfile | grep -w FREEMIX_CUTOFF | cut -d '=' -f2 )
    refdir=$( cat $runfile | grep -w REFDIR | cut -d '=' -f2 )
    sites=$refdir/omni.sites.vcf.gz

    if [ ! -d $rootdir ]
    then
       MSG="$rootdir output root directory not found. verifySample failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $RealignOutput ]
    then
       MSG="$rootdir output realign directory not found. verifySample failed"
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
    mergedfile=${sample}.realigned.calmd.bam

    if [ ! -s $goodSamples ]
    then
        truncate -s 0 $goodsamples
    fi
    if [ ! -s $badSamples ]
    then
        truncate -s 0 $badsamples
    fi
    if [ $freemix_cutoff -eq $freemix_cutoff 2>/dev/null ]
    then
        echo -e "ok val, it is numeric"
    else
	MSG="FREEMIX_CUTFF=$freemix_cutff invalid value. verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi

    echo -e "##########################################################################"
    echo -e "#############  merging all chr files for sample=$sample  #################"
    echo -e "##########################################################################"
    `echo date`
    cd $RealignOutput

    chrbamfiles=`find ./ -name "*${sample}.realigned.calmd.bam"`
    bamsList=$( echo $chrbamfiles | sed "s/\.\// /g" | tr "\n" " " )
    if [ `expr ${#bamsList}` -lt 1 ]
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
        $javadir/java -Xmx1024m -Xms1024m -jar $picardir/MergeSamFiles.jar $bamsList OUTPUT=$mergedfile USE_THREADING=true 
        exitcode=$?
        `echo date`
        if [ $exitcode -ne 0 ]
        then
	    MSG="picard MergeSamFiles command failed exitcode=$exitcode verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
    fi

    echo -e "making sure  $mergedfile was created"

    if [ ! -s $mergedfile ]
    then
	MSG="picard MergeSamFiles command produced an empty file while trying to merge $bamsList verifySample failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
    fi

    echo -e "##########################################################################"
    echo -e "#############  running verifyBamID for sample=$sample    #################"
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
    if [ ! -s $outfile ]
    then
	MSG="verifyBamID produced empty file  exitcode=$exitcode verifySample failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    echo -e "##########################################################################"
    echo -e "#############  filtering for sample=$sample              #################"
    echo -e "##########################################################################"
    `echo date`

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
######################  editing stopped here   ######################################################################
    echo -e "produce a parsing statement here"
    freemix=$( cat $outfile | grep "the right line" | cut -f 7 ) 

    echo -e "checking that we actually grabbed something that is a number"

    if [ $freemix -eq $freemix 2>/dev/null ]
    then
        echo -e "ok val it is a number"
    else
	MSG="$outfile parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    if [ $freemix  -lt  $freemix_cutoff ]
    then
	echo -e "##########################################################################"
        echo -e "sample-$sample passed verifyBamID filter\nadd to list of good-samples"
	echo -e "##########################################################################"
        echo $sample >> $goodSamples
        filterflag="PASSED"
    else
	echo -e "############################################################################"
        echo -e "sample-$sample DID NOT passed verifyBamID filter\nadd to list of bad-samples"
	echo -e "############################################################################"
        echo $sample >> $badSamples
        filterflag="FAILED"
    fi
    echo `date`
