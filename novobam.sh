#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# -gt 14 ]
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else

	set -x
	echo `date`
        scriptfile=$0
        alignerdir=$1
        params=$2
        thr=$3
        ref=$4
	outputdir=$5
        samfile=$6
        bamfile=$7
        samdir=$8
        paired=$9
        BAM=${10}
        elog=${11}
        olog=${12}
        email=${13}
        qsubfile=${14}

        parameters=$( echo $params | tr "_" " " )
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        cd $outputdir

        if [ $paired == "1" ]
        then
	    $alignerdir/novoalign -F BAMPE -f $BAM -d $ref -o SAM $parameters > $samfile
        else
	    $alignerdir/novoalign -F BAMSE -f $BAM -d $ref -o SAM $parameters > $samfile
        fi
        exitcode=$?
	echo `date`

        if [ $exitcode -ne 0 ]
        then
           MSG="novoalign command failed.  exitcode=$exitcode alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit $exitcode;
        fi

        if [ ! -s $samfile ]
        then
           MSG="$samfile aligned BAM file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

	$samdir/samtools view -bS -o $bamfile $samfile
        exitcode=$?
	echo `date`

        if [ $exitcode -ne 0 ]
        then
           MSG="samtools view command failed.  exitcode=$exitcode  novobam alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit $exitcode;
        fi

	if [ ! -s $bamfile ]
	then
	    MSG="$bamfile BAM file not created. sam2bam step failed during  novobam alignment."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi       


        $samdir/samtools index $bamfile
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
           MSG="samtools index command failed.  exitcode=$exitcode  novobam alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit $exitcode;
        fi
        echo `date`

fi
