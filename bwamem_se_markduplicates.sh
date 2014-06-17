#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 15 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        aligndir=$1
        parms=$2
        ref=$3
        outputdir=$4
	bamfile=$5
        Rone=$6
        scriptdir=$7
        samdir=$8
        samblasterdir=$9
        picardir=${10}
        elog=${11}
        olog=${12}
        email=${13}
        qsubfile=${14}
        RGparms=${15}

        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        parameters=$( echo $parms | tr "_" " " )



        header=$( echo $RGparms  | tr ":" "\t" )
        rgheader=$( echo -n -e "@RG\t" )$( echo $header  | tr "=" ":" )



########## step 1: align, mark duplicates, convert to bam - on the fly
        cd $outputdir
        $aligndir/bwa mem $parameters -R "${rgheader}" $ref $Rone | $samblasterdir/samblaster |  $samdir/samtools view -bS -> ${bamfile}.wdups

        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem samblaster samtools command failed on $Rone.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s ${bamfile}.wdups ]
        then
            MSG="${bamfile}.wdups aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        echo `date`


########## step 2: indexing merged bam file    
        $samdir/samtools index ${bamfile}.wdups
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
             MSG="samtools index command failed.  exitcode=$exitcode. bwamem_pe_markduplicates stopped "
             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
             exit $exitcode;
        fi
        echo `date`

        #$samdir/samtools flagstat $tmpfilewdups > $tmpfilewdups.flagstat
        $samdir/samtools view -H ${bamfile}.wdups > ${bamfile}.wdups.header
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
             MSG="samtools view command failed.  exitcode=$exitcode. bwamem_pe_markduplicates stopped "
             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
             exit $exitcode;
        fi

        echo `date`
        java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
            INPUT=${bamfile}.wdups \
            OUTPUT=${bamfile}.wdups.flagstat \
            VALIDATION_STRINGENCY=SILENT

         exitcode=$?
         if [ $exitcode -ne 0 ]
         then
             MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode . bwamem_pe_markduplicates stopped "
             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
             exit $exitcode;
         fi
        echo `date`



fi
