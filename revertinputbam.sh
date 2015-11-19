#!/bin/bash
######################################
#  script to revertsam on input bamfiles that need to be realigned-recal
#
######################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 9 ];
then
	MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        umask 0027

        scriptfile=$0
        infile=$1
        outfile=$2
        picardir=$3
        samdir=$4
        javadir=$5
        elog=$6
        olog=$7
        email=$8
        qsubfile=$9
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $infile ]
        then
	    MSG="$infile input bam file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $picardir ]
        then
	    MSG="PICARDIR=$picardir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      

        if [ ! -d $samdir ]
        then
	    MSG="SAMDIR=$samdir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      
	outputdir=`dirname $outfile`
        cd $outputdir
        prefix=`basename $infile`
        samfile=$outputdir/${prefix}.sam
	echo `date`
        $javadir/java -Xmx1024m -Xms1024m -jar $picardir/RevertSam.jar \
            COMPRESSION_LEVEL=0 \
            INPUT=$infile \
            OUTPUT=$samfile \
            VALIDATION_STRINGENCY=SILENT

        exitcode=$?
	echo `date`

        if [ $exitcode -ne 0 ]
        then
	    MSG="revertsam command failed. exitcode=$exitcode "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
        if [ ! -s $samfile ]
        then
	    MSG="$samfile file not created. revertsam failed"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

	$samdir/samtools view -bS -o $outfile $samfile
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
	    MSG="samtools view command failed. exitcode=$exitcode "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

        if [ ! -s $outfile ]
        then
	    MSG="$outfile file not created. revertsam failed"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
fi
