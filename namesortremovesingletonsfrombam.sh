#!/bin/sh
######################################
#  script to sort reads by name and remove singletons
#
######################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 11 ];
then
	MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        infile=$1
        outfile=$2
        picardir=$3
        samdir=$4
        sortdir=$5
        scriptdir=$6
        thr=$7
        elog=$9
        olog=${10}
        email=${11}
        qsubfile=${12}
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $infile ]
        then
	    MSG="$infile input bam file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $sortdir ]
        then
	    MSG="$sortdir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      

	outputdir=`dirname $outfile`
        sufix=`basename $outfile`
        prefix=$( echo $RANDOM )
        flagstat=${prefix}_${sufix}_preprocessing.flagstat
        temfile=${prefix}_${sufix}
        cd $outputdir
	echo `date`
        echo "sorting bam by readname..."
            $sortdir/novosort --namesort --threads $thr $infile > $temfile".name_sorted"

            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="sortbyname command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

            if [ ! -s $temfile".name_sorted" ]
            then
		MSG="${temfile}.name_sorted  -preprocessing of bam file failed to create bam file sorted by name"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
	    echo `date`

            $samdir/samtools view $temfile".name_sorted" > $temfile".name_sorted.sam"
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="bam2sam command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

            $samdir/samtools view -H $temfile".name_sorted" > $temfile".name_sorted.header"
            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="header generation failed. exitcode=$exitcode "
                echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                exit $exitcode;
            fi


            awk '{print $1}' $temfile".name_sorted.sam" > $temfile".name_sorted.sam.names_only"
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="extracting name tags from sam file, command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

            uniq -c $temfile".name_sorted.sam.names_only" > $temfile".name_sorted.sam.names_only.counted"
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="counting tags, command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
            awk '{if ($1==2) print $2}' $temfile".name_sorted.sam.names_only.counted" > $temfile".name_sorted.sam.names_only.counted.occur_only_twice"
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="finding unique tags, command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

            $scriptdir/FindDifferencesInOutput.pl  $temfile".name_sorted.sam.names_only.counted.occur_only_twice" $temfile".name_sorted.sam"  $temfile".reads_not_twice" $temfile".preprocessed_no_header.sam"
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="FindDifferencesInOutput program failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
            cat $temfile".name_sorted.header" $temfile".preprocessed_no_header.sam" > $temfile".preprocessed.sam"
            $samdir/samtools view -bS $temfile".preprocessed.sam" > $outfile
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="sam2bam command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi            

        if [ ! -s $outfile ]
        then
	    MSG="$outfile  -preprocessing of bam file failed to create file"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
	echo `date`

	java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
        INPUT=$outfile \
        OUTPUT=$flagstat \
        VALIDATION_STRINGENCY=SILENT

        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
	    MSG="CollectAlignmentSummaryMetrics command failed on outfile $outfile. exitcode=$exitcode "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

fi
