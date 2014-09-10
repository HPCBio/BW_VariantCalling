#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# -gt 15 ]
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
    set -x
    echo `date`
    scriptfile=$0
    runfile=$1
    picardir=$2
    samdir=$3
    javamodule=$4
    outputdir=$5
    bamfile=$6
    infile=$7
    outfile=$8
    rgparms=$9
    flag=${10}
    chr=${11}
    elog=${12}
    olog=${13}
    email=${14}
    qsubfile=${15}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

    memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )


    if [ ! -d $outputdir ]
    then
       MSG="$outputdir realign directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       MSG="$picardir picard directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="$picardir samtools directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    if [ -z $javamodule ]
    then
	MSG="A value must be specified for JAVAMODULE in configuration file"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    else
        #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
            `module load $javamodule`
    fi          
    if [ ! -s $bamfile ]
    then
	MSG="$bamfile bam file to be sorted was not found"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    cd $outputdir
    echo `date`
    if [ ! -s $infile ]
    then
        $memprof $samdir/samtools view -b $bamfile $chr > $infile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="samtools view command failed exitcode=$exitcode. $infile bam file to be sorted within region:[$chr] was not created"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
        if [ ! -s $infile ]
        then
	    MSG="$infile bam file to be sorted within region:[$chr] was not created"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
	echo `date`
        $memprof $samdir/samtools index $infile
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="samtools index command failed exitcode=$exitcode. $infile bam file to be sorted within region:[$chr] was not created"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
	echo `date`
    fi

    tmpfile=tmp.wrg.$infile
    parameters=$( echo $rgparms | tr ":" " " )
    sortflag=$( echo $flag | tr '[a-z]' '[A-Z]' )

    ## before sorting, we need to make sure the bam file has readgroup info

    if [ $sortflag == "NCSA" ]
    then
       echo "alignment was done inhouse. we need to add_readgroup info"
       $memprof java -Xmx6g -Xms512m -jar $picardir/AddOrReplaceReadGroups.jar \
	   INPUT=$infile \
	   OUTPUT=$tmpfile \
	   MAX_RECORDS_IN_RAM=null \
	   TMP_DIR=$outputdir \
	   SORT_ORDER=unsorted \
           $parameters \
	   VALIDATION_STRINGENCY=SILENT

       exitcode=$?
       if [ $exitcode -ne 0 ]
       then
	    MSG="addorreplacereadgroup command failed exitcode=$exitcode  sortnode failed "
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
       fi
       echo `date`
    else
	echo "alignment was done at Mayo. checking if readgroup info is present"
	$memprof $samdir/samtools view -H $infile > $infile.header
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
	    MSG="samtools view command failed exitcode=$exitcode  sortnode failed "
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi
	echo `date`
	match=$( cat $file.header | grep '^@RG' )
	lenmatch=`expr ${#match}`
	if [ $lenmatch -gt 0 ]
	then
            echo "readgroup info found in input file."
            cp $infile $tmpfile
	else
            echo "readgroup info NOT found in input file. Adding it now..."
	    $memprof java -Xmx6g -Xms512m -jar $picardir/AddOrReplaceReadGroups.jar \
		   INPUT=$infile \
		   OUTPUT=$tmpfile \
		   MAX_RECORDS_IN_RAM=null \
		   TMP_DIR=$outputdir \
		   SORT_ORDER=unsorted \
		   $parameters \
		   VALIDATION_STRINGENCY=SILENT
            exitcode=$?
            if [ $exitcode -ne 0 ]
            then
		MSG="addorreplacereadgroup command failed exitcode=$exitcode  sortnode failed "
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi
	    echo `date`

	fi
    fi

    if [ ! -s $tmpfile ]
    then
	MSG="$tmpfile bam file not created. add_readGroup step failed. sortnode failed "
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi
    echo `date`

    $memprof java -Xmx6g -Xms512m -jar $picardir/SortSam.jar \
	INPUT=$tmpfile \
	OUTPUT=$outfile \
	TMP_DIR=$outputdir \
	SORT_ORDER=coordinate \
	MAX_RECORDS_IN_RAM=null \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="sortsam command failed exitcode=$exitcode  sortnode failed "
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`

    if [ ! -s $outfile ]
    then
	MSG="$outfile sort bam file not created.  sortnode failed "
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit 1;
    fi
    $memprof $samdir/samtools index $outfile
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="samtools index command failed exitcode=$exitcode  sortnode failed "
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`
    $memprof $samdir/samtools view -H $outfile > $outfile.header
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="samtools viewcommand failed exitcode=$exitcode  sortnode failed "
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`
fi
