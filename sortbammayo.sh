#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 10 ] 
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else
    set -x
    echo `date`
    umask 0027
    outputdir=$1
    inbamfile=$2
    sortedplain=$3
    outfilewdups=$4
    outfilenodups=$5
    runfile=$6
    elog=$7
    olog=$8
    email=$9
    scriptfile=$0
    qsubfile=${10}
    LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"



    if [ ! -s $runfile ]
    then
       MSG="$runfile configuration file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    set +x; echo -e "\n\n#############      parameter checking ###############\n\n" >&2; set -x;

    sample=`basename $outputdir .bam`
    picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
    samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
    markdup=$( cat $runfile | grep -w ^MARKDUP | cut -d '=' -f2 )
    deldup=$( cat $runfile | grep -w ^REMOVE_DUP | cut -d '=' -f2 )
    revertsam=$( cat $runfile | grep -w ^REVERTSAM | cut -d '=' -f2 )
    javadir=$( cat $runfile | grep -w ^JAVADIR | cut -d '=' -f2 )
    sID=$sample
    sPU=$sample
    sSM=$sample
    sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )       
    sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )        
    sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
    RGparms=$( echo "RGID=${sID} RGLB=${sLB} RGPL=${sPL} RGPU=${sPU} RGSM=${sSM} RGCN=${sCN}" )

    if [ ! -d $outputdir ]
    then
       MSG="$outputdir output directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $picardir ]
    then
       MSG="$picardir picard directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ ! -d $samdir ]
    then
       MSG="$samdir samtools directory not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi
    if [ -z $javadir ]
    then
       MSG="A value for JAVAdir has to be specified in configuration file."
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    #else
        #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
    #        `module load $javamodule`
    fi

    if [ $revertsam != "1" -a $revertsam != "0" -a $revertsam != "YES" -a $revertsam != "NO" ]
    then
           MSG="Invalid value for REVERTSAM=$revertsam"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$email""
            exit 1;
    else
        if [ $revertsam == "YES" ]
        then
                $revertsam="1"
        fi
        if [ $revertsam == "NO" ]
        then
                $revertsam="0"
        fi
    fi


    cd $outputdir

    if [ ! -s $inbamfile ]
    then
       MSG="$outputdir/$inbamfile BAM file not found"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
       exit 1;
    fi

    set +x; echo -e "\n\n#############    params are ok                                 ###############\n\n" >&2; set -x;


    bamfile=`basename $inbamfile`
    tmpfile=temp.$bamfile
    tmpfilerevertsam=tmprsam.$bamfile

    if [ $revertsam == "1" ]
    then
    
         set +x; echo -e "\n\n#############    step 1: reversam (optional)                   ###############\n\n" >&2; set -x;
  
        echo "revertsam then addreadgroup..."
	$javadir/java -Xmx1024m -Xms1024m -jar $picardir/RevertSam.jar \
        COMPRESSION_LEVEL=0 \
        INPUT=$inbamfile \
        OUTPUT=$tmpfilerevertsam \
	VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="revertsam command failed exitcode=$exitcode  sortbammayo failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`		

	if [ ! -s $tmpfilerevertsam ]
	then
	    MSG="$tmpfilerevertsam bam file not created.  sortbammayo failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi    

        set +x; echo -e "\n\n#############    step 2: addreadgroup               ###############\n\n" >&2; set -x;


	$javadir/java -Xmx1024m -Xms1024m -jar $picardir/AddOrReplaceReadGroups.jar \
	    INPUT=$tmpfilerevertsam \
	    OUTPUT=$tmpfile \
	    MAX_RECORDS_IN_RAM=null \
	    TMP_DIR=$outputdir \
	    SORT_ORDER=unsorted \
	    $RGparms \
	    VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="addorreplacereadgroup command failed exitcode=$exitcode  sortbammayo failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`		

    else
        set +x; echo -e "\n\n#############    step 1: reversam (skipped)        ###############\n\n" >&2; set -x;
    
        set +x; echo -e "\n\n#############    step 2: addreadgroup              ###############\n\n" >&2; set -x;


	$javadir/java -Xmx1024m -Xms1024m -jar $picardir/AddOrReplaceReadGroups.jar \
	    INPUT=$inbamfile \
	    OUTPUT=$tmpfile \
	    MAX_RECORDS_IN_RAM=null \
	    TMP_DIR=$outputdir \
	    SORT_ORDER=unsorted \
	    $RGparms \
	    VALIDATION_STRINGENCY=SILENT
	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="addorreplacereadgroup command failed exitcode=$exitcode  sortbammayo failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`		
    fi

    if [ ! -s $tmpfile ]
    then
	MSG="$tmpfile bam file not created. add_readGroup step failed sortbammayo failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit 1;
    fi    
    echo `date`

    set +x; echo -e "\n\n#############    step 3: sortsam                 ###############\n\n" >&2; set -x;

    $javadir/java -Xmx1024m -Xms1024m -jar $picardir/SortSam.jar \
	INPUT=$tmpfile \
	OUTPUT=$sortedplain \
	TMP_DIR=$outputdir \
	SORT_ORDER=coordinate \
	MAX_RECORDS_IN_RAM=null \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="sortsam command failed exitcode=$exitcode  sortbammayo failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit $exitcode;
    fi
    echo `date`		

    if [ ! -s $sortedplain ]
    then
	MSG="$sortedplain sorted bam file not created.  sortbammayo failed"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	exit 1;
    fi

    $samdir/samtools index $sortedplain
    $samdir/samtools flagstat $sortedplain > $sortedplain.flagstat
    $samdir/samtools view -H $sortedplain > $sortedplain.header
    rm $tmpfile

    echo `date`
        
    set +x; echo -e "\n\n#############    step 4: mardups  step               ###############\n\n" >&2; set -x;


    if [ $markdup == "YES" -a $deldup != "TRUE" ]
    then

        set +x; echo -e "\n\n#############    step 4: mark but keep   duplicates       ###############\n\n" >&2; set -x;

        echo "marking duplicates in sorted bam file"
        $javadir/java -Xmx1024m -Xms1024m -jar $picardir/MarkDuplicates.jar \
	    INPUT=$sortedplain \
	    OUTPUT=$outfilewdups \
	    TMP_DIR=$outputdir \
	    METRICS_FILE=$outfilewdups.dup.metrics \
	    MAX_RECORDS_IN_RAM=null \
	    CREATE_INDEX=true \
	    VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	if [ $exitcode -ne 0 ]
	then
	    MSG="markduplicates command failed exitcode=$exitcode  sortbammayo failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	echo `date`			    
	
	if [ ! -s $outfilewdups ]
	then
	    MSG="$outfilewdups file not created. markDuplicates step failed sortbammayo failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""	    
	    exit 1;
	fi
        echo "indexing bam file w marked duplicates"
	$samdir/samtools index $outfilewdups
	$samdir/samtools flagstat $outfilewdups > $outfilewdups.flagstat
	$samdir/samtools view -H $outfilewdups > $outfilewdups.header
    else
        set +x; echo -e "\n\n#############    step 4: mark and remove   duplicates       ###############\n\n" >&2; set -x;

	if [ $deldup == "TRUE" ]
	then
            echo "removing marked duplicates in sorted bam file"
            $javadir/java -Xmx1024m -Xms1024m -jar $picardir/MarkDuplicates.jar \
		INPUT=$sortedplain \
		OUTPUT=$outfilenodups \
		TMP_DIR=$outputdir \
		METRICS_FILE=$outfilenodups.dup.metrics \
		MAX_RECORDS_IN_RAM=null \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=true \
		VALIDATION_STRINGENCY=SILENT

	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="markduplicates command failed exitcode=$exitcode  sortbammayo failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi
	    echo `date`			    

	    if [ ! -s $outfilenodups ]
	    then
		MSG="$outfilenodups file not created. RemoveDuplicates step failed sortbammayo failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
            echo "indexing bam file w removed duplicates"
	    $samdir/samtools index $outfilenodups
	    $samdir/samtools flagstat $outfilenodups > $outfilenodups.flagstat
	    $samdir/samtools view -H $outfilenodups > $outfilenodups.header

	    echo `date`
	else
            echo "remove duplicates; job not requested"
            `cp $sortedplain $outfilewdups`
            `cp $sortedplain.flagstat $outfilewdups.flagstat`
            `cp $sortedplain.header $outfilewdups.header`
	    echo `date`
        fi
    fi
fi
