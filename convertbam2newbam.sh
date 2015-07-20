#!/bin/sh
######################################
#  script to sort reads by name and remove singletons
#
######################################
# written in collaboration with Mayo Bioinformatics core group
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 9 ];
then
	MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        outputdir=$1
        tmpdir=$2
        infile=$3
        outfilename=$4
        runfile=$5
        elog=$6
        olog=$7
        email=$8
        qsubfile=$9
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        sortdir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2 )

        if [ ! -s $infile ]
        then
	    MSG="$infile input bam file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $sortdir ]
        then
	    MSG="$sortdir  sorttool directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      

        if [ $revertsam != "YES" -a $revertsam != "1" -a $revertsam != "NO" -a $revertsam != "0" ]
        then
	    MSG="Invalid value for REVERTSAM=$revertsam"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      

        outfile=$outputdir/$outfilename
        flagstat=${outfilename}_preprocessing.flagstat
        temfile=${outfilename}_preprocessing
        cd $tmpdir
	echo `date`

        ############################
        # sort by name
        ############################       

        echo "sorting bam by readname..."
        $sortdir/novosort --namesort --threads $thr $infile > $temfile".name_sorted"

        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="sortbyname command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi

        if [ ! -s $temfile".name_sorted" ]
        then
		MSG="${temfile}.name_sorted  -preprocessing of bam file failed to create bam file sorted by name"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	fi
	echo `date`

        $samdir/samtools view $temfile".name_sorted" > $temfile".name_sorted.sam"
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="bam2sam command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi

        $samdir/samtools view -H $temfile".name_sorted" > $temfile".name_sorted.header"
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
                MSG="header generation failed. exitcode=$exitcode "
                echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit $exitcode;
        fi

        ############################
        # remove singletons
        ############################       
        echo "removing singletons from bam"
        awk '{print $1}' $temfile".name_sorted.sam" > $temfile".name_sorted.sam.names_only"
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="extracting name tags from sam file, command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi

        uniq -c $temfile".name_sorted.sam.names_only" > $temfile".name_sorted.sam.names_only.counted"
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="counting tags, command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi

        awk '{if ($1==2) print $2}' $temfile".name_sorted.sam.names_only.counted" > $temfile".name_sorted.sam.names_only.counted.occur_only_twice"
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="finding unique tags, command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi

        $scriptdir/FindDifferencesInOutput.pl  $temfile".name_sorted.sam.names_only.counted.occur_only_twice" $temfile".name_sorted.sam"  $temfile".reads_not_twice" $temfile".preprocessed_no_header.sam"
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="FindDifferencesInOutput program failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi

        cat $temfile".name_sorted.header" $temfile".preprocessed_no_header.sam" > $temfile".preprocessed.sam"


        ############################
        # revertsam (optional)
        ############################
        echo "final conversion sam2bam; may include revertsam"

        if [ $revertsam == "YES" -o $revertsam == "1" ]
        then
	    java -Xmx6g -Xms512m -jar $picardir/RevertSam.jar \
                 COMPRESSION_LEVEL=0 \
                 INPUT=$temfile".preprocessed.sam" \
                 OUTPUT=$temfile".preprocessed.reverted.sam" \
                 VALIDATION_STRINGENCY=SILENT

            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="revertsam command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi                           
            $samdir/samtools view -bS $temfile".preprocessed.reverted.sam" > $outfile
            exitcode=$?
	    echo `date`
        else          
            $samdir/samtools view -bS $temfile".preprocessed.sam" > $outfile
            exitcode=$?
	    echo `date`
        fi
        if [ $exitcode -ne 0 ]
        then
		MSG="sam2bam command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi            

        if [ ! -s $outfile ]
        then
	    MSG="$outfile  -preprocessing of bam file failed to create file"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
	echo `date`

        $samdir/samtools index $outfile
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
	    MSG="indexing command failed for file $outfile exitcode=$exitcode "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

	java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
        INPUT=$outfile \
        OUTPUT=$flagstat \
        VALIDATION_STRINGENCY=SILENT

        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
	    MSG="CollectAlignmentSummaryMetrics command failed on outfile $outfile. exitcode=$exitcode "
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi

fi
