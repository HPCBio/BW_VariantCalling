#!/bin/sh
######################################
#  script to convert bam files back to fastq as pre requisite to alignment
#
######################################
# written in collaboration with Mayo Bioinformatics core group
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=lmainzer@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 9 ];
then
	MSG= "parameter mismatch"
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

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        sortdir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2 )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        bam2fastqparms=$( cat $runfile | grep -w BAM2FASTQPARMS | cut -d '=' -f2- )
        bam2fastqflag=$( cat $runfile | grep -w BAM2FASTQFLAG | cut -d '=' -f2 )
        javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
        provenance=$( cat $runfile | grep -w PROVENANCE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

        if [ ! -d $picardir ]
        then
	    MSG="PICARDIR=$picardir  directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      

        if [ -z $javamodule ]
        then
	    MSG="Value for JAVAMODULE must be specified in configuration file"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        else 
            #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
            `module load $javamodule`
        fi
        if [ $provenance != "SINGLE_SOURCE" -a $provenance != "MULTI_SOURCE" ]
        then
	    MSG="Invalid value for PROVENANCE=$provenance"
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
        flagstat=${outfilename}_preproc.flagstat
        temfile=${outfilename}_preproc
        cleanbam=${outfilename}.cleaned.bam

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

        $samdir/samtools view -bS $temfile".preprocessed.sam" > $temfile".preprocessed.bam"
        exitcode=$?
	echo `date`
        if [ $exitcode -ne 0 ]
        then
		MSG="sam2bam command failed. exitcode=$exitcode "
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
        fi
        ############################
        # revertsam (optional)
        ############################
        echo "optional conversion: revertsam"

        if [ $revertsam == "YES" -o $revertsam == "1" ]
        then
	    java -Xmx6g -Xms512m -jar $picardir/RevertSam.jar \
                 COMPRESSION_LEVEL=0 \
                 INPUT=$temfile".preprocessed.bam" \
                 OUTPUT=$temfile".preprocessed.reverted.bam" \
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

            if [ ! -s ${temfile}.preprocessed.reverted.bam ]
            then
		MSG="revertsam failed to create file ${temfile}.proprocessed.reverted.bam"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
            `mv ${temfile}.preprocessed.reverted.bam $cleanbam`
            exitcode=$?
	    echo `date`
        else
            `mv ${temfile}.preprocessed.bam $cleanbam`
            exitcode=$?
	    echo `date`
        fi
	
        if [ $exitcode -ne 0 ]
        then
	    MSG="no bam file produced for the final conversion step samtofastq. exitcode=$exitcode "
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        fi            


        if [ ! -s $cleanbam ]
        then
	    MSG="no bam file produced for the final conversion step samtofastq"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi

        #########################
        #  bam2fastq
        #########################
        echo "final conversion: samtofastq"

        if [ $provenance == "SINGLE_SOURCE" ]
        then
            echo "bam files come from a single source. split by sample"
	    java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
                          INPUT=$cleanbam \
                          TMP_DIR=$tmpdir \
                          OUTPUT_PER_RG=true \
                          OUTPUT_DIR=$tmpdir \
                          $bam2fastqparms \
                          VALIDATION_STRINGENCY=SILENT 
	    exitcode=$?
	    echo `date`
	    if [ $exitcode -ne 0 ]
	    then
	       MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
	       exit $exitcode;
	    fi
            if [ `find -name "*.fastq" | wc -l` -lt 1 ]
            then
		MSG="samtofastq command did not produced fastq files"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
	    echo `date`
            `rename _1.fastq _1.${outfilename}.fastq *_1.fastq`
            `rename _2.fastq _2.${outfilename}.fastq *_2.fastq`
        fi

        if [ $provenance == "MULTI_SOURCE" ]
        then
            echo "bam files from multi-sources. split by R1 and or R2"
	    R1=${outfilename}_R1.fastq
	    R2=${outfilename}_R2.fastq

	    if [ $paired == "1" ]
            then
                echo "paired-end reads will be generated"
		java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
		  		FASTQ=$R1 \
				SECOND_END_FASTQ=$R2 \
				INPUT=$cleanbam \
				TMP_DIR=$tmpdir \
				$bam2fastqparms \
				VALIDATION_STRINGENCY=SILENT 
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
		   MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		   exit $exitcode;
		fi
		echo `date`
		if [ ! -s $R1 -a ! -s $R2 ]
		then
		   MSG="$R1 $R2 FASTQ files not created. bam2fastq conversion failed."
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		   exit 1;
		 fi
            else
                echo "single reads will be generated"
		java -Xmx6g -Xms512m  -jar $picardir/SamToFastq.jar \
				FASTQ=$R1 \
				INPUT=$cleanbam \
				TMP_DIR=$tmpdir \
				$bam2fastqparms \
				VALIDATION_STRINGENCY=SILENT
		exitcode=$?
		echo `date`
		if [ $exitcode -ne 0 ]
		then
		   MSG="samtofastq command failed.  exitcode=$exitcode  bam2fastq conversion failed"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
		   exit $exitcode;
		fi
		echo `date`
		if [ ! -s $R1 ]
		then
		   MSG="$R1 Empty fasq file. bam2fastq conversion failed."
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		   exit 1;
		fi
	    fi
	fi

        `mv *.fastq $outputdir`
	echo `date`
fi
