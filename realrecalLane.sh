#!/bin/bash
# written in collaboration with Mayo Bioinformatics core group
#	
#  script to realign and recalibrate by lane in multiplexed samples 
########################################################
#redmine=hpcbio-redmine@igb.illinois.edu
set -x
if [ $# != 13 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

	echo `date`
	scriptfile=$0
        lane=$1
        realrecaldir=$2
        outputfile=$3
        chr=$4
        inputfile=$5
        RGparms=$6
        realparms=$7
        recalparms=$8
        runfile=$9
	elog=${10}
	olog=${11}
	email=${12}
        qsubfile=${13}
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
        #dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        #kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
        #targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
        realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
	outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
        thr=`expr $threads "-" 1`
        #thr=`expr $threads "/" 2`

        # cleaning up the lists
        #chrinfiles=$( echo $chrinfiles | tr ":" " " )
#        chrinputfiles=$( echo $chrinputfiles | tr ":" " " ) # not used anymore?
        RGparms=$( echo $RGparms | tr "::" ":" | sed "s/:/\tRG/g" | sed "s/ID\=/RGID\=/" )
        #rgheader=$( echo $RGparms | tr "=" ":" | tr " " "\t" )
        real2parms=$( echo $realparms | tr ":" " " | sed "s/known /-known /g" )
        realparms=$( echo $realparms | tr ":" " " | sed "s/known /--known /g" )
        recalparms=$( echo $recalparms | tr ":" " " | sed "s/knownSites /--knownSites /g")

        if [ ! -d $picardir ]
        then
	    MSG="$picardir picard directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge
"mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""	    exit 1;
        fi
        if [ ! -d $samdir ]
        then
	    MSG="$samdir samtools directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ -z $javadir ]
        then
	    MSG="Value for JAVADIR must be specified in configuration file"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        #else
            #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
        #    `module load $javamodule`
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $inputfile ]
        then
	    MSG="$inputfile aligned bam file not found for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        if [ ! -d $realrecaldir ]
        then
	    MSG="$rearecaldir realign directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        infile=`basename $inputfile`
        echo "#################################################################################"
	echo "################   STEP1: extract chromosome from aligned bam ###################"
        echo "#################################################################################"

        cd $realrecaldir
	$samdir/samtools view -b -h $inputfile $chr > presorted_norg.${chr}.$infile  
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome samtools command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	if [ ! -s presorted_norg.${chr}.$infile ]
	then
	    MSG="split by chromosome samtools command failed it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi

	
        $javadir/java -Xmx1024m -Xms1024m -jar $picardir/AddOrReplaceReadGroups.jar \
	    INPUT=presorted_norg.${chr}.$infile \
	    OUTPUT=presorted_wrg.${chr}.$infile \
	    TMP_DIR=$realrecaldir \
	    SORT_ORDER=coordinate \
	    $RGparms

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="picard addreadgroup command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
        
	$javadir/java -Xmx1024m -Xms1024m -jar $picardir/SortSam.jar \
	    INPUT=presorted_wrg.${chr}.$infile \
	    OUTPUT=sorted_wrg.${chr}.$infile \
	    TMP_DIR=$realrecaldir \
	    SORT_ORDER=coordinate 

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome sortsam command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	if [ ! -s sorted_wrg.${chr}.$infile ]
	then
	    MSG="split by chromosome SortSam command failed. it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome addreadgroups command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
	if [ ! -s sorted_wrg.${chr}.$infile ]
	then
	    MSG="split by chromosome addreadgroups command failed. it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi

        echo "split by chromosome was generated=$chr.$infile"

	$samdir/samtools index sorted_wrg.${chr}.$infile  
	$samdir/samtools view -H  sorted_wrg.${chr}.$infile > sorted_wrg.${chr}.${infile}.header  
	echo `date`

        echo "#################################################################################"
	echo "##############################   STEP2: realign              ####################"
        echo "#################################################################################"
		
        $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    -I sorted_wrg.$chr.$infile \
	    -T RealignerTargetCreator \
            -nt $thr \
	    -o realign.$chr.$lane.list $realparms

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="realignertargetcreator command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi

	
	if [ ! -s realign.$chr.$lane.list ]
	then
	    MSG="realign.$chr.$lane.list realignertargetcreator file not created. realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
	fi
	echo `date`

	$javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    -I sorted_wrg.$chr.$infile \
	    -T IndelRealigner \
            -L $chr \
	    -o realign.$chr.$lane.realigned.bam \
	    -targetIntervals realign.$chr.$lane.list $realignparams $real2parms

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="indelrealigner command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
	fi


	if [ ! -s realign.$chr.$lane.realigned.bam ]
	then
	    MSG="realigned.bam indelrealigner file not created. realignment-recalibration stopped for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi

        echo "#################################################################################"
	echo "##############################   STEP3: recalibration       #####################"
        echo "#################################################################################"

        if [ "$recalibrator" == "BQSR" ]
        then
	    echo "recalibrator == BQSR"

            $memprov $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
                -R $refdir/$ref \
                $recalparms \
                -I realign.$chr.$lane.realigned.bam \
                -T BaseRecalibrator \
                --out recal.$chr.$lane.recal_report.grp \
                -nct $thr

            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="BQSR command failed exitcode=$exitcode  recalibration stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit $exitcode;
            fi


            if [ ! -s recal.$chr.$lane.recal_report.grp ]
            then
                MSG="recal.$chr.$lane.recal_report.grp recalibration report file not created, realignment-recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
            echo `date`


            $memprov $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
                -R $refdir/$ref \
                -I realign.$chr.$lane.realigned.bam \
                -T PrintReads \
                -BQSR recal.$chr.$lane.recal_report.grp \
                --out recal.$chr.$lane.real.recal.bam \
                -nct $thr

            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="BQSR PrintReads command failed exitcode=$exitcode  recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit $exitcode;
            fi


            if [ ! -s recal.$chr.$lane.real.recal.bam ]
            then
                MSG="recal.$chr.real.recal.bam recalibrated file not created, realignment-recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
            $samdir/samtools index recal.$chr.$lane.real.recal.bam
            cp recal.$chr.$lane.real.recal.bam $outputfile
            cp recal.$chr.$lane.real.recal.bam.bai $outputfile.bai

        else
	    echo "recalibrator =! BQSR"
   	    $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		-R $refdir/$ref \
		$recalparms \
		-I realign.$chr.$lane.realigned.bam \
		-T CountCovariates \
		-cov ReadGroupCovariate \
		-cov QualityScoreCovariate \
		-cov CycleCovariate \
		-cov DinucCovariate \
		-recalFile recal.$chr.$lane.recal_data.csv 

	    exitcode=$?
	    echo `date`		
	    if [ $exitcode -ne 0 ]
	    then
		MSG="countcovariates command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi

	    if [ ! -s recal.$chr.$lane.recal_data.csv ]
	    then
		MSG="recal.$chr.$lane.recal_data.csv countcovariates file not created realignment-recalibration stopped for lane=$lane"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
	    echo `date`

	    $javadir/java -Xmx8g -Xms1024m  -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		-R $refdir/$ref \
		-L $chr \
		-I realign.$chr.$lane.realigned.bam \
		-T TableRecalibration \
		--out recal.$chr.$lane.real.recal.bam \
		-recalFile recal.$chr.$lane.recal_data.csv 

	    exitcode=$?
	    if [ $exitcode -ne 0 ]
	    then
		MSG="tablerecalibration command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit $exitcode;
	    fi
	    echo `date`		


	    if [ ! -s recal.$chr.$lane.real.recal.bam ]
	    then
		MSG="$chr.real.recal.bam tablerecalibration file not created realignment-recalibration stopped for lane=$lane"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi


            $samdir/samtools index recal.$chr.$lane.real.recal.bam
            cp recal.$chr.$lane.real.recal.bam $outputfile
            cp recal.$chr.$lane.real.recal.bam.bai $outputfile.bai
        fi # end choosing between BQSR and CountCovariates/TableRecalibration

	$samdir/samtools flagstat $outputfile > $outputfile.flagstat
	exitcode=$?
	echo `date`		
	if [ $exitcode -ne 0 ]
	then
	    MSG="samtools flagstat command failed with outputfile=$outputfile exitcode=$exitcode for lane=$lane"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
	fi
