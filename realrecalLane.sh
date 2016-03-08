#!/bin/bash
#	
#  script to realign and recalibrate by lane in multiplexed samples 
########################################################
redmine=hpcbio-redmine@igb.illinois.edu
set -x
if [ $# != 14 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

	echo `date`
        umask 0027
	scriptfile=$0
        lane=$1
        realrecaldir=$2
        outputfile=$3
        chr=$4
        inputfile=$5
        RGparms=$6
        region=$7
        realparms=$8
        recalparms=$9
        runfile=${10}
	elog=${11}
	olog=${12}
	email=${13}
        qsubfile=${14}
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        RealignOutputLogs=`dirname $elog`

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi
        

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;
        
        javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
        novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
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
        region=$( echo $region | tr ":" " " | sed "s/knownSites /--knownSites /g")

        if [ ! -d $picardir ]
        then
	    MSG="$picardir picard directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
        fi
        if [ ! -d $samdir ]
        then
	    MSG="$samdir samtools directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi

        if [ -z $javadir ]
        then
	    MSG="Value for JAVADIR must be specified in configuration file"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        #else
            #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
        #    `module load $javamodule`
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi
        if [ ! -s $inputfile ]
        then
	    MSG="$inputfile aligned bam file not found for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi


        if [ ! -d $realrecaldir ]
        then
	    MSG="$rearecaldir realign directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit 1;
        fi

        infile=`basename $inputfile`
        
        set +x; echo -e "\n\n" >&2;
        echo "#################################################################################" >&2
        echo "################  STEP1: extract chromosome from aligned bam  ###################" >&2
        echo "#################################################################################" >&2
        echo -e "\n\n" >&2; set -x;
        
        cd $realrecaldir
	$samdir/samtools view -bu -@ $thr -h $inputfile $chr > presorted_norg.${chr}.$infile  
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome samtools command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi
	if [ ! -s presorted_norg.${chr}.$infile ]
	then
	    MSG="split by chromosome samtools command failed it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi

	ln -s presorted_norg.${chr}.$infile  presorted_wrg.${chr}.$infile

        $novodir/novosort --threads $thr --tmpdir $realrecaldir -m 16g presorted_wrg.${chr}.$infile > sorted_wrg.${chr}.$infile 
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome sortsam command failed exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi
	if [ ! -s sorted_wrg.${chr}.$infile ]
	then
	    MSG="split by chromosome SortSam command failed. it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi
	$samdir/samtools index sorted_wrg.${chr}.$infile  
	$samdir/samtools view -H  sorted_wrg.${chr}.$infile > sorted_wrg.${chr}.${infile}.header  
	echo `date`


        set +x; echo -e "\n\n" >&2;
        echo "#################################################################################" >&2
        echo "################  STEP2: realign                              ###################" >&2
        echo "#################################################################################" >&2
        echo -e "\n\n" >&2; set -x;
				
		
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
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi

	
	if [ ! -s realign.$chr.$lane.list ]
	then
	    MSG="realign.$chr.$lane.list realignertargetcreator file not created. realignment-recalibration stopped for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
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
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi


	if [ ! -s realign.$chr.$lane.realigned.bam ]
	then
	    MSG="realigned.bam indelrealigner file not created. realignment-recalibration stopped for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
            exit 1;
        fi

        set +x; echo -e "\n\n" >&2;
        echo "#################################################################################" >&2
        echo "################  STEP3: recalibration                        ###################" >&2
        echo "#################################################################################" >&2
        echo -e "\n\n" >&2; set -x;

        if [ "$recalibrator" == "BQSR" ]
        then
            set +x; echo -e "\n\n" >&2;       
            echo "#################################################################################" >&2
	    echo "                            recalibrator == BQSR" >&2
            echo "#################################################################################" >&2
            echo -e "\n\n" >&2; set -x;         

            $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
                -R $refdir/$ref \
                $recalparms \
		$region \
                -I realign.$chr.$lane.realigned.bam \
                -T BaseRecalibrator \
                --out recal.$chr.$lane.recal_report.grp \
                -nct $thr

            exitcode=$?
            echo `date`
            if [ $exitcode -ne 0 ]
            then
                MSG="BQSR command failed exitcode=$exitcode  recalibration stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
                exit $exitcode;
            fi


            if [ ! -s recal.$chr.$lane.recal_report.grp ]
            then
                MSG="recal.$chr.$lane.recal_report.grp recalibration report file not created, realignment-recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
                exit 1;
            fi
            echo `date`


            $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
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
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
                exit $exitcode;
            fi


            if [ ! -s recal.$chr.$lane.real.recal.bam ]
            then
                MSG="recal.$chr.real.recal.bam recalibrated file not created, realignment-recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
                exit 1;
            fi

            $samdir/samtools index recal.$chr.$lane.real.recal.bam
            ln -s $realrecaldir/recal.$chr.$lane.real.recal.bam $outputfile
            ls -s $realrecaldir/recal.$chr.$lane.real.recal.bam.bai $outputfile.bai

        else
            set +x; echo -e "\n\n" >&2;       
            echo "#################################################################################" >&2
	    echo "                            recalibrator IS NOT BQSR" >&2
            echo "#################################################################################" >&2
            echo -e "\n\n" >&2; set -x;         
            
   	    $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		-R $refdir/$ref \
		$recalparms \
		$region \
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
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
		exit $exitcode;
	    fi

	    if [ ! -s recal.$chr.$lane.recal_data.csv ]
	    then
		MSG="recal.$chr.$lane.recal_data.csv countcovariates file not created realignment-recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
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
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
		exit $exitcode;
	    fi
	    echo `date`		


	    if [ ! -s recal.$chr.$lane.real.recal.bam ]
	    then
		MSG="$chr.real.recal.bam tablerecalibration file not created realignment-recalibration stopped for lane=$lane"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
		exit 1;
	    fi


            $samdir/samtools index recal.$chr.$lane.real.recal.bam
            ln -s $realrecaldir/recal.$chr.$lane.real.recal.bam $outputfile
            ln -s $realrecaldir/recal.$chr.$lane.real.recal.bam.bai $outputfile.bai
        fi # end choosing between BQSR and CountCovariates/TableRecalibration

	$samdir/samtools flagstat $outputfile > $outputfile.flagstat
	exitcode=$?
	echo `date`		
	if [ $exitcode -ne 0 ]
	then
	    MSG="samtools flagstat command failed with outputfile=$outputfile exitcode=$exitcode for lane=$lane"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecalLane.${lane}.Anisimov.msg
	    exit $exitcode;
	fi
