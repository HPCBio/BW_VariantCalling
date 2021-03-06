#!/bin/bash
#	
#  script to realign and recalibrate the aligned file(s) 
########################################################
redmine=hpcbio-redmine@igb.illinois.edu
set -x
if [ $# != 14 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
else					

	echo `date`
        umask 0027
	scriptfile=$0
        realrecaldir=$1
        outputfile=$2	
        chr=$3
        inputfile=$4
        RGparms=$5
        region=$6
        realparms=$7
        recalparms=$8
        runfile=$9
        flag=${10}
	elog=${11}
	olog=${12}
	email=${13}
        qsubfile=${14}
        
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


        RealignOutputLogs=`dirname $elog`

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
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
        sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
        novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
        dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
        realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
	outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
        thr=`expr $threads "-" 1`
        chrinfiles=$( echo $chrinfiles | tr ":" " " )
        RGparms=$( echo $RGparms | tr "::" ":" | sed "s/:/\tRG/g" | sed "s/ID\=/RGID\=/" )
        region=$( echo $region | tr ":" " " | sed "s/knownSites /--knownSites /g")
        real2parms=$( echo $realparms | tr ":" " " | sed "s/known /-known /g" )
        realparms=$( echo $realparms | tr ":" " " | sed "s/known /--known /g" )
        recalparms=$( echo $recalparms | tr ":" " " | sed "s/knownSites /--knownSites /g")
        skipRealign="NO"
        MillsIndels=$refdir/$indeldir/Mills_and_1000G_gold_standard.indels.b37.vcf
        OneKIndels=$refdir/$indeldir/1000G_phase1.indels.b37.vcf

        
        if [ ! -d $picardir ]
        then
	    MSG="$picardir picard directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge
"mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""	    exit 1;
        fi
        if [ ! -d $samdir ]
        then
	    MSG="$samdir samtools directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi
        if [ ! -d $sambambadir ]
        then
	    MSG="$sambambadir sambamba directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi

        if [ -z $javadir ]
        then
	    MSG="Value for JAVADIR must be specified in configuration file"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi
        if [ ! -s $refdir/$dbSNP ]
        then
	    MSG="$dbSNP DBSNP file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi
        if [ ! -s $MillsIndels ]
        then
	    MSG="$MillsIndels indels file  not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi
        if [ ! -s $OneKIndels ]
        then
	    MSG="$OneKIndelsIndels indels file  not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi

        if [ ! -s $inputfile ]
        then
	    MSG="$inputfile aligned bam file not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${lane}.Anisimov.msg
	    exit 1;
        fi


        if [ ! -d $realrecaldir ]
        then
	    MSG="$realrecaldir realign directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    exit 1;
        fi

        cd $realrecaldir

        set +x; echo -e "\n\n" >&2;
        echo "#################################################################################" >&2
        echo "################  STEP1: extract chromosome from aligned bam  ###################" >&2
        echo "#################################################################################" >&2
        echo -e "\n\n" >&2; set -x;
        
        cd $realrecaldir
	$samdir/samtools view -bu -@ $thr -h $inputfile $chr > presorted_wrg.${chr}.$infile  
	#$sambambadir/sambamba view -f bam -h -t $thr $inputfile $chr > presorted_wrg.${chr}.$infile  
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome samtools command failed exitcode=$exitcode  realignment-recalibration stopped for $infile"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
	    exit $exitcode;
	fi
	if [ ! -s presorted_wrg.${chr}.$infile ]
	then
	    MSG="split by chromosome samtools command failed it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for $infile"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
	    exit $exitcode;
	fi


        $novodir/novosort --index --threads $thr --tmpdir $realrecaldir -m 16g -o sorted_wrg.${chr}.$infile  presorted_wrg.${chr}.$infile 
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="split by chromosome sortsam command failed exitcode=$exitcode  realignment-recalibration stopped for $infile"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
	    exit $exitcode;
	fi
	if [ ! -s sorted_wrg.${chr}.$infile ]
	then
	    MSG="split by chromosome SortSam command failed. it produced an empty file exitcode=$exitcode  realignment-recalibration stopped for $infile"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
	    exit $exitcode;
	fi
	#$samdir/samtools index sorted_wrg.${chr}.$infile  
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
	    -o realign.$chr.$infile.list $realparms

	echo `date`

        echo -e "########### checking to see if the target list is empty                   #############"	
	if [ -z realign.$chr.$infile.list ]
	then
	    ## the target list is empty, we need to skip IndelAligner command
	    skipRealign="YES"
	    #MSG="realign.$chr.$infile.list realignertargetcreator file not created. realignment-recalibration stopped for $infile"
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
            #exit 1;
	fi

        if [ $skipRealign != "YES" ]
        then
             echo -e "########### the target list is NOT empty. Proceed with IndelRealigner #############"
	     $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    -I sorted_wrg.$chr.$infile \
	    -T IndelRealigner \
            -L $chr \
	    -o realign.$chr.$infile.realigned.bam \
	    -targetIntervals realign.$chr.$infile.list $realignparams $real2parms

 	    exitcode=$?
	    echo `date`
	    if [ $exitcode -ne 0 ]
	    then
	        MSG="indelrealigner command failed exitcode=$exitcode  realignment-recalibration stopped"
                echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
	        exit $exitcode;
	    fi
        else
             echo -e "########### the target list is empty. Skip IndelRealigner #############"
             cp sorted_wrg.$chr.$infile realign.$chr.$infile.realigned.bam
        fi

	if [ ! -s realign.$chr.$infile.realigned.bam ]
	then
	    MSG="realigned.bam indelrealigner file not created. realignment-recalibration stopped for $infile"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $RealignOutputLogs/FAILED_realrecal.${infile}.Anisimov.msg
            exit 1;
        fi


fi

