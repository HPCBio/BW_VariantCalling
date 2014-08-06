#!/bin/sh
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 14 ]
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
	bamprefix=$5
        Rone=$6
        Rtwo=$7
        runfile=$8
        elog=${9}
        olog=${10}
        email=${11}
        qsubfile=${12}
        RGparms=${13}
        AlignOutputLogs=${14}

        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"



        if [ ! -s $runfile ]
        then
           MSG="$runfile configuration file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )

        samprocessing=$( cat $runfile | grep -w SAMPROCESSING | cut -d "=" -f2 )

        samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
        samblasterdir=$( cat $runfile | grep -w SAMBLASTERDIR | cut -d "=" -f2 )
        sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d "=" -f2 )
        novodir=$( cat $runfile | grep -w NOVODIR | cut -d "=" -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d "=" -f2 )
        alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )


        header=$( echo $RGparms  | tr ":" "\t" )
        rgheader=$( echo -n -e "@RG\t" )$( echo $header  | tr "=" ":" )



########## step 1: align, mark duplicates, convert to bam - on the fly
        cd $outputdir
        if [ $samprocessing == "SAMTOOLS" ]
        then
           echo `date`
           $aligndir/bwa mem $alignparms -R "${rgheader}" $ref $Rone $Rtwo | $samblasterdir/samblaster |  $samdir/samtools view -bS -> ${bamprefix}.wdups
           exitcode=$?
           echo `date`
        elif [ $samprocessing == "SAMBAMBA" ]
        then 
           echo `date`
           $aligndir/bwa mem $alignparms -R "${rgheader}" $ref $Rone $Rtwo | $samblasterdir/samblaster -o ${bamprefix}.sam.wdups
           exitcode=$?
           echo `date`
           $sambambadir/sambamba view -p -t 32 -f bam -S ${bamprefix}.sam.wdups -o ${bamprefix}.wdups
           moreexitcode=$?
           exitcode=$(( $exitcode + $moreexitcode ))
           echo `date`
        fi
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem samblaster command failed on $Rone.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit $exitcode;
        fi
        if [ ! -s ${bamprefix}.wdups ]
        then
            MSG="${bamprefix}.wdups aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit 1;
        fi
        echo `date`






########## step 2: sort by coordinate and index on the fly: required for creating an indexed bam, 
##########      which in turn is required for extracting alignments by chromosome
    $novodir/novosort --tmpdir $outputdir --threads $threads --index ${bamprefix}.wdups -o ${bamprefix}.wdups.sorted.bam
    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
         MSG="novosort command failed.  exitcode=$exitcode mergenovo stopped"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" 
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
         cp $qsubfile $AlignOutputLogs/FAILEDjobs/
         exit 1;
    fi
    echo `date`


########### step 2: indexing merged bam file    
#        $samdir/samtools index ${bamprefix}.wdups
#        exitcode=$?
#        if [ $exitcode -ne 0 ]
#        then
#             MSG="samtools index command failed.  exitcode=$exitcode. bwamem_pe_markduplicates stopped "
#             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
#             exit $1;
#        fi
#        echo `date`

        #$samdir/samtools flagstat $tmpfilewdups > $tmpfilewdups.flagstat
        $samdir/samtools view -H ${bamprefix}.wdups.sorted.bam > ${bamprefix}.wdups.sorted.bam.header
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
             MSG="samtools view command failed.  exitcode=$exitcode. bwamem_pe_markduplicates stopped "
             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
             cp $qsubfile $AlignOutputLogs/FAILEDjobs/
             exit $exitcode;
        fi

        echo `date`
        java -Xmx6g -Xms512m -jar $picardir/CollectAlignmentSummaryMetrics.jar \
            INPUT=${bamprefix}.wdups.sorted.bam \
            OUTPUT=${bamprefix}.wdups.sorted.bam.flagstat \
            VALIDATION_STRINGENCY=SILENT

         exitcode=$?
         if [ $exitcode -ne 0 ]
         then
             MSG="collectalignmentsummarymetrics command failed.  exitcode=$exitcode . bwamem_pe_markduplicates stopped "
             #echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
             echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
             cp $qsubfile $AlignOutputLogs/FAILEDjobs/
             exit $exitcode;
         fi
        echo `date`



fi
