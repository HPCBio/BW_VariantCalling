#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#	
#  script to realign and recalibrate the aligned file(s) 
########################################################
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 13 ];
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x

	echo `date`
	scriptfile=$0
        realrecaldir=$1
        outputfile=$2	
        chr=$3
        chrinfiles=$4
        #chrinputfiles=$5 # not used anymore?
        region=$5
        realparms=$6
        recalparms=$7
        runfile=$8
        flag=${9}
	elog=${10}
	olog=${11}
	email=${12}
        qsubfile=${13}
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
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
        dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
        targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
        realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
	outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
        #thr=`expr $threads "-" 1`
        thr=`expr $threads "/" 2`

        # cleaning up the lists
        chrinfiles=$( echo $chrinfiles | tr ":" " " )
#        chrinputfiles=$( echo $chrinputfiles | tr ":" " " ) # not used anymore?
        region=$( echo $region | tr ":" " " )
        realparms=$( echo $realparms | tr ":" " " )
        recalparms=$( echo $recalparms | tr ":" " " )

        if [ ! -d $picardir ]
        then
	    MSG="$picardir picard directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge
"mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""	    exit 1;
        fi
        if [ ! -d $samdir ]
        then
	    MSG="$samdir samtools directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ -z $javadir ]
        then
	    MSG="Value for JAVADIR must be specified in configuration file"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        #else
            #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
        #    `module load $javamodule`
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        inputdir=$outputrootdir/realign
        if [ ! -d $inputdir ]
        then
	    MSG="$inputdir realign directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        cd $realrecaldir

        if [ $flag == 1 ]
        then
		echo "realign then recalibrate"
		echo "realigning..."
		
            	$memprof $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $chrinfiles \
		    -T RealignerTargetCreator \
		    -o realign.$chr.bam.list $realparms $region

		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="realignertargetcreator command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`
	
		if [ ! -s realign.$chr.bam.list ]
		then
		    MSG="realign.$chr.bam.list realignertargetcreator file not created. realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit 1;
		fi
		echo `date`


		$memprof $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $chrinfiles \
		    -T IndelRealigner \
                    -L $chr \
		    -o realign.$chr.realigned.bam \
		    -targetIntervals realign.$chr.bam.list $realignparams $realparms

		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="indelrealigner command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`

		if [ ! -s realign.$chr.realigned.bam ]
		then
		    MSG="realigned.bam indelrealigner file not created. realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit 1;
                fi
		$memprof mv realign.$chr.realigned.bai realign.$chr.realigned.bam.bai
		$memprof cp realign.$chr.realigned.bam realign.$chr.real.cleaned.bam
		$memprof cp realign.$chr.realigned.bam.bai realign.$chr.real.cleaned.bam.bai
                $memprof $samdir/samtools flagstat realign.$chr.real.cleaned.bam > realign.$chr.real.cleaned.flagstat

		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="samtools flagstat command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		

		echo "recalibrating..."

                if [ "$recalibrator" == "BQSR" ]
                then
                   $memprov $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
                       -R $refdir/$ref \
                       $recalparms \
                       -I realign.$chr.real.cleaned.bam \
                       $region  \
                       -T BaseRecalibrator \
                       --out recal.$chr.recal_report.grp \
                       -nct 32

                   exitcode=$?
                   if [ $exitcode -ne 0 ]
                   then
                       MSG="BQSR command failed exitcode=$exitcode  recalibration stopped"
                       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                       exit $exitcode;
                   fi
                   echo `date`

                   if [ ! -s recal.$chr.recal_report.grp ]
                   then
                       MSG="recal.$chr.recal_report.grp recalibration report file not created, realignment-recalibration stopped"
                       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                       exit 1;
                   fi
                   echo `date`


                   $memprov $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
                       -R $refdir/$ref \
                       -I realign.$chr.real.cleaned.bam \
                       $region  \
                       -T PrintReads \
                       -BQSR recal.$chr.recal_report.grp \
                       --out recal.$chr.real.recal.bam \
                       -nct 32

                   exitcode=$?
                   if [ $exitcode -ne 0 ]
                   then
                       MSG="BQSR PrintReads command failed exitcode=$exitcode  recalibration stopped"
                       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                       exit $exitcode;
                   fi
                   echo `date`


                   if [ ! -s recal.$chr.real.recal.bam ]
                   then
                       MSG="recal.$chr.real.recal.bam recalibrated file not created, realignment-recalibration stopped"
                       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                       exit 1;
                   fi
                   $memprof cp recal.$chr.real.recal.bam $outputfile
                   $memprof cp recal.$chr.real.recal.bai $outputfile.bai

                else	
   		   $memprof $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		       -R $refdir/$ref \
		       $recalparms \
		       -I realign.$chr.real.cleaned.bam \
		       $region  \
		       -T CountCovariates \
		       -cov ReadGroupCovariate \
		       -cov QualityScoreCovariate \
		       -cov CycleCovariate \
		       -cov DinucCovariate \
		       -recalFile recal.$chr.recal_data.csv 

		   exitcode=$?
		   if [ $exitcode -ne 0 ]
		   then
		       MSG="countcovariates command failed exitcode=$exitcode  realignment-recalibration stopped"
		       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       exit $exitcode;
		   fi
		   echo `date`		

		   if [ ! -s recal.$chr.recal_data.csv ]
		   then
		       MSG="recal.$chr.recal_data.csv countcovariates file not created realignment-recalibration stopped"
		       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       exit 1;
		   fi
		   echo `date`

		   $memprof $javadir/java -Xmx1024m -Xms1024m  -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		       -R $refdir/$ref \
		       -L $chr \
		       -I realign.$chr.real.cleaned.bam \
		       -T TableRecalibration \
		       --out recal.$chr.real.recal.bam \
		       -recalFile recal.$chr.recal_data.csv 

		   exitcode=$?
		   if [ $exitcode -ne 0 ]
		   then
		       MSG="tablerecalibration command failed exitcode=$exitcode  realignment-recalibration stopped"
		       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       exit $exitcode;
		   fi
		   echo `date`		


		   if [ ! -s recal.$chr.real.recal.bam ]
		   then
		       MSG="$chr.real.recal.bam tablerecalibration file not created realignment-recalibration stopped"
		       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		       exit 1;
		   fi
		   $memprof cp recal.$chr.real.recal.bam $outputfile
		   $memprof cp recal.$chr.real.recal.bai $outputfile.bai
                fi # end choosing between BQSR and CountCovariates/TableRecalibration

		$memprof $samdir/samtools flagstat $outputfile > $outputfile.flagstat
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="samtools flagstat command failed exitcode=$exitcode "
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		

        else
		echo "recalibrate then realign. not common practice"
		echo "recalibrating"

		$memprof $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    $recalparms \
		    $chrinfiles \
		    $region \
		    -T CountCovariates \
		    -cov ReadGroupCovariate \
		    -cov QualityScoreCovariate \
		    -cov CycleCovariate \
		    -cov DinucCovariate \
		    -recalFile recal.$chr.recal_data.csv 

		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="countcovariates command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		
	
		if [ ! -s recal.$chr.recal_data.csv ]
		then
		    MSG="recal.$chr.recal_data.csv countcovariates file not created realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
		fi
	
		echo `date`

		$memprof $javadir/java -Xmx1024m -Xms1024m  -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -L $chr \
		    $chrinfiles \
		    -T TableRecalibration \
		    --out recal.$chr.recal.bam \
		    -recalFile recal.$chr.recal_data.csv 

		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="tablerecalibration command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		

		if [ ! -s recal.$chr.recal.bam ]
		then
		    MSG="recal.$chr.recal.bam tablerecalibration file not created realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
                fi
                $memprof cp recal.$chr.recal.bam recal.$chr.recal.cleaned.bam
                $memprof cp recal.$chr.recal.bai recal.$chr.recal.cleaned.bam.bai
		$memprof $samdir/samtools flagstat recal.$chr.recal.cleaned.bam > recal.$chr.recal.cleaned.flagstat
		
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="samtools flagstat command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date` 		

		echo "realigning"
		$memprof $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -I recal.$chr.recal.bam \
		    -o realign.$chr.recal.bam.list \
		    -T RealignerTargetCreator $realparms $region
	
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="realignertargetcreator command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		
		if [ ! -s realign.$chr.recal.bam.list ]
		then
		    MSG="realign.$chr.recal.bam.list realigntargetcreator file not created. realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                    exit 1;
		fi
		echo `date`

		$memprof $javadir/java -Xmx1024m -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
		    -R $refdir/$ref \
		    -I recal.$chr.recal.bam \
		    -T IndelRealigner \
                    -L $chr \
		    -targetIntervals realign.$chr.recal.bam.list $realignparms $realparms \
		    -o realign.$chr.recal.realigned.bam

		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="indelrealigner command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		
		if [ ! -s realign.$chr.recal.realigned.bam ]
		then
		    MSG="realign.$chr.recal.realigned.bam indelrealigner file not created. realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit 1;
		fi
		$memprof cp realign.$chr.recal.realigned.bam $outputfile
		$memprof cp realign.$chr.recal.realigned.bam.bai $outputfile.bai
		$memprof $samdir/samtools flagstat $outputfile > $outputfile.flagstat
		exitcode=$?
		if [ $exitcode -ne 0 ]
		then
		    MSG="samtools flagstat command failed exitcode=$exitcode  realignment-recalibration stopped"
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		    exit $exitcode;
		fi
		echo `date`		
        fi

fi
