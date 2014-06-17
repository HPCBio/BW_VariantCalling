#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#
#  script to realign and recalibrate the aligned file(s)
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 9 ]
then
    MSG="parameter mismatch."
    echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
    #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
    exit 1;
else
	set -x
	echo `date`
	scriptfile=$0
        realigndir=$1
        realignlogdir=$2
        aligndir=$3
        runfile=$4
        flag=$5
	elog=$6
	olog=$7
	email=$8
        qsubfile=$9
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        tabixdir=$( cat $runfile | grep -w TABIXDIR | cut -d '=' -f2 )
        vcftoolsdir=$( cat $runfile | grep -w VCFTOOLSDIR | cut -d '=' -f2 )
        dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        kgenome=$( cat $runfile | grep -w KGENOME | cut -d '=' -f2 )
        targetkit=$( cat $runfile | grep -w ONTARGET | cut -d '=' -f2 )
        realignparams=$( cat $runfile | grep -w REALIGNPARMS | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 )
        chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
        indices=$( echo $chrindex | sed 's/^/chr/' | sed 's/:/ chr/g' )
        sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
        sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
        sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
        skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
        cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )

        if [ $type == "GENOME" -o $type == "WHOLE_GENOME" -o $type == "WHOLEGENOME" -o $type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $type == "EXOME" -o $type == "WHOLE_EXOME" -o $type == "WHOLEEXOME" -o $type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for TYPE=$type in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		#echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi

        if [ $cleanupflag != "1" -a $cleanupflag != "0" -a $cleanupflag != "YES" -a $cleanupflag != "NO" ]
        then
           MSG="Invalid value for REMOVETEMPFILES=$cleanupflag"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit 1;
        else
            if [ $cleanupflag == "1" ]
            then
                $cleanupflag="YES"
            fi
            if [ $cleanupflag == "0" ]
            then
                $cleanupflag="NO"
            fi
        fi

        if [ $skipvcall != "1" -a $skipvcall != "0" -a $skipvcall != "YES" -a $skipvcall != "NO" ]
        then
           MSG="Invalid value for SKIPVCALL=$skipvcall"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit 1;
        else
            if [ $skipvcall == "1" ]
            then
                $skipvcall="YES"
            fi
            if [ $skipvcall == "0" ]
            then
                $skipvcall="NO"
            fi
        fi

        if [ -z $epilogue ]
        then
           MSG="Value for EPILOGUE must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
        fi
        if [ -z $javamodule ]
        then
           MSG="Value for JAVAMODULE must be specified in configuration file"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        if [ ! -d $picardir ]
        then
	    MSG="$picardir picard directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
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

        if [ -s $refdir/$dbSNP ]
        then
	    realparms="-known:$refdir/$dbSNP"
            recalparms="--knownSites:$refdir/$dbSNP"
        fi
        if [ -s $refdir/$kgenome ]
        then
	    realparms=$realparms":-known:$refdir/$kgenome"
	    recalparms=$recalparms":--knownSites:$refdir/$kgenome"
        fi
        if [ ! -d $realignlogdir ]
        then
	    MSG="$realignlogdir realignlog directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $realigndir ]
        then
	    MSG="$realigndir realign directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        # get aligned bam files when alignment is done inhouse
        if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ]
        then
            echo "inhouse alignment. gathering aligned bam files"
            cd $aligndir
            allfiles=`find ./ -name "*.wdups.sorted.bam"`
	    if [ `expr ${#allfiles}` -lt 1 ]
	    then
		MSG="No bam file(s) found to perform realign-recalibrate at $aligndir"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi
            listfiles=""
	    sep=" "
            for file in $allfiles
            do
		newname=$aligndir/$file
		listfiles=$newname${sep}${listfiles}
            done
	else
            if [ $analysis == "REALIGN_ONLY" -o $analysis == "REALIGNONLY" ]
            then
		echo "3rd param, aligndir,  has the list of aligned bam files"
		listfiles=$( echo $aligndir | tr ":" " " )
            else
		MSG="Invalid value for analysis=$analysis"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
	fi

        pipeid=$( cat $outputrootdir/logs/MAINpbs )
        outputdir=$realigndir

        #if [ $skipvcall == "NO" ]
        #then
        #    vardir=$outputrootdir/variant
        #    varlogdir=$outputrootdir/logs/variant
        #    if [ ! -d $vardir ]
        #    then
	#	mkdir -p $vardir
	#	mkdir -p $varlogdir
	#    fi
	#    if [ ! -d $varlogdir ]
        #    then
	#	mkdir -p $varlogdir
        #    else
	#	`rm $varlogdir/*`
        #    fi
        #fi

        #generating regions and intervals files in BED format
        echo `date`
        for chr in $indices
        do
            i=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
            if [ -d $targetkit ]
            then
		if [ `cat $targetkit/${chr}.bed | wc -l` -gt 0 ]
                then
                    region[$i]="-L:$targetkit/${chr}.bed"
                else
		    region[$i]="-L:$chr"
                fi
            else
		region[$i]="-L:$chr"
            fi
        done
        echo `date`

        igv_files=""
        vcf_files=""
        sep=":"

        ####################
        # main loops start here
        # outer loop by chromosome; inner loop by sample
        #####################

        cd $outputdir
	for chr in $indices
        do
            echo "generating real-recal calls fro chr=${chr} ..."
	    echo `date`
	    inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )

            for bam in $listfiles
            do
                bamfile=`basename $bam`
		sample=`basename $bamfile .bam`
		sID=$sample
		sPU=$sample
		sSM=$sample
		RGparms=$( echo "RGID=${sID}:RGLB=${sLB}:RGPU=${sPU}:RGSM=${sSM}:RGPL=${sPL}:RGCN=${sCN}" )

                qsub1=$realignlogdir/qsub.sort.$bamfile.$chr
                echo "#PBS -V" > $qsub1
                echo "#PBS -A $pbsprj" >> $qsub1
                echo "#PBS -N ${pipeid}_sort_${bamfile}_$chr" >> $qsub1
                echo "#PBS -l epilogue=$epilogue" >> $qsub1
		echo "#PBS -l walltime=$pbscpu" >> $qsub1
		echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
		echo "#PBS -o $realignlogdir/log.sort.$bamfile.$chr.ou" >> $qsub1
		echo "#PBS -e $realignlogdir/log.sort.$bamfile.$chr.in" >> $qsub1
                echo "#PBS -q $pbsqueue" >> $qsub1
                echo "#PBS -m ae" >> $qsub1
                echo "#PBS -M $email" >> $qsub1
                echo "aprun -n 1 -d $thr $scriptdir/sortnode.sh $picardir $samdir $javamodule $outputdir $bam $bamfile.$chr.bam $bamfile.$chr.sorted.bam $RGparms $flag $chr $realignlogdir/log.sort.$bamfile.$chr.in $realignlogdir/log.sort.$bamfile.$chr.ou $email $realignlogdir/qsub.sort.$bamfile.$chr" >> $qsub1
                `chmod a+r $qsub1`
                sortjobid=`qsub $qsub1`
                # new line to avoid hiccup
                #`qhold -h u $sortjobid`
                echo $sortjobid >> $outputrootdir/logs/REALSORTEDpbs.$bamfile
                echo $sortjobid >> $outputrootdir/logs/REALSORTED.$bamfile_$chr
		chrinfiles[$inx]=${chrinfiles[$inx]}":-I:$outputdir/$bamfile.$chr.sorted.bam"
		chrinputfiles[$inx]=${chrinputfiles[$inx]}":INPUT=$outputdir/$bamfile.$chr.sorted.bam"
	    done
	    echo `date`


            sortid=$( cat $outputrootdir/logs/REALSORTED.$bamfile_$chr | sed "s/\..*//g" | tr "\n" ":" )
            outputfile=$chr.realrecal.$bamfile.output.bam
            igv_files=${igv_files}":INPUT=${outputfile}"
	    echo "realign-recalibrate for interval:$chr..."
	    qsub2=$realignlogdir/qsub.realrecal.$bamfile.$chr
	    echo "#PBS -V" > $qsub2
	    echo "#PBS -A $pbsprj" >> $qsub2
	    echo "#PBS -N ${pipeid}_realrecal_$bamfile.$chr" >> $qsub2
	    echo "#PBS -l epilogue=$epilogue" >> $qsub2
	    echo "#PBS -l walltime=$pbscpu" >> $qsub2
	    echo "#PBS -l nodes=1:ppn=$thr" >> $qsub2
	    echo "#PBS -o $realignlogdir/log.realrecal.$bamfile.$chr.ou" >> $qsub2
	    echo "#PBS -e $realignlogdir/log.realrecal.$bamfile.$chr.in" >> $qsub2
	    echo "#PBS -q $pbsqueue" >> $qsub2
	    echo "#PBS -m ae" >> $qsub2
	    echo "#PBS -M $email" >> $qsub2
	    echo "#PBS -W depend=afterok:$sortid" >> $qsub2
	    echo "aprun -n 1 -d $thr $scriptdir/realrecalold.sh $outputdir $outputfile $chr ${chrinfiles[$inx]} ${chrinputfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $realignlogdir/log.realrecal.$bamfile.$chr.in $realignlogdir/log.realrecal.$bamfile.$chr.ou $email $realignlogdir/qsub.realrecal.$bamfile.$chr" >> $qsub2
	    `chmod a+r $qsub2`
	    recaljobid=`qsub $qsub2`
            # new line to avoid hiccup
            #`qhold -h u $recaljobid`
	    echo $recaljobid >> $outputrootdir/logs/REALRECALpbs.$bamfile


            if [ $skipvcall == "NO" ]
            then
                vardir=$( echo $outputdir | sed 's/realign/variant/' )
                varlogdir=$( echo $realignlogdir | sed 's/realign/variant/' )
  
		vcf_files=${vcf_files}":${outputfile}"
                echo "variant calling call.."
		qsub3=$varlogdir/qsub.vcallgatk.$bamfile.$chr
		echo "#PBS -V" > $qsub3
		echo "#PBS -A $pbsprj" >> $qsub3
		echo "#PBS -N ${pipeid}_vcall_$chr" >> $qsub3
		echo "#PBS -l epilogue=$epilogue" >> $qsub3
		echo "#PBS -l walltime=$pbscpu" >> $qsub3
		echo "#PBS -l nodes=1:ppn=$thr" >> $qsub3
		echo "#PBS -o $varlogdir/log.vcallgatk.$bamfile.$chr.ou" >> $qsub3
		echo "#PBS -e $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $qsub3
		echo "#PBS -q $pbsqueue" >> $qsub3
		echo "#PBS -m ae" >> $qsub3
		echo "#PBS -M $email" >> $qsub3
		echo "#PBS -W depend=afterok:$recaljobid" >> $qsub3
		echo "aprun -n 1 -d $thr $scriptdir/vcallgatk.sh $vardir $outputdir $outputfile $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$bamfile.$chr.in $varlogdir/log.vcallgatk.$bamfile.$chr.ou $email $varlogdir/qsub.vcallgatk.$bamfile.$chr" >> $qsub3
		`chmod a+r $qsub3`
		vcalljobid=`qsub $qsub3`
		echo $vcalljobid >> $outputrootdir/logs/VCALLGATKpbs.$bamfile
            else
                     echo "variant calling will not be run"
            fi
	    echo `date`
            echo "bottom of the loop"
     done
     echo `date`

     echo "clean up and produce summary table"
     if [ $skipvcall == "NO" ]
     then
	 listjobids=$( cat $outputrootdir/logs/VCALLGATKpbs.$bamfile | sed "s/\..*//g" | tr "\n" ":" )
     else
	 listjobids=$( cat $outputrootdir/logs/REALRECALpbs.$bamfile | sed "s/\..*//g" | tr "\n" ":" )
     fi

     # preparing output files for realignment and/or variant calling

     qsub5=$outputrootdir/logs/qsub.igvbam.$bamfile
     echo "#PBS -V" > $qsub5
     echo "#PBS -A $pbsprj" >> $qsub5
     echo "#PBS -N ${pipeid}_igvbam.$bamfile" >> $qsub5
     echo "#PBS -l epilogue=$epilogue" >> $qsub5
     echo "#PBS -l walltime=$pbscpu" >> $qsub5
     echo "#PBS -l nodes=1:ppn=$thr" >> $qsub5
     echo "#PBS -o $outputrootdir/logs/log.igvbam.$bamfile.ou" >> $qsub5
     echo "#PBS -e $outputrootdir/logs/log.igvbam.$bamfile.in" >> $qsub5
     echo "#PBS -q $pbsqueue" >> $qsub5
     echo "#PBS -m ae" >> $qsub5
     echo "#PBS -M $email" >> $qsub5
     echo "#PBS -W depend=afterok:$listjobids" >> $qsub5
     echo "aprun -n 1 -d $thr $scriptdir/igvbam.sh $outputdir $igv_files $runfile $outputrootdir/logs/log.igvbam.$bamfile.in $outputroot/logs/log.igvbam.$bamfile.ou $email $outputrootdir/logs/qsub.igvbam.$bamfile"  >> $qsub5
     `chmod a+r $qsub5`
     igvjobid=`qsub $qsub5`
     echo $igvjobid >> $outputrootdir/logs/IGVBAMpbs.$bamfile

     if [ $cleanupflag == "YES" ]
     then
	 qsub6=$outputrootdir/logs/qsub.cleanup.realn.$bamfile
	 echo "#PBS -V" > $qsub6
	 echo "#PBS -A $pbsprj" >> $qsub6
	 echo "#PBS -N ${pipeid}_cleanup_realn.$bamfile" >> $qsub6
	 echo "#PBS -l epilogue=$epilogue" >> $qsub6
	 echo "#PBS -l walltime=$pbscpu" >> $qsub6
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub6
	 echo "#PBS -o $outputrootdir/logs/log.cleanup.realn.$bamfile.ou" >> $qsub6
	 echo "#PBS -e $outputrootdir/logs/log.cleanup.realn.$bamfile.in" >> $qsub6
	 echo "#PBS -q $pbsqueue" >> $qsub6
	 echo "#PBS -m ae" >> $qsub6
	 echo "#PBS -M $email" >> $qsub6
	 echo "#PBS -W depend=afterok:$igvjobid" >> $qsub6
	 echo "$scriptdir/cleanup.sh $outputrootdir $analysis $outputrootdir/logs/log.cleanup.realn.$bamfile.in $outputrootdir/logs/log.cleanup.realn.$bamfile.ou $email $outputrootdir/logs/qsub.cleanup.realn.$bamfile"  >> $qsub6
	 `chmod a+r $qsub6`
	 cleanjobid=`qsub $qsub6`
	 echo $cleanjobid >> $outputrootdir/logs/CLEANUPpbs.$bamfile
     fi

     if [ $skipvcall == "NO" ]
     then
	 qsub1=$varlogdir/qsub.mergevcf.$bamfile
	 echo "#PBS -V" > $qsub1
	 echo "#PBS -A $pbsprj" >> $qsub1
	 echo "#PBS -N ${pipeid}_mergevcf.$bamfile" >> $qsub1
	 echo "#PBS -l epilogue=$epilogue" >> $qsub1
	 echo "#PBS -l walltime=$pbscpu" >> $qsub1
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	 echo "#PBS -o $varlogdir/log.mergevcf.$bamfile.ou" >> $qsub1
	 echo "#PBS -e $varlogdir/log.mergevcf.$bamfile.in" >> $qsub1
	 echo "#PBS -q $pbsqueue" >> $qsub1
	 echo "#PBS -m ae" >> $qsub1
	 echo "#PBS -M $email" >> $qsub1
	 echo "#PBS -W depend=afterok:$listjobids" >> $qsub1
	 echo "$scriptdir/mergevcf.sh $vardir $vcf_files $tabixdir $vcftoolsdir $varlogdir/log.mergevcf.$bamfile.in $varlogdir/log.mergevcf.$bamfile.ou $email $varlogdir/qsub.mergevcf.$bamfile"  >> $qsub1
	 `chmod a+r $qsub1`
	 mergevcfjobid=`qsub $qsub1`
	 echo $mergevcfjobid >> $outputrootdir/logs/MERGEVCFpbs

	 if [ $cleanupflag == "YES" ]
	 then
	     qsub6=$outputrootdir/logs/qsub.cleanup.vcall.$bamfile
	     echo "#PBS -V" > $qsub6
	     echo "#PBS -A $pbsprj" >> $qsub6
	     echo "#PBS -N ${pipeid}_cleanup_vcall.$bamfile" >> $qsub6
	     echo "#PBS -l epilogue=$epilogue" >> $qsub6
	     echo "#PBS -l walltime=$pbscpu" >> $qsub6
	     echo "#PBS -l nodes=1:ppn=1" >> $qsub6
	     echo "#PBS -o $outputrootdir/logs/log.cleanup.vcall.$bamfile.ou" >> $qsub6
	     echo "#PBS -e $outputrootdir/logs/log.cleanup.vcall.$bamfile.in" >> $qsub6
	     echo "#PBS -q $pbsqueue" >> $qsub6
	     echo "#PBS -m ae" >> $qsub6
	     echo "#PBS -M $email" >> $qsub6
	     echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub6
	     echo "$scriptdir/cleanup.sh $outputrootdir VCALL_ONLY $outputrootdir/logs/log.cleanup.vcall.$bamfile.in $outputrootdir/logs/log.cleanup.vcall.$bamfile.ou $email $outputrootdir/logs/qsub.cleanup.vcall.$bamfile"  >> $qsub6
	     `chmod a+r $qsub6`
	     cleanjobid=`qsub $qsub6`
	     echo $cleanjobid >> $outputrootdir/logs/CLEANUPpbs.$bamfile
	 fi
     fi

     # new lines to avoid hiccups and missing jobs in summary

     heldjobs_realsorted=$( cat $outputrootdir/logs/REALSORTEDpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
     heldjobs_realrecal=$( cat $outputrootdir/logs/REALRECALpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
     heldjobs_vcall=$( cat $outputrootdir/logs/VCALLGATKpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )
     cleanupjobs=$( cat $outputrootdir/logs/CLEANUPpbs.$bamfile | sed "s/\..*//g" | tr "\n" " " )

     #`qrls -h u $heldjobs_realsorted`
     #`qrls -h u $heldjobs_realrecal`
     #`qrls -h u $heldjobs_vcall`
     `sleep 10s`

     qsub4=$outputrootdir/logs/qsub.summary.realn.allok.$bamfile
     echo "#PBS -V" > $qsub4
     echo "#PBS -A $pbsprj" >> $qsub4
     echo "#PBS -N ${pipeid}_summaryok.$bamfile" >> $qsub4
     echo "#PBS -l epilogue=$epilogue" >> $qsub4
     echo "#PBS -l walltime=$pbscpu" >> $qsub4
     echo "#PBS -l nodes=1:ppn=1" >> $qsub4
     echo "#PBS -o $outputrootdir/logs/log.summary.realn.$bamfile.ou" >> $qsub4
     echo "#PBS -e $outputrootdir/logs/log.summary.realn.$bamfile.in" >> $qsub4
     echo "#PBS -q $pbsqueue" >> $qsub4
     echo "#PBS -m ae" >> $qsub4
     echo "#PBS -M $email" >> $qsub4
     if [ `expr ${#cleanupjobs}` -gt 0 ]
     then 
	 echo "#PBS -W depend=afterok:$cleanupjobs" >> $qsub4
     else
         if [ $skipvcall == "YES" ]
         then
	     echo "#PBS -W depend=afterok:$igvjobid" >> $qsub4
         else
	     echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub4
         fi
     fi
     echo "$scriptdir/summary.sh $outputrootdir $email exitok"  >> $qsub4
     `chmod a+r $qsub4`
     lastjobid=`qsub $qsub4`
     echo $lastjobid >> $outputrootdir/logs/SUMMARYpbs.$bamfile

     if [ `expr ${#lastjobid}` -lt 1 ]
     then
         echo "at least one job aborted"
	 qsub5=$outputrootdir/logs/qsub.summary.realn.afterany.$bamfile
	 echo "#PBS -V" > $qsub5
	 echo "#PBS -A $pbsprj" >> $qsub5
	 echo "#PBS -N ${pipeid}_summary_afterany.$bamfile" >> $qsub5
	 echo "#PBS -l epilogue=$epilogue" >> $qsub5
	 echo "#PBS -l walltime=$pbscpu" >> $qsub5
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub5
	 echo "#PBS -o $outputrootdir/logs/log.summary.realn.afterany.$bamfile.ou" >> $qsub5
	 echo "#PBS -e $outputrootdir/logs/log.summary.realn.afterany.$bamfile.in" >> $qsub5
	 echo "#PBS -q $pbsqueue" >> $qsub5
	 echo "#PBS -m ae" >> $qsub5
	 echo "#PBS -M $email" >> $qsub5
	 echo "#PBS -W depend=afterany:$listjobids" >> $qsub5
	 echo "$scriptdir/summary.sh $outputtopdir $email exitnotok"  >> $qsub5
	 `chmod a+r $qsub5`
	 badjobid=`qsub $qsub5`
	 echo $badjobid >> $outputrootdir/logs/SUMMARYpbs.$bamfile
     fi
     `chmod -R 770 $outputroordir/logs`
     echo `date`
fi
