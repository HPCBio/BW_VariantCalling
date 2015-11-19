#!/bin/bash
#
#  script to realign and recalibrate the aligned file(s)
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 5 ]
then
    MSG="parameter mismatch."
    echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
    exit 1;
else
        umask 0027
	set -x
	echo `date`
	scriptfile=$0
        runfile=$1
	elog=$2
	olog=$3
	email=$4
        qsubfile=$5
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
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
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        samplefileinfo=$( cat $runfile | grep -w SAMPLEFILENAMES | cut -d '=' -f2 )
        chrindex=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
        indices=$( echo $chrindex | sed 's/^/chr/' | sed 's/:/ chr/g' )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
        cleanupflag=$( cat $runfile | grep -w REMOVETEMPFILES | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUOTHERWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUOTHEREXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for INPUTTYPE=$input_type in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi

        if [ $cleanupflag != "1" -a $cleanupflag != "0" -a $cleanupflag != "YES" -a $cleanupflag != "NO" ]
        then
           MSG="Invalid value for REMOVETEMPFILES=$cleanupflag"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$email""
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
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$email""
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
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 740 $epilogue`
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
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
	    MSG="$samplefileinfo file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        vardir=$outputrootdir/variant
        varlogdir=$outputrootdir/logs/variant
        pipeid=$( cat $outputrootdir/logs/pbs.CONFIGURE )
        vcf_files=""

        if [ ! -d $vardir ]
        then
            mkdir -p $vardir
            mkdir -p $varlogdir
	fi
	if [ ! -d $varlogdir ]
        then
            mkdir -p $varlogdir
        else
	    `rm $varlogdir/*`
        fi

        #generating regions and intervals files in BED format
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

        # main loop
        while read sampledetail
        do
            echo "processing next line of input file..."
            if [ `expr ${#sampledetail}` -gt 0 ]
            then
		echo "processing $sampledetail"
                bam=$( echo $sampledetail | grep ^BAM | cut -d ':' -f2 | cut -d '=' -f2 )

		if [ `expr ${#bam}` -lt 1  ]
		then
		    MSG="parsing of line in SAMPLEFILENAMES failed. variant calling stopped"
                    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                    exit 1;
                else 
                  if [ ! -s $bam ]
                  then
		      MSG="Empty BAM file $bam. variant calling stopped"
                      echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                      exit 1;
                  fi
		fi

		prefix=`basename $bam`
		dirname=`dirname $bam`
		for chr in $indices
		do
		    vcf_files=${vcf_files}":${prefix}"
		    inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/MT/27/' )

                     echo "invoking variant calling..."
		     qsub3=$varlogdir/qsub.vcallgatk.${preffix}.$chr
		     #echo "#PBS -V" > $qsub3
		     echo "#PBS -A $pbsprj" >> $qsub3
		     echo "#PBS -N ${pipeid}_vcallgatk_${preffix}_$chr" >> $qsub3
		     echo "#PBS -l epilogue=$epilogue" >> $qsub3
		     echo "#PBS -l walltime=$pbscpu" >> $qsub3
		     echo "#PBS -l nodes=1:ppn=$thr" >> $qsub3
		     echo "#PBS -o $varlogdir/log.vcallgatk.${preffix}.$chr.ou" >> $qsub3
		     echo "#PBS -e $varlogdir/log.vcallgatk.${preffix}.$chr.in" >> $qsub3
		     echo "#PBS -q $pbsqueue" >> $qsub3
		     echo "#PBS -m ae" >> $qsub3
		     echo "#PBS -M $email" >> $qsub3
                     echo "source /opt/modules/default/init/bash" >> $qsub3
		     echo "$scriptdir/vcallgatk.sh $vardir $dirname $prefix $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.${preffix}.$chr.in $varlogdir/log.vcallgatk.${preffix}.$chr.ou $email $varlogdir/qsub.vcallgatk.${preffix}.$chr" >> $qsub3
		     #`chmod a+r $qsub3`
		     vcalljobid=`qsub $qsub3`
		     echo $vcalljobid >> $outputrootdir/logs/pbs.VARCALL
		done
            else
                 echo "skipping empty line"
            fi
	    echo `date`
     done < $samplefileinfo

     #`chmod -R 770 $vardir/`
     listjobids=$( cat $outputrootdir/logs/pbs.VARCALL | sed "s/\.[a-z]*//g" | tr "\n" ":" )
     lastjobid=""
     mergevcfjobid=""
     cleanjobid=""

     # producing a single VCF file only if analysis=REALIGNMENT
     # and SKIPVCALL=NO

     if [ $skipvcall == "NO" ]
     then
	 qsub1=$varlogdir/qsub.mergevcf
	 echo "#PBS -V" > $qsub1
	 echo "#PBS -A $pbsprj" >> $qsub1
	 echo "#PBS -N ${pipeid}_mergevcf" >> $qsub1
	 echo "#PBS -l epilogue=$epilogue" >> $qsub1
	 echo "#PBS -l walltime=$pbscpu" >> $qsub1
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub1
	 echo "#PBS -o $varlogdir/log.mergevcf.ou" >> $qsub1
	 echo "#PBS -e $varlogdir/log.mergevcf.in" >> $qsub1
	 echo "#PBS -q $pbsqueue" >> $qsub1
	 echo "#PBS -m ae" >> $qsub1
	 echo "#PBS -M $email" >> $qsub1
	 echo "#PBS -W depend=afterok:$listjobids" >> $qsub1
	 echo "$scriptdir/mergevcf.sh $vardir $vcf_files $tabixdir $vcftoolsdir $varlogdir/log.mergevcf.in $varlogdir/log.mergevcf.ou $email $varlogdir/qsub.mergevcf"  >> $qsub1
	 #`chmod a+r $qsub1`
	 mergevcfjobid=`qsub $qsub1`
	 echo $mergevcfjobid >> $outputrootdir/logs/pbs.MERGEVCF
     fi

     # removing temporary files

     if [ $cleanupflag == "YES" ]
     then
	 qsub6=$outputrootdir/logs/qsub.cleanup.vcall
	 echo "#PBS -V" > $qsub6
	 echo "#PBS -A $pbsprj" >> $qsub6
	 echo "#PBS -N ${pipeid}_cleanup_vcall" >> $qsub6
	 echo "#PBS -l epilogue=$epilogue" >> $qsub6
	 echo "#PBS -l walltime=$pbscpu" >> $qsub6
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub6
	 echo "#PBS -o $outputrootdir/logs/log.cleanup.vcall.ou" >> $qsub6
	 echo "#PBS -e $outputrootdir/logs/log.cleanup.vcall.in" >> $qsub6
	 echo "#PBS -q $pbsqueue" >> $qsub6
	 echo "#PBS -m ae" >> $qsub6
	 echo "#PBS -M $email" >> $qsub6
         if [ `expr ${#mergevcfjobid}` -lt 1 ]
         then
	     echo "#PBS -W depend=afterok:$listjobids" >> $qsub6
         else
	     echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub6
         fi
	 echo "$scriptdir/cleanup.sh $outputrootdir $analysis $outputrootdir/logs/log.cleanup.vcall.in $outputrootdir/logs/log.cleanup.vcall.ou $email $outputrootdir/logs/qsub.cleanup.vcall"  >> $qsub6
	 #`chmod a+r $qsub6`
	 cleanjobid=`qsub $qsub6`
	 echo $cleanjobid >> $outputrootdir/logs/pbs.CLEANUP
     fi

     # producing summary report
     `sleep 30s`
     qsub4=$outputrootdir/logs/qsub.summary.vcall.allok
     echo "#PBS -V" > $qsub4
     echo "#PBS -A $pbsprj" >> $qsub4
     echo "#PBS -N ${pipeid}_summaryok" >> $qsub4
     echo "#PBS -l epilogue=$epilogue" >> $qsub4
     echo "#PBS -l walltime=$pbscpu" >> $qsub4
     echo "#PBS -l nodes=1:ppn=1" >> $qsub4
     echo "#PBS -o $outputrootdir/logs/log.summary.vcall.ou" >> $qsub4
     echo "#PBS -e $outputrootdir/logs/log.summary.vcall.in" >> $qsub4
     echo "#PBS -q $pbsqueue" >> $qsub4
     echo "#PBS -m ae" >> $qsub4
     echo "#PBS -M $email" >> $qsub4
     if [ `expr ${#cleanjobid}` -gt 0 ]
     then
	 echo "#PBS -W depend=afterok:$cleanjobid" >> $qsub4
     else
         if [ `expr ${#mergevcfjobid}` -lt 1 ]
         then
	     echo "#PBS -W depend=afterok:$listjobids" >> $qsub4
         else
	     echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub4
         fi         
     fi
     echo "$scriptdir/summary.sh $runfile $email exitok"  >> $qsub4
     #`chmod a+r $qsub4`
     lastjobid=`qsub $qsub4`
     echo $lastjobid >> $outputrootdir/logs/pbs.SUMMARY

     if [ `expr ${#lastjobid}` -lt 1 ]
     then
         echo "at least one job aborted"
	 qsub5=$outputrootdir/logs/qsub.summary.vcall.afterany
	 echo "#PBS -V" > $qsub5
	 echo "#PBS -A $pbsprj" >> $qsub5
	 echo "#PBS -N ${pipeid}_summary_afterany" >> $qsub5
	 echo "#PBS -l epilogue=$epilogue" >> $qsub5
	 echo "#PBS -l walltime=$pbscpu" >> $qsub5
	 echo "#PBS -l nodes=1:ppn=1" >> $qsub5
	 echo "#PBS -o $outputrootdir/logs/log.summary.vcall.afterany.ou" >> $qsub5
	 echo "#PBS -e $outputrootdir/logs/log.summary.vcall.afterany.in" >> $qsub5
	 echo "#PBS -q $pbsqueue" >> $qsub5
	 echo "#PBS -m ae" >> $qsub5
	 echo "#PBS -M $email" >> $qsub5
	 echo "#PBS -W depend=afterany:$listjobids" >> $qsub5
	 echo "$scriptdir/summary.sh $runfile $email exitnotok"  >> $qsub5
	 #`chmod a+r $qsub5`
	 badjobid=`qsub $qsub5`
	 echo $badjobid >> $outputrootdir/logs/pbs.SUMMARY
     fi
     #`chmod -R 770 $outputroordir/logs`
     echo `date`
fi

echo -e "############## now making PBS logs group readable"

find $outputrootdir -name logs -type d | awk 'print "chmod -R g=rwx "$1}' | sh -x
