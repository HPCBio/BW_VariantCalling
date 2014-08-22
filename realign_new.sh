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
   schedule=$( cat $runfile | grep -w SCHEDULE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
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
   #indices=$( echo $chrindex | sed 's/^/chr/' | sed 's/:/ chr/g' )
   indices=$( echo $chrindex | sed 's/:/ /g' )
   sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
   sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
   sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
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



   #################make sure the checks work everywhere: != vs ne
        if [ $schedule ne "ANISIMOV" -a $schedule ne "QSUB" -a $schedule ne "SERIAL" ]
        then
           MSG="Invalid value for SCHEDULE=$schedule"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
           exit 1; 
   elif [ $schedule eq "ANISIMOV" ]
        then
           truncate -s 0 $realignlogdir/$realrecal.AnisimovJoblist
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
   #   mkdir -p $vardir
   #   mkdir -p $varlogdir
   #    fi
   #    if [ ! -d $varlogdir ]
        #    then
   #   mkdir -p $varlogdir
        #    else
   #   `rm $varlogdir/*`
        #    fi
        #fi

        #generating regions and intervals files in BED format
        echo `date`
   i=1
        for chr in $indices
        do
            #i=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )
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
       (($i++)) ########### make sure this works
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
   inx=1
   for chr in $indices
        do
            echo "generating real-recal calls fro chr=${chr} ..."
       echo `date`
       #inx=$( echo $chr | sed 's/chr//' | sed 's/X/25/' | sed 's/Y/26/' | sed 's/M/27/' )

            for bam in $listfiles
            do
                bamfile=`basename $bam`
      sample=`basename $bamfile .bam`
      sID=$sample
      sPU=$sample
      sSM=$sample
      RGparms=$( echo "RGID=${sID}:RGLB=${sLB}:RGPU=${sPU}:RGSM=${sSM}:RGPL=${sPL}:RGCN=${sCN}" )

                qsub_sortnode=$realignlogdir/qsub.sort.$bamfile.$chr
                echo "#PBS -V" > $qsub_sortnode
                echo "#PBS -A $pbsprj" >> $qsub_sortnode
                echo "#PBS -N ${pipeid}_sort_${bamfile}_$chr" >> $qsub_sortnode
      echo "#PBS -l walltime=$pbscpu" >> $qsub_sortnode
      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_sortnode
      echo "#PBS -o $realignlogdir/log.sort.$bamfile.$chr.ou" >> $qsub_sortnode
      echo "#PBS -e $realignlogdir/log.sort.$bamfile.$chr.in" >> $qsub_sortnode
                echo "#PBS -q $pbsqueue" >> $qsub_sortnode
                echo "#PBS -m ae" >> $qsub_sortnode
                echo "#PBS -M $email" >> $qsub_sortnode

      # schedule using a qsub or via the launcher/serial
      if [ $schedule eq "QSUB" ]
      then
                   echo "aprun -n 1 -d $thr $scriptdir/sortnode.sh $picardir $samdir $javamodule $outputdir $bam $bamfile.$chr.bam $bamfile.$chr.sorted.bam $RGparms $flag $chr $realignlogdir/log.sort.$bamfile.$chr.in $realignlogdir/log.sort.$bamfile.$chr.ou $email $realignlogdir/qsub.sort.$bamfile.$chr" >> $qsub_sortnode
                   `chmod a+r $qsub_sortnode`
                   sortjobid=`qsub $qsub_sortnode`
                   # new line to avoid hiccup
                   #`qhold -h u $sortjobid`
                   echo $sortjobid >> $outputrootdir/logs/REALSORTEDpbs.$bamfile
                   echo $sortjobid >> $outputrootdir/logs/REALSORTED.$bamfile_$chr
           else 
         # this block constructs the list of jobs to perform in serial on that chromosome, on that sample
                   # clear out the jobs file for that chromosome, for that sample   
              truncate -s 0 $realrecal.$bamfile.$chr.serialjobs
         # add this jobs file into joblist for the launcher
              echo "$realignlogdir realrecal.$bamfile.$chr.serialjobs" >> $realignlogdir/$realrecal.AnisimovJoblist

                   echo "nohup $scriptdir/sortnode.sh $picardir $samdir $javamodule $outputdir $bam $bamfile.$chr.bam $bamfile.$chr.sorted.bam $RGparms $flag $chr $realignlogdir/log.sort.$bamfile.$chr.in $realignlogdir/log.sort.$bamfile.$chr.ou $email $realignlogdir/qsub.sort.$bamfile.$chr > $realignlogdir/log.sort.$bamfile.$chr.in" >> $realignlogdir/realrecal.$bamfile.$chr.serialjobs
      fi

      chrinfiles[$inx]=${chrinfiles[$inx]}":-I:$outputdir/$bamfile.$chr.sorted.bam"
      chrinputfiles[$inx]=${chrinputfiles[$inx]}":INPUT=$outputdir/$bamfile.$chr.sorted.bam"
       done # end loop over samples
       echo `date`


            sortid=$( cat $outputrootdir/logs/REALSORTED.$bamfile_$chr | sed "s/\..*//g" | tr "\n" ":" )
            outputfile=$chr.realrecal.$bamfile.output.bam
            igv_files=${igv_files}":INPUT=${outputfile}"
       echo "realign-recalibrate for interval:$chr..."
       qsub_realrecalold=$realignlogdir/qsub.realrecal.$bamfile.$chr
       echo "#PBS -V" > $qsub_realrecalold
       echo "#PBS -A $pbsprj" >> $qsub_realrecalold
       echo "#PBS -N ${pipeid}_realrecal_$bamfile.$chr" >> $qsub_realrecalold
       echo "#PBS -l walltime=$pbscpu" >> $qsub_realrecalold
       echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_realrecalold
       echo "#PBS -o $realignlogdir/log.realrecal.$bamfile.$chr.ou" >> $qsub_realrecalold
       echo "#PBS -e $realignlogdir/log.realrecal.$bamfile.$chr.in" >> $qsub_realrecalold
       echo "#PBS -q $pbsqueue" >> $qsub_realrecalold
       echo "#PBS -m ae" >> $qsub_realrecalold
       echo "#PBS -M $email" >> $qsub_realrecalold
       echo "#PBS -W depend=afterok:$sortid" >> $qsub_realrecalold

       if [ $schedule eq "QSUB" ]
       then
          echo "aprun -n 1 -d $thr $scriptdir/realrecalold.sh $outputdir $outputfile $chr ${chrinfiles[$inx]} ${chrinputfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $realignlogdir/log.realrecal.$bamfile.$chr.in $realignlogdir/log.realrecal.$bamfile.$chr.ou $email $realignlogdir/qsub.realrecal.$bamfile.$chr" >> $qsub_realrecalold
          `chmod a+r $qsub_realrecalold`
          recaljobid=`qsub $qsub_realrecalold`
               # new line to avoid hiccup
               #`qhold -h u $recaljobid`
          echo $recaljobid >> $outputrootdir/logs/REALRECALpbs.$bamfile
       else 
          # this block constructs the list of jobs to perform in serial on that chromosome, on that sample
          echo "nohup $scriptdir/realrecalold.sh $outputdir $outputfile $chr ${chrinfiles[$inx]} ${chrinputfiles[$inx]} ${region[$inx]} $realparms $recalparms $runfile $flag $realignlogdir/log.realrecal.$bamfile.$chr.in $realignlogdir/log.realrecal.$bamfile.$chr.ou $email $realignlogdir/qsub.realrecal.$bamfile.$chr > $realignlogdir/log.realrecal.$bamfile.$chr.in" >> $realignlogdir/realrecal.$bamfile.$chr.serialjobs
       fi

            if [ $skipvcall == "NO" ]
            then
                vardir=$( echo $outputdir | sed 's/realign/variant/' )
                varlogdir=$( echo $realignlogdir | sed 's/realign/variant/' )
  
      vcf_files=${vcf_files}":${outputfile}"
                echo "variant calling call.."
      qsub_vcallgatk=$varlogdir/qsub.vcallgatk.$bamfile.$chr
      echo "#PBS -V" > $qsub_vcallgatk
      echo "#PBS -A $pbsprj" >> $qsub_vcallgatk
      echo "#PBS -N ${pipeid}_vcall_$chr" >> $qsub_vcallgatk
      echo "#PBS -l walltime=$pbscpu" >> $qsub_vcallgatk
      echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_vcallgatk
      echo "#PBS -o $varlogdir/log.vcallgatk.$bamfile.$chr.ou" >> $qsub_vcallgatk
      echo "#PBS -e $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $qsub_vcallgatk
      echo "#PBS -q $pbsqueue" >> $qsub_vcallgatk
      echo "#PBS -m ae" >> $qsub_vcallgatk
      echo "#PBS -M $email" >> $qsub_vcallgatk
      echo "#PBS -W depend=afterok:$recaljobid" >> $qsub_vcallgatk
           if [ $schedule eq "QSUB" ]
      then
         echo "aprun -n 1 -d $thr $scriptdir/vcallgatk.sh $vardir $outputdir $outputfile $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$bamfile.$chr.in $varlogdir/log.vcallgatk.$bamfile.$chr.ou $email $varlogdir/qsub.vcallgatk.$bamfile.$chr" >> $qsub_vcallgatk
         `chmod a+r $qsub_vcallgatk`
         vcalljobid=`qsub $qsub_vcallgatk`
         echo $vcalljobid >> $outputrootdir/logs/VCALLGATKpbs.$bamfile
      else
# something is messed up with redirection here
         echo "nohup $scriptdir/vcallgatk.sh $vardir $outputdir $outputfile $chr ${region[$inx]} $runfile $varlogdir/log.vcallgatk.$bamfile.$chr.in $varlogdir/log.vcallgatk.$bamfile.$chr.ou $email $varlogdir/qsub.vcallgatk.$bamfile.$chr >> $qsub_vcallgatk > $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $realignlogdir/realrecal.$bamfile.$chr.serialjobs
           fi
            else
                     echo "variant calling will not be run"
            fi


            (($inx)) ########### make sure this works
       echo `date`
            echo "bottom of the loop over chromosomes"
     done
     echo `date`


     

     if [ $schedule eq "ANISIMOV" ]
     then

################### figure out $numsamples !!!!!!!!!!!!!!!!!!#

        echo "scheduling Anisimov Launcher joblist"
   qsub_aisimov=$varlogdir/qsub.vcallgatk.$bamfile.$chr
   echo "#PBS -V" > $qsub_aisimov
   echo "#PBS -A $pbsprj" >> $qsub_aisimov
   echo "#PBS -N ${pipeid}_vcall_$chr" >> $qsub_aisimov
   echo "#PBS -l walltime=$pbscpu" >> $qsub_aisimov
   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_aisimov
   echo "#PBS -o $varlogdir/log.vcallgatk.$bamfile.$chr.ou" >> $qsub_aisimov
   echo "#PBS -e $varlogdir/log.vcallgatk.$bamfile.$chr.in" >> $qsub_aisimov
   echo "#PBS -q $pbsqueue" >> $qsub_aisimov
   echo "#PBS -m ae" >> $qsub_aisimov
   echo "#PBS -M $email" >> $qsub_aisimov

        echo "aprun -n $numsamples -N 1 -d $thr ~anisimov/scheduler/scheduler.x $realignlogdir/$realrecal.AnisimovJoblist /bin/bash > $realignlogdir/RealRecalVarCallAnisimov.joblist.log" >> $qsub_anisimov
$realignlogdir/$realrecal.AnisimovJoblist
 
        qsub $qsub_anisimov

     fi

     echo "clean up and produce summary table"
     if [ $skipvcall == "NO" ]
     then
    listjobids=$( cat $outputrootdir/logs/VCALLGATKpbs.$bamfile | sed "s/\..*//g" | tr "\n" ":" )
     else
    listjobids=$( cat $outputrootdir/logs/REALRECALpbs.$bamfile | sed "s/\..*//g" | tr "\n" ":" )
     fi

     # preparing output files for realignment and/or variant calling

     qsub_igvbam=$outputrootdir/logs/qsub.igvbam.$bamfile
     echo "#PBS -V" > $qsub_igvbam
     echo "#PBS -A $pbsprj" >> $qsub_igvbam
     echo "#PBS -N ${pipeid}_igvbam.$bamfile" >> $qsub_igvbam
     echo "#PBS -l walltime=$pbscpu" >> $qsub_igvbam
     echo "#PBS -l nodes=1:ppn=$thr" >> $qsub_igvbam
     echo "#PBS -o $outputrootdir/logs/log.igvbam.$bamfile.ou" >> $qsub_igvbam
     echo "#PBS -e $outputrootdir/logs/log.igvbam.$bamfile.in" >> $qsub_igvbam
     echo "#PBS -q $pbsqueue" >> $qsub_igvbam
     echo "#PBS -m ae" >> $qsub_igvbam
     echo "#PBS -M $email" >> $qsub_igvbam
     echo "#PBS -W depend=afterok:$listjobids" >> $qsub_igvbam
     echo "aprun -n 1 -d $thr $scriptdir/igvbam.sh $outputdir $igv_files $runfile $outputrootdir/logs/log.igvbam.$bamfile.in $outputroot/logs/log.igvbam.$bamfile.ou $email $outputrootdir/logs/qsub.igvbam.$bamfile"  >> $qsub_igvbam
     `chmod a+r $qsub_igvbam`
     igvjobid=`qsub $qsub_igvbam`
     echo $igvjobid >> $outputrootdir/logs/IGVBAMpbs.$bamfile

     if [ $cleanupflag == "YES" ]
     then
    qsub_cleanup=$outputrootdir/logs/qsub.cleanup.realn.$bamfile
    echo "#PBS -V" > $qsub_cleanup
    echo "#PBS -A $pbsprj" >> $qsub_cleanup
    echo "#PBS -N ${pipeid}_cleanup_realn.$bamfile" >> $qsub_cleanup
    echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
    echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
    echo "#PBS -o $outputrootdir/logs/log.cleanup.realn.$bamfile.ou" >> $qsub_cleanup
    echo "#PBS -e $outputrootdir/logs/log.cleanup.realn.$bamfile.in" >> $qsub_cleanup
    echo "#PBS -q $pbsqueue" >> $qsub_cleanup
    echo "#PBS -m ae" >> $qsub_cleanup
    echo "#PBS -M $email" >> $qsub_cleanup
    echo "#PBS -W depend=afterok:$igvjobid" >> $qsub_cleanup
    echo "$scriptdir/cleanup.sh $outputrootdir $analysis $outputrootdir/logs/log.cleanup.realn.$bamfile.in $outputrootdir/logs/log.cleanup.realn.$bamfile.ou $email $outputrootdir/logs/qsub.cleanup.realn.$bamfile"  >> $qsub_cleanup
    `chmod a+r $qsub_cleanup`
    cleanjobid=`qsub $qsub_cleanup`
    echo $cleanjobid >> $outputrootdir/logs/CLEANUPpbs.$bamfile
     fi

     if [ $skipvcall == "NO" ]
     then
    qsub_mergevcf=$varlogdir/qsub.mergevcf.$bamfile
    echo "#PBS -V" > $qsub_mergevcf
    echo "#PBS -A $pbsprj" >> $qsub_mergevcf
    echo "#PBS -N ${pipeid}_mergevcf.$bamfile" >> $qsub_mergevcf
    echo "#PBS -l walltime=$pbscpu" >> $qsub_mergevcf
    echo "#PBS -l nodes=1:ppn=1" >> $qsub_mergevcf
    echo "#PBS -o $varlogdir/log.mergevcf.$bamfile.ou" >> $qsub_mergevcf
    echo "#PBS -e $varlogdir/log.mergevcf.$bamfile.in" >> $qsub_mergevcf
    echo "#PBS -q $pbsqueue" >> $qsub_mergevcf
    echo "#PBS -m ae" >> $qsub_mergevcf
    echo "#PBS -M $email" >> $qsub_mergevcf
    echo "#PBS -W depend=afterok:$listjobids" >> $qsub_mergevcf
    echo "$scriptdir/mergevcf.sh $vardir $vcf_files $tabixdir $vcftoolsdir $varlogdir/log.mergevcf.$bamfile.in $varlogdir/log.mergevcf.$bamfile.ou $email $varlogdir/qsub.mergevcf.$bamfile"  >> $qsub_mergevcf
    `chmod a+r $qsub_mergevcf`
    mergevcfjobid=`qsub $qsub_mergevcf`
    echo $mergevcfjobid >> $outputrootdir/logs/MERGEVCFpbs

    if [ $cleanupflag == "YES" ]
    then
        qsub_cleanup=$outputrootdir/logs/qsub.cleanup.vcall.$bamfile
        echo "#PBS -V" > $qsub_cleanup
        echo "#PBS -A $pbsprj" >> $qsub_cleanup
        echo "#PBS -N ${pipeid}_cleanup_vcall.$bamfile" >> $qsub_cleanup
        echo "#PBS -l walltime=$pbscpu" >> $qsub_cleanup
        echo "#PBS -l nodes=1:ppn=1" >> $qsub_cleanup
        echo "#PBS -o $outputrootdir/logs/log.cleanup.vcall.$bamfile.ou" >> $qsub_cleanup
        echo "#PBS -e $outputrootdir/logs/log.cleanup.vcall.$bamfile.in" >> $qsub_cleanup
        echo "#PBS -q $pbsqueue" >> $qsub_cleanup
        echo "#PBS -m ae" >> $qsub_cleanup
        echo "#PBS -M $email" >> $qsub_cleanup
        echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub_cleanup
        echo "$scriptdir/cleanup.sh $outputrootdir VCALL_ONLY $outputrootdir/logs/log.cleanup.vcall.$bamfile.in $outputrootdir/logs/log.cleanup.vcall.$bamfile.ou $email $outputrootdir/logs/qsub.cleanup.vcall.$bamfile"  >> $qsub_cleanup
        `chmod a+r $qsub_cleanup`
        cleanjobid=`qsub $qsub_cleanup`
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

     qsub_summary=$outputrootdir/logs/qsub.summary.realn.allok.$bamfile
     echo "#PBS -V" > $qsub_summary
     echo "#PBS -A $pbsprj" >> $qsub_summary
     echo "#PBS -N ${pipeid}_summaryok.$bamfile" >> $qsub_summary
     echo "#PBS -l walltime=$pbscpu" >> $qsub_summary
     echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
     echo "#PBS -o $outputrootdir/logs/log.summary.realn.$bamfile.ou" >> $qsub_summary
     echo "#PBS -e $outputrootdir/logs/log.summary.realn.$bamfile.in" >> $qsub_summary
     echo "#PBS -q $pbsqueue" >> $qsub_summary
     echo "#PBS -m ae" >> $qsub_summary
     echo "#PBS -M $email" >> $qsub_summary
     if [ `expr ${#cleanupjobs}` -gt 0 ]
     then 
    echo "#PBS -W depend=afterok:$cleanupjobs" >> $qsub_summary
     else
         if [ $skipvcall == "YES" ]
         then
        echo "#PBS -W depend=afterok:$igvjobid" >> $qsub_summary
         else
        echo "#PBS -W depend=afterok:$mergevcfjobid" >> $qsub_summary
         fi
     fi
     echo "$scriptdir/summary.sh $outputrootdir $email exitok"  >> $qsub_summary
     `chmod a+r $qsub_summary`
     lastjobid=`qsub $qsub_summary`
     echo $lastjobid >> $outputrootdir/logs/SUMMARYpbs.$bamfile

     if [ `expr ${#lastjobid}` -lt 1 ]
     then
         echo "at least one job aborted"
    qsub_summary=$outputrootdir/logs/qsub.summary.realn.afterany.$bamfile
    echo "#PBS -V" > $qsub_summary
    echo "#PBS -A $pbsprj" >> $qsub_summary
    echo "#PBS -N ${pipeid}_summary_afterany.$bamfile" >> $qsub_summary
    echo "#PBS -l walltime=$pbscpu" >> $qsub_summary
    echo "#PBS -l nodes=1:ppn=1" >> $qsub_summary
    echo "#PBS -o $outputrootdir/logs/log.summary.realn.afterany.$bamfile.ou" >> $qsub_summary
    echo "#PBS -e $outputrootdir/logs/log.summary.realn.afterany.$bamfile.in" >> $qsub_summary
    echo "#PBS -q $pbsqueue" >> $qsub_summary
    echo "#PBS -m ae" >> $qsub_summary
    echo "#PBS -M $email" >> $qsub_summary
    echo "#PBS -W depend=afterany:$listjobids" >> $qsub_summary
    echo "$scriptdir/summary.sh $outputtopdir $email exitnotok"  >> $qsub_summary
    `chmod a+r $qsub_summary`
    badjobid=`qsub $qsub_summary`
    echo $badjobid >> $outputrootdir/logs/SUMMARYpbs.$bamfile
     fi
     `chmod -R 770 $outputroordir/logs`
     echo `date`
fi
