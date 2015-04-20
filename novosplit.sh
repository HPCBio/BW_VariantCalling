#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=lmainzer@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# -gt 15 ]
then
	MSG="parameter mismatch."
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        alignerdir=$1
        params=$2
        ref=$3
	outputdir=$4
        samfile=$5
        bamfile=$6
        scriptdir=$7
        runfile=$8
        paired=$9
        R1=${10}

        parameters=$( echo $params | tr "_" " " )
        samprocessing=$( cat $runfile | grep -w SAMPROCESSING | cut -d "=" -f2 )
        sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d "=" -f2 )
        memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
        threads=$( cat $runfile | grep -w PBSTHREADS | cut -d "=" -f2 )


        cd $outputdir
        #######################     PAIRED ALIGNMENT      ####################
        if [ $paired -eq 1 ]
        then
           R2=${11}
           elog=${12}
           olog=${13}
           email=${14}
           qsubfile=${15}
           LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

           if [ $samprocessing == "SAMTOOLS" ]
           then
              echo `date`

              $memprof $alignerdir/novoalign -d $ref -f $R1 $R2 -o SAM $parameters | $samdir/samtools view -bS -o $bamfile 
              novoalign_exitcode=${PIPESTATUS[0]}
              samtools_exitcode=${PIPESTATUS[1]}
              echo `date`

              if [ $novoalign_exitcode -ne 0 ]
              then
                 MSG="novoalign command failed.  exitcode=$novoalign_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $novoalign_exitcode;
              fi
              if [ $samtools_exitcode -ne 0 ]
              then
                 MSG="sam to bam conversion failed.  exitcode=$samtools_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $samtools_exitcode;
              fi
              if [ ! -s $bamfile ]
              then
                  MSG="$outputdir/$bamfile BAM file not created. sam2bam step failed during  novosplit alignment."
                  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  exit 1;
              fi
              echo `date`




           elif [ $samprocessing == "SAMBAMBA" ]
           then
              echo `date`

              $memprof $alignerdir/novoalign -d $ref -f $R1 $R2 -o SAM $parameters | ${sambambadir}/sambamba view -f bam -h --sam-input /dev/stdin -t $threads --output-filename $bamfile
              novoalign_exitcode=${PIPESTATUS[0]}
              sambamba_exitcode=${PIPESTATUS[1]}
              echo `date`

              if [ $novoalign_exitcode -ne 0 ]
              then
                 MSG="novoalign command failed.  exitcode=$novoalign_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $novoalign_exitcode;
              fi
              if [ $sambamba_exitcode -ne 0 ]
              then
                 MSG="sam to bam conversion failed.  exitcode=$sambamba_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $samtools_exitcode;
              fi
              if [ ! -s $bamfile ]
              then
                  MSG="$outputdir/$bamfile BAM file not created. sam2bam step failed during  novosplit alignment."
                  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  exit 1;
              fi
              echo `date`
           fi





        #######################     SINGLE-ENDED ALIGNMENT      ####################
        else
           elog=${11}
           olog=${12}
           email=${13}
           qsubfile=${14}
           LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

           if [ $samprocessing == "SAMTOOLS" ]
           then
              echo `date`

              $memprof $alignerdir/novoalign -d $ref -f $R1 -o SAM $parameters | $samdir/samtools view -bS -o $bamfile
              novoalign_exitcode=${PIPESTATUS[0]}
              samtools_exitcode=${PIPESTATUS[1]}
              echo `date`

              if [ $novoalign_exitcode -ne 0 ]
              then
                 MSG="novoalign command failed.  exitcode=$novoalign_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $novoalign_exitcode;
              fi
              if [ $samtools_exitcode -ne 0 ]
              then
                 MSG="sam to bam conversion failed.  exitcode=$samtools_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $samtools_exitcode;
              fi
              if [ ! -s $bamfile ]
              then
                  MSG="$outputdir/$bamfile BAM file not created. sam2bam step failed during  novosplit alignment."
                  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  exit 1;
              fi
              echo `date`



           elif [ $samprocessing == "SAMBAMBA" ]
           then
              echo `date`

              $memprof $alignerdir/novoalign -d $ref -f $R1 -o SAM $parameters | ${sambambadir}/sambamba view -f bam -h --sam-input /dev/stdin -t $threads --output-filename $bamfile
              novoalign_exitcode=${PIPESTATUS[0]}
              sambamba_exitcode=${PIPESTATUS[1]}
              echo `date`

              if [ $novoalign_exitcode -ne 0 ]
              then
                 MSG="novoalign command failed.  exitcode=$novoalign_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $novoalign_exitcode;
              fi
              if [ $sambamba_exitcode -ne 0 ]
              then
                 MSG="sam to bam conversion failed.  exitcode=$sambamba_exitcode  novosplit alignment failed"
                 #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                 exit $samtools_exitcode;
              fi
              if [ ! -s $bamfile ]
              then
                  MSG="$outputdir/$bamfile BAM file not created. sam2bam step failed during  novosplit alignment."
                  #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                  exit 1;
              fi
              echo `date`
           fi



        fi # end paired-ended/single-ended condition 


fi
