#!/bin/bash
######################################
#  script to convert bam files back to fastq as pre requisite to alignment
#
######################################
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 8 ];
then
	MSG= "parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        umask 0027
        scriptfile=$0
        inputdir=$1
        newfqfiles=$2
        runfile=$3
        samplefileinfo=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        set +x; echo -e "\n#########   check params \n" >&2; set -x

        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        provenance=$( cat $runfile | grep -w PROVENANCE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 )
        if [ ! -d $inputdir ]
        then
	    MSG="INPUTDIR=$inputdir directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
	    MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ $provenance != "SINGLE_SOURCE" -a $provenance != "MULTI_SOURCE" ]
        then
	    MSG="Invalid value for PROVENANCE=$provenance"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        set +x; echo -e "\n#########   checking inputdir for fastq files tha need merging by sample \n" >&2; set -x
 

        cd $inputdir

        if [ `find -name "*.fastq" | wc -l` -lt 1 ]
        then
	    MSG="no fastq files created after bam2fastq conversion"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        newsamples=""
        sep=" "
        newnames="\n"
        bam_samples=""
        allbams2fq=$inputdir/fqfiles4merging
        touch $allbams2fq
        fqfiles=$( echo $newfqfiles | tr ":" " " )

        if [ $provenance == "SINGLE_SOURCE" ]
        then

           set +x; echo -e "\n#########   CASE1: provenance is single source \n" >&2; set -x  
        
        
           echo "fastq files in $inputdir need further processing..."
           for fqfile in $fqfiles
           do
               thesamples=`find -name "*_1.${fqfile}*.fastq" | sed 's/.\///g' | sed 's/_1.[0-9]*.[a-z]*.fastq//'`
               echo -e $thesamples >> $allbams2fq
	  done
	  
          set +x; echo -e "\n#########   constructing list of files that need conversion \n" >&2; set -x
	  
          bam_samples=$( cat $allbams2fq | tr " " "\n" | sort | uniq -c | sed 's/^ *//g'  | tr " " ":")
          for line in $bam_samples
          do
               nline=$( echo $line | sed 's/:://g' | sed 's/^://' )
	       frq=$( echo $nline | cut -d ':' -f1 )
	       b2fsample=$( echo $nline | cut -d ':' -f2 | sed 's/_1.[0-9]*.[a-z]*.fastq//' )
               if [ $frq == "1" ]
               then
               
                  set +x; echo -e "\n#########   this file is unique. collect R1 and R2 \n" >&2; set -x
               
                  readone=${b2fsample}_r1.fastq
                  touch $readone
                  cat ${b2fsample}_1.*.fastq >> $readone
                  exitcode=$?
                  if [ -s $readone -a $exitcode -eq 0 ]
                  then
                      readtwo=${b2fsample}_r2.fastq
                      touch $readtwo
                      cat ${b2fsample}_2.*.fastq >> $readtwo
                      exitcode=$?
                      if [ $exitcode -eq 0 ]
                      then
                          if [ -s $readtwo ]
                          then
			      newnames="FASTQ:${b2fsample}=$inputdir/$readone $inputdir/${readtwo}\n$newnames"
			  else
			      newnames="FASTQ:${b2fsample}=$inputdir/${readone}\n$newnames"
                          fi
		      fi
                      newsamples=${b2fsample}${sep}$newsamples
                  else
                      MSG="$b2fsample Empty fastq file.  exitcode=$exitcode updateconfi failed"
		      echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		      exit 1;
                  fi
               else
                  if [ $frq != "" ]
                  then

                      set +x; echo -e "\n#########   multisamples with same name. concatenate them... \n" >&2; set -x
               
                      readone=${b2fsample}_r1.fastq
                      touch $readone
                      cat ${b2fsample}_1.*.fastq >> $readone
                      if [ -s $readone ]
                      then
                          readtwo=${b2fsample}_r2.fastq
                          touch $readtwo
			  cat ${b2fsample}_2.*.fastq >> $readtwo
                          if [ -s $readtwo ]
                          then
			      newnames="FASTQ:${b2fsample}=$inputdir/$readone $inputdir/${readtwo}\n$newnames"
			  else
			      newnames="FASTQ:${b2fsample}=$inputdir/${readone}\n$newnames"
			  fi
			  newsamples=${b2fsample}${sep}$newsamples
                          #`rm ${b2fsample}_1.*fastq`
                          #`rm ${b2fsample}_2.*fastq`
                      else
			  MSG="$readone Empty fastq file. updateconfig failed"
			  echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
			  exit 1;
                      fi
                  fi
               fi           
           done
	fi
        if [ $provenance == "MULTI_SOURCE" ]
        then

               set +x; echo -e "\n#########   CASE2: provenance is multi-source. fastq files in input dir do not need further processing \n" >&2; set -x  
        
               for fqfile in $fqfiles 
               do
		   if [ `find -name "*${fqfile}*" | wc -l` -gt 1 ]
                   then
                       # paired-end reads
		       R1=${fqfile}_R1.fastq
                       R2=${fqfile}_R2.fastq
		       newnames="FASTQ:${fqfile}=$inputdir/$R1 $inputdir/${R2}\n$newnames"
                       newsamples=${fqfile}${sep}${newsamples}
                   else
                       # single reads
		       R1=${fqfile}_R1.fastq
		       newnames="FASTQ:${fqfile}=$inputdir/$R1\n$newnames"
                       newsamples=${fqfile}${sep}${newsamples}
                   fi
               done
        
        fi

        set +x; echo -e "\n#########   Udpating configuration files \n" >&2; set -x  

        if [ `expr ${#newnames}` -gt 0 ]
        then
            directory=`dirname $samplefileinfo`
            oldfile=`basename $samplefileinfo`
            oldfilename=$directory/$oldfile.old
            newfilename=$directory/$oldfile
            mv $samplefileinfo $oldfilename
            echo -e $newnames >> $newfilename
            if [ ! -s $newfilename ]
            then
		MSG="Empty file samplefilenames after bam2fastq conversion. updateconfig failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        else
	    MSG="Empty file samplefilenames after bam2fastq conversion. updateconfig failed"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        
        set +x; echo -e "\n#########   Udpating run files \n" >&2; set -x  
        
        directory=`dirname $runfile`
        oldfile=`basename $runfile`
        oldrun=$directory/$oldfile.old
        newrun=$directory/$oldfile
        mv $runfile $oldrun
        `perl $scriptdir/lned.pl $oldrun $newrun SAMPLENAMES "${newsamples}"`
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
 	    MSG="lned.pl failed.  exitcode=$exitcode. update config files after bam2fastq conversion failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        else
            if [ ! -s $newrun ]
            then
		MSG="Empty configuration files after bam2fastq conversion"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        fi
	echo `date`
fi
