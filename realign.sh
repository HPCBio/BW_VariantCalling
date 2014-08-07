#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#
# realign.sh
# Second module in the GGPS analysis pipeline
redmine=hpcbio-redmine@igb.illinois.edu
#redmine=grendon@illinois.edu
if [ $# != 5 ]
then
	MSG="parameter mismatch. "
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        runfile=$1
        elog=$2
        olog=$3
        email=$4
        qsubfile=$5
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        if [ !  -s $runfile ]
        then
           MSG="$runfile configuration file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi


echo -e "\n\n\n##################################### PARSING RUN INFO FILE ########################################\n\n\n"


        sampledir=$( cat $runfile | grep -w SAMPLEDIR | cut -d '=' -f2 )
	outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        thr=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        type=$( cat $runfile | grep -w TYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        igvdir=$( cat $runfile | grep -w IGVDIR | cut -d '=' -f2 )
        extradir=$outputdir/extractreads
        realrecalflag=$( cat $runfile | grep -w REALIGNORDER | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
        paired=$( cat $runfile | grep -w PAIRED | cut -d '=' -f2 )
        rlen=$( cat $runfile | grep -w READLENGTH | cut -d '=' -f2 )
        multisample=$( cat $runfile | grep -w MULTISAMPLE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        samples=$( cat $runfile | grep -w SAMPLENAMES | cut -d '=' -f2 | tr ":" "\n")
        region=$( cat $runfile | grep -w CHRINDEX | cut -d '=' -f2 )
        resortflag=$( cat $runfile | grep -w RESORTBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
        revertsam=$( cat $runfile | grep -w REVERTSAM | cut -d '=' -f2  )
        indices=$( echo $region | sed 's/^/chr/' | sed 's/:/ chr/g' )
        epilogue=$( cat $runfile | grep -w EPILOGUE | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )



####################################################################################################
#####################################                       ########################################
#####################################  CREATE  DIRECTORIES  ########################################
#####################################                       ########################################
####################################################################################################

echo -e "\n\n\n#####################################  CREATE  DIRECTORIES  ########################################\n\n\n"



        RealignOutputDir=$outputdir/realign
        if [ -d $RealignOutputDir ]
        then
           echo "$RealignOutputDir already exists"
        else
           mkdir -p $RealignOutputDir
        fi

        TopOutputLogs=$outputdir/logs
        if [ -d $TopOutputLogs ]
        then
           echo "$TopOutputLogs already exists"
           pbsids=""
        else
           mkdir -p $TopOutputLogs
        fi

        RealignOutputLogs=$TopOutputLogs/realign
        if [ ! -d $RealignOutputLogs ]
        then
            mkdir $RealignOutputLogs
        fi
        `chmod -R 770 $RealignOutputLogs/`
        # where messages about failures will go
        truncate -s 0 $RealignOutputLogs/FAILEDmessages
        if [ ! -d $RealignOutputLogs/FAILEDjobs ]
        then
            mkdir $RealignOutputLogs/FAILEDjobs
        else
            rm -r $RealignOutputLogs/FAILEDjobs/*
        fi
        `chmod -R 770 $RealignOutputLogs/FAILEDjobs`



        if [ $type == "GENOME" -o $type == "WHOLE_GENOME" -o $type == "WHOLEGENOME" -o $type == "WGS" ]
        then
            pbscpu=$( cat $runfile | grep -w PBSCPUALIGNWGEN | cut -d '=' -f2 )
            pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN | cut -d '=' -f2 )
        else
            if [ $type == "EXOME" -o $type == "WHOLE_EXOME" -o $type == "WHOLEEXOME" -o $type == "WES" ]
            then
		pbscpu=$( cat $runfile | grep -w PBSCPUALIGNEXOME | cut -d '=' -f2 )
		pbsqueue=$( cat $runfile | grep -w PBSQUEUEEXOME | cut -d '=' -f2 )
            else
		MSG="Invalid value for TYPE=$type in configuration file."
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
                exit 1;
            fi
        fi

        if [ $resortflag != "1" -a $resortflag != "0" -a $resortflag != "YES" -a $resortflag != "NO" ]
        then
           MSG="Invalid value for RESORTBAM=$resortflag"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email""
            exit 1;
        else
            if [ $resortflag == "1" ]
            then
                $resortflag="YES"
            fi
            if [ $resortflag == "0" ]
            then
                $resortflag="NO"
            fi
        fi

        if [ -z $epilogue ]
        then
           MSG="Invalid value for EPILOGUE=$epilogue in configuration file"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        else
           `chmod 750 $epilogue`
        fi
        if [ ! -d $scriptdir ]
        then
           MSG="$scriptdir script directory not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $refdir ]
        then
           MSG="$refdir directory of reference genome was not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -d $picardir ]
        then
           MSG="$picardir picard directory was not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
           MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi


        # check that sampledir exists
        if [ ! -d $sampledir ]
        then
           MSG="SAMPLEDIR=$sampledir file not found"
           echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
           #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
           exit 1;
        fi

        # construct a list of SampleNames, check that files actually exist
        if [ ! -e $TopOutputLogs/SAMPLENAMES.list ]
        then
           numsamples=0
           truncate -s 0 $TopOutputLogs/SAMPLENAMES.tmp.list
           for fastqfile in $sampledir/*
           do
               # strip path, which read (left/right), and extension from input files
               # and put that info into the SampleNames file
               SampleName=$( basename $fastqfile | sed 's/_read.\?\.[^.]*$//' )
               echo -e "$SampleName" >> $TopOutputLogs/SAMPLENAMES.tmp.list
               let numsamples+=1
           done
           # paired-ended fastq will produce duplicate lines in the SampleNames file, so remove the duplicates
           uniq  $TopOutputLogs/SAMPLENAMES.tmp.list >  $TopOutputLogs/SAMPLENAMES.list
           rm  $TopOutputLogs/SAMPLENAMES.tmp.list
        else
           numsamples=`wc -l $TopOutputLogs/SAMPLENAMES.list | cut -d ' ' -f 1`
        fi

        if [ $numsamples -lt 1 ]
            then
              MSG="No samples found in SAMPLEDIR=$sampledir."
              echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
              #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
              exit 1;
            fi
        if [ $numsamples -gt 1 -a $multisample == "YES" ]
        then
            echo "multiple samples to be aligned."
        else
           if [ $numsamples -eq 1 -a $multisample == "NO" ]
           then
              echo "single sample to be aligned."
           fi
        fi


        pipeid=$( cat $output_logs/MAINpbs )
        
        if [ ! -d $outputdir ]
        then
	    mkdir -p $outputdir
        fi
        if [ ! -d $output_logs ]
        then
	    mkdir  -p $output_logs
	fi  

		if [ $realrecalflag != "1" -a $realrecalflag != "0" ]
		then
		    echo "realign-recalibration order flag is not set properly. Default value [1] will be assiged to it"
		    realrecalflag="1"
		fi
		

		listfiles="";
		sep=":";
		JOBSmayo=""
		JOBSncsa=""

		# finding all aligned BAMs to be realigned-recalibrated

		if [ $analysis == "REALIGN" -o $analysis == "REALIGNMENT" ]
		then
		    echo "alignment was done inhouse. no need to resort"
		    echo "We need to wait until the alignment jobs enter the queue"

		    while [ ! -s $output_logs/ALIGN_NCSA_jobids ]
		    do
			`sleep 60s`
		    done

		    JOBSncsa=$( cat $output_logs/ALIGN_NCSA_jobids | sed "s/\..*//g" | tr "\n" ":" | sed "s/::/:/g" )
		else
		    echo "we need to check entries in samplefileinfo before launching other realignment analyses"
		    while read sampledetail
		    do
			echo "processing next line in file ..."
			if [ `expr ${#sampledetail}` -lt 7 ]
			then
			    echo "skip empty line"
			else
			    echo "preprocessing for realignment $sampledetail"
			    bamfile=$( echo $sampledetail | cut -d '=' -f2 )
			    sampletag=$( echo $sampledetail | cut -d '=' -f1 | cut -d ':' -f2 )
			    if [ `expr ${#bamfile}` -lt 1 ]
			    then
				MSG="parsing SAMPLEFILENAMES file failed. realignment failed."
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
				exit 1;
			    fi
			    if [ `expr ${#sampletag}` -lt 1 ]
			    then
				MSG="parsing SAMPLEFILENAMES file failed. realignment failed."
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
				exit 1;
			    fi
			
			    if [ ! -s $bamfile ]
			    then
				MSG="parsing SAMPLEFILENAMES file failed. realignment failed"
				echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
				exit 1;
			    fi
			fi
                    done <  $TopOutputLogs/SAMPLENAMES.list
                    # end loop over input samples

		fi   

		if [ $resortflag == "YES" -a $analysis == "REALIGN_ONLY" ]
		then
		    echo "alignment was NOT done inhouse. Need to resort bam files. Checking input files"
                    while read SampleName
                    do
                       echo "processing next sample"
                          if [ `expr ${#SampleName}` -lt 7 ]
                          then
                             echo "skipping empty line"
                          else

                          echo "realigning: $SampleName"

			    prefix=`basename $SampleName .wrg.sorted.bam`
			    outputalign=$outputdir/align/$prefix
			    outputlogs=$output_logs/align/$prefix
			    tmpbamfile=$SampleName
			    sortedplain=${prefix}.wrg.sorted.bam
			    sorted=${prefix}.wdups.sorted.bam
			    sortednodups=${prefix}.nodups.sorted.bam

			    if [ ! -d $outputalign ]
			    then
				mkdir -p $outputalign
                                if [ ! -d $outputlogs ]
			        then
			            mkdir -p $outputlogs
                                else
                                    `rm -r $outputlogs/*`
			        fi
		            fi

		   qsub1=$outputlogs/qsub.sortbammayo.$prefix
		   echo "#PBS -V" > $qsub1
		   echo "#PBS -A $pbsprj" >> $qsub1
		   echo "#PBS -N ${pipeid}_sortbamayo_${prefix}" >> $qsub1
                   echo "#PBS -l epilogue=$epilogue" >> $qsub1
		   echo "#PBS -l walltime=$pbscpu" >> $qsub1
		   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub1
		   echo "#PBS -o $outputlogs/log.sortbammayo.${prefix}.ou" >> $qsub1
		   echo "#PBS -e $outputlogs/log.sortbammayo.${prefix}.in" >> $qsub1
		   echo "#PBS -q $pbsqueue" >> $qsub1
		   echo "#PBS -m ae" >> $qsub1
		   echo "#PBS -M $email" >> $qsub1
		   echo "aprun -n 1 -d $thr $scriptdir/sortbammayo.sh $outputalign $tmpbamfile $sortedplain $sorted $sortednodups $runfile $outputlogs/log.sortbammayo.${prefix}.in $outputlogs/log.sortbammayo.${prefix}.ou $email $outputlogs/qsub.sortbammayo.${prefix}" >> $qsub1
		   `chmod a+r $qsub1`
                   sortid=`qsub $qsub1`
                   #`qhold -h u $sortid`
               	   echo $sortid >> $output_logs/REALSORTEDMAYOpbs
		   echo "extracting reads"
		   qsub2=$outputlogs/qsub.extractreadsbam.$prefix
		   echo "#PBS -V" > $qsub2
		   echo "#PBS -A $pbsprj" >> $qsub2
		   echo "#PBS -N ${pipeid}_extrbam${prefix}" >> $qsub2
                   echo "#PBS -l epilogue=$epilogue" >> $qsub2
		   echo "#PBS -l walltime=$pbscpu" >> $qsub2
		   echo "#PBS -l nodes=1:ppn=$thr" >> $qsub2
		   echo "#PBS -o $outputlogs/log.extractreadsbam.${prefix}.ou" >> $qsub2
		   echo "#PBS -e $outputlogs/log.extractreadsbam.${prefix}.in" >> $qsub2
		   echo "#PBS -q $pbsqueue" >> $qsub2
		   echo "#PBS -m ae" >> $qsub2
		   echo "#PBS -M $email" >> $qsub2
                   echo "#PBS -W depend=afterok:$sortid" >> $qsub2
                   echo "aprun -n 1 -d $thr $scriptdir/extract_reads_bam.sh $outputalign $sorted $runfile $outputlogs/log.extractreadsbam.${prefix}.in $outputlogs/log.extractreadsbam.${prefix}.ou $outputlogs/qsub.extractreadsbam.$prefix $igvdir $extradir" >> $qsub2
		   `chmod a+r $qsub2`               
                   extraid=`qsub $qsub2`
                   #`qhold -h u $extraid`
                   echo $extraid >> $output_logs/EXTRACTREADSpbs
                   listfiles=$outputalign/$sorted${sep}${listfiles}
		fi
            done <  $TopOutputLogs/SAMPLENAMES.list
            # end loop over input samples

            cp $output_logs/REALSORTEDMAYOpbs $output_logs/ALN_MAYO_jobids
            JOBSmayo=$( cat $output_logs/ALN_MAYO_jobids | sed "s/\..*//g" | tr "\n" ":" | sed "s/::/:/g" )
        fi

        if [ $resortflag == "NO" -a $analysis == "REALIGN_ONLY" ]
        then
	    echo "alignment was NOT done inhouse. BAM files will not be resorted"
            if [ $revertsam == "0" -o $revertsam == "NO" ]
            then
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!
#################################THIS IS WRONG AND HAS TO BE FIXED!!!!!!!!!!!!!!!!!!!!!!!!!!!

                echo "input is aligned bam that is suitable for realignment and recalibration... no need for preprocessing"
                while read sampledetail
                do
		    bam=$( echo $sampledetail | cut -d "=" -f2 )
		    listfiles=${bam}${sep}${listfiles}
                done <  $TopOutputLogs/SAMPLENAMES.list
                # end loop over input samples

            else
                MSG="invalid value for preprocessing this kind of input: aligned bam. set RESORTBAM=YES and rerun the pipeline"
	        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
		exit 1;
           fi
        fi


        # grab job ids for align and for preprocessing done in this module

        alignids=$( cat $output_logs/ALIGNEDpbs | sed "s/\..*//" | tr "\n" " " )
        mergeids=$( cat $output_logs/MERGEDpbs | sed "s/\..*//" | tr "\n" " " )
        sortedmayoids=$( cat $output_logs/REALSORTEDMAYOpbs | sed "s/\..*//" | tr "\n" " " )
        extraids=$( cat $output_logs/EXTRACTREADSpbs | sed "s/\..*//" | tr "\n" " " )



        if [ $multisample == "NO" ]
        then

           while read SampleName
           do
              echo "processing next sample"
              if [ `expr ${#SampleName}` -lt 7 ]
              then
                 echo "skipping empty line"
              else

              echo "realigning: $SampleName"


                 listfiles=$outputdir/align/${SampleName}

                 realigndir=$outputdir/realign/$SampleName
                 if [ ! -d $realigndir ]
                 then
                     mkdir -p $realigndir
                 else
                     echo "$realigndir already exists. resetting it"
                     `rm -r $realigndir/*`
                 fi

                 if [ $skipvcall == "NO" ]
                 then
                     vartopdir=$outputdir/variant
                     varlogdir=$outputdir/logs/variant
                     if [ ! -d $vartopdir ]
                     then
                         mkdir -p $vartopdir
                     fi
                     if [ ! -d $varlogdir ]
                     then
                         mkdir -p $varlogdir
                     else
                         `rm $varlogdir/*`
                     fi
                     vardir=$outputdir/variant/$SampleName
                     if [ ! -d $vardir ]
                     then
                         mkdir -p $vardir
                     fi

                 fi




 	         qsub3=$RealignOutputLogs/qsub.${SampleName}.realign_new
                 echo "#PBS -V" > $qsub3
         	 echo "#PBS -A $pbsprj" >> $qsub3
	         echo "#PBS -N ${pipeid}_realign_new.${SampleName}" >> $qsub3
                 echo "#PBS -l epilogue=$epilogue" >> $qsub3
	         echo "#PBS -l walltime=$pbscpu" >> $qsub3
         	 echo "#PBS -l nodes=1:ppn=1" >> $qsub3
	         echo "#PBS -o $RealignOutputLogs/log.${SampleName}.realign_new.ou" >> $qsub3
         	 echo "#PBS -e $RealignOutputLogs/log.${SampleName}.realign_new.in" >> $qsub3
	         echo "#PBS -q $pbsqueue" >> $qsub3
         	 echo "#PBS -m ae" >> $qsub3
	         echo "#PBS -M $email" >> $qsub3
                 if [ `expr ${#JOBSmayo}` -gt 0 ]
                 then
                     echo "#PBS -W depend=afterok:$JOBSmayo" >> $qsub3
                 else
                     echo "#PBS -W depend=afterok:$JOBSncsa" >> $qsub3
                 fi
                 echo "$scriptdir/realign_new.sh $realigndir $RealignOutputLogs $listfiles $runfile $realrecalflag $RealignOutputLogs/log.${SampleName}.realign_new.in $RealignOutputLogs/log.${SampleName}.realign_new.ou $email $RealignOutputLogs/qsub.${SampleName}.realign_new" >> $qsub3
	         `chmod a+r $qsub3`               
                 realrecaljob=`qsub $qsub3`
                 #`qhold -h u $realrecaljob` 
                 echo $realrecaljob >> $output_logs/RECALLpbs
              fi
           done <  $TopOutputLogs/SAMPLENAMES.list
           # end loop over input samples

        else


           realigndir=$outputdir/realign
           if [ ! -d $realigndir ]
           then
               mkdir -p $realigndir
           else
               echo "$realigndir already exists. resetting it"
               `rm -r $realigndir/*`
           fi

           if [ $skipvcall == "NO" ]
           then
               varlogdir=$outputdir/logs/variant
               if [ ! -d $varlogdir ]
               then
                   mkdir -p $varlogdir
               else
                   `rm $varlogdir/*`
               fi
               vardir=$outputdir/variant
               if [ ! -d $vardir ]
               then
                   mkdir -p $vardir
               fi

           fi


           listfiles=$outputdir/align
           qsub3=$RealignOutputLogs/qsub.realign_new
           echo "#PBS -V" > $qsub3
           echo "#PBS -A $pbsprj" >> $qsub3
           echo "#PBS -N ${pipeid}_realign_new" >> $qsub3
           echo "#PBS -l epilogue=$epilogue" >> $qsub3
           echo "#PBS -l walltime=$pbscpu" >> $qsub3
           echo "#PBS -l nodes=1:ppn=1" >> $qsub3
           echo "#PBS -o $RealignOutputLogs/log.realign_new.ou" >> $qsub3
           echo "#PBS -e $RealignOutputLogs/log.realign_new.in" >> $qsub3
           echo "#PBS -q $pbsqueue" >> $qsub3
           echo "#PBS -m ae" >> $qsub3
           echo "#PBS -M $email" >> $qsub3
           if [ `expr ${#JOBSmayo}` -gt 0 ]
           then
               echo "#PBS -W depend=afterok:$JOBSmayo" >> $qsub3
           else
               echo "#PBS -W depend=afterok:$JOBSncsa" >> $qsub3
           fi
           echo "$scriptdir/realign_new.sh $realigndir $RealignOutputLogs $listfiles $runfile $realrecalflag $RealignOutputLogs/log.realign_new.in $RealignOutputLogs/log.realign_new.ou $email $RealignOutputLogs/qsub.realign_new" >> $qsub3
           `chmod a+r $qsub3`
           realrecaljob=`qsub $qsub3`
           #`qhold -h u $realrecaljob` 
           echo $realrecaljob >> $output_logs/RECALLpbs
 
        fi
                 #`qrls -h u $alignids`
                 #`qrls -h u $mergeids`
                 #`qrls -h u $sortedmayoids`
                 #`qrls -h u $extraids`
                 #`qrls -h u $realrecaljob`

        echo "done realig/recalibrating  all bam files."
        echo `date`
	`chmod -R 770 $outputdir/`
	`chmod -R 770 $output_logs/`
fi
