#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 14 ]
then
        MSG="parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
	
set -x
echo `date`
umask 0027
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
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

echo "####################################################################################################"
echo "#####################################                       ########################################"
echo "##################################### PARSING RUN INFO FILE ########################################"
echo "##################################### AND SANITY CHECK      ########################################"
echo "####################################################################################################"

rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d "=" -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d "=" -f2 )
markduptool=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d "=" -f2 | tr '[a-z]' '[A-Z]' )
samprocessing=$( cat $runfile | grep -w SAMPROCESSING | cut -d "=" -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d "=" -f2 )
samblasterdir=$( cat $runfile | grep -w SAMBLASTERDIR | cut -d "=" -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d "=" -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d "=" -f2 )
thr=$( cat $runfile | grep -w PBSTHREADS | cut -d "=" -f2 )
alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
header=$( echo $RGparms  | tr ":" "\t" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${header}"  | tr "=" ":" )
dup_cutoff=$( cat $runfile | grep -w DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w MAP_CUTOFF | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )

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
if [ ! -d $samblasterdir ]
then
    MSG="$samblasterdir samblaster directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ ! -d $sambambadir ]
then
    MSG="$sambambadir sambamba directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ ! -d $novodir ]
then
    MSG="$novodir novocraft directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi
if [ `expr ${#markduptool}` -lt 1 ]
then
    MSG="MARKDUPLICATESTOOL=$markduptool invalid value"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
    exit 1;
fi

threads=`expr $thr "-" 1`

all_exitcodes=0

if [ $markduptool != "PICARD" ]
then
    echo -e "#######################################   CASE1: MARKDUPLICATESTOOL != PICARD ##################################################"
    echo -e "####################################### step 1: align, mark duplicates, convert to bam - on the fly ############################"

    cd $outputdir
        if [ $samprocessing == "SAMTOOLS" ]
        then
           echo `date`
           $aligndir/bwa mem $alignparms -R "${rgheader}" $ref $Rone $Rtwo | $samblasterdir/samblaster |  $samdir/samtools view -bSu -@ $threads -> ${bamprefix}.wdups
           exitcode=$?
           echo `date`
           all_exitcodes=$(( $exitcode + $all_exitcodes ))
        elif [ $samprocessing == "SAMBAMBA" ]
        then 
           echo `date`
           $aligndir/bwa mem $alignparms -R "${rgheader}" $ref $Rone $Rtwo | $samblasterdir/samblaster -o ${bamprefix}.sam.wdups
           exitcode=$?
           echo `date`
           all_exitcodes=$(( $exitcode + $all_exitcodes ))
           
           $memprof $sambambadir/sambamba view -t $thr -f bam -S ${bamprefix}.sam.wdups -o ${bamprefix}.wdups
           exitcode=$?
           echo `date`
           all_exitcodes=$(( $exitcode + $all_exitcodes ))           
        fi
        
        echo -e "cheching to see if any of the commands in this block failed"        
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem/samblaster one of more commands failed on $Rone.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit $exitcode;
        fi
        if [ ! -s ${bamprefix}.wdups ]
        then
            MSG="${bamprefix}.wdups aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit 1;
        fi
        echo `date`
fi


if [ $markduptool == "PICARD" ]
then
    echo -e "#######################################   CASE2: MARKDUPLICATESTOOL == PICARD ##################################################"
    echo -e "####################################### step 1: align, then convert sam to bam, then sort, then mark duplicates  #################"

        cd $outputdir
        echo `date`
        $aligndir/bwa mem -M $alignparms -R "${rgheader}" $ref $Rone $Rtwo |  $samdir/samtools view -bSu -@ $threads -> ${bamprefix}.tmp.bam
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem command failed on $Rone $Rtwo   exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit $exitcode;
        fi        
        if [ ! -s ${bamprefix}.tmp.bam ]
        then
            MSG="${bamprefix}.tmp.bam aligned bam file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit 1;
        fi        

        ### sometimes we may have a BAM file with NO alignmnets, just the header
#        numAlignments=$( $samdir/samtools view -c ${bamprefix}.tmp.bam ) 
        numAlignments=$( $sambambadir/sambamba view -c -t $thr ${bamprefix}.tmp.bam ) 
        echo `date`
        if [ $numAlignments -eq 0 ]
        then
            MSG="bwa mem command did not produce alignments. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            exit 1;
        fi

        $novodir/novosort --index --tmpdir $outputdir --threads $threads -m 16g --kt --compression 1 -o ${bamprefix}.sorted.bam ${bamprefix}.tmp.bam
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="novosort  failed on ${bamprefix}.tmp.bam  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit $exitcode;
        fi        
        if [ ! -s ${bamprefix}.sorted.bam ]
        then
            MSG="${bamprefix}.srted.bam aligned file not created. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
            cp $qsubfile $AlignOutputLogs/FAILEDjobs/
            exit 1;
        fi

        $samdir/samtools view -H ${bamprefix}.sorted.bam > ${bamprefix}.sorted.bam.header

        #######################################################################################################
        #######################################################################################################
        ## here we insert the code for performing QC with percent mapped reads ... as in mergeLanesPerSample.sh
        #######################################################################################################
        #######################################################################################################

        echo `date`

        if [ -s $rootdir/SAMPLENAMES_multiplexed.list ]
        then
             echo "#######################################################################################################"
             echo "########    This is a with-QC sample. We need to perform QC and write results to a txt file     #######"
             echo "#######################################################################################################"
             echo "#######################################################################################################"
             echo "########   QC rules: duplication cutoff <= $dup_cutoff AND mapped_reads cutoff >= $map_cutoff   #######"
             echo "########                                                                                        #######"
             echo "#######################################################################################################"
             ### create the output file for the results of QC and reset variables            
             QCfile=$rootdir/QC_Results.txt
             QCresult=""
             cd $outputdir
             cd ../
	     sample=`basename $PWD`
             cd $outputdir

             ### generate the files with the pertinent stats
             prestats=${bamprefix}.sorted.bam.flagstat
             poststats=${bamprefix}.sorted.bam.properlyMapped.bam.flagstat
             
             $samdir/samtools flagstat ${bamprefix}.sorted.bam > $prestats
#             $samdir/samtools view -bu -F 12 ${bamprefix}.sorted.bam > ${bamprefix}.sorted.bam.properlyMapped.bam 
             $sambambadir/sambamba view -f bam -t $thr -F "proper_pair" ${bamprefix}.sorted.bam > ${bamprefix}.sorted.bam.properlyMapped.bam 
             $samdir/samtools flagstat  ${bamprefix}.sorted.bam.properlyMapped.bam > $poststats


             ####### making sure that stats were produced 
             if [ ! -s $prestats ]
             then
 	         MSG="samtools flagstat command produced an empty file with ${bamprefix}.sorted.bam"
 	         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
 	         exit $exitcode;
             fi
             if [ ! -s $poststats ]
              then
 	         MSG="samtools flagstat command produced an empty file with ${bamprefix}.sorted.bam.properlyMapped.bam"
 	         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
 	         exit $exitcode;
             fi

  	     echo "#####################################################################################################"
 	     echo "#####################################################################################################"
             echo "             now extracting the stats of interest from the  flagstat files "
 	     echo "#####################################################################################################"
 	     echo "#####################################################################################################"
           

             tot_mapped=$( cat $poststats | grep "mapped (" | cut -d ' ' -f1 )
             tot_reads=$( cat $prestats | grep "in total" | cut -d ' ' -f1 )
             tot_dups=$( cat $poststats | grep "duplicates" | cut -d ' ' -f1 )

             #now testing if these variables have numbers

             if [ $tot_dups -eq $tot_dups 2>/dev/null ]
             then
               echo -e "ok val"
             else
	       MSG="$prestats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
	     fi
             if [ $tot_reads -eq $tot_reads 2>/dev/null ]
             then
               echo -e "ok val"
             else
	       MSG="$prestats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
	     fi

             if [ $tot_mapped -eq $tot_mapped 2>/dev/null ]
             then
               echo -e "ok val"
             else
	       MSG="$prestats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
             fi

 	     echo "#####################################################################################################"
 	     echo "#####################################################################################################"
             echo "      now calculating percentages to be used for the filtering step and making sure they have numbers"
 	     echo "#####################################################################################################"
 	     echo "#####################################################################################################"
           

             perc_dup=$(( tot_dups * 100 / tot_reads ))
             perc_mapped=$(( tot_mapped * 100 / tot_reads ))
             
             #now testing if these variables have numbers
             if [ $perc_dup -eq $perc_dup 2>/dev/null ]
             then
               echo -e "ok val"
             else
	       MSG="$stats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
             fi

             if [ $perc_mapped -eq $perc_mapped 2>/dev/null ]
             then
               echo -e "ok val"
             else
	       MSG="$stats samtools flagstat file parsed incorrectly"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	       exit $exitcode;
             fi

 	     echo "#####################################################################################################"
 	     echo "#####################################################################################################"
             echo "                     now applying the filters and populating the output reports  "
             echo "                     WE ARE NOT HALTING EXECUTION IF SAMPLE FAILS THIS QC TEST   "
 	     echo "#####################################################################################################"
 	     echo "#####################################################################################################"


             if [ $perc_dup -lt $dup_cutoff ]
             then
               echo -e "$sample passed first filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
               if [ $perc_mapped -gt $map_cutoff ]
               then
                  echo -e "$sample passed second filter percent_mapped with value $perc_mapped, minimum cutoff is $map_cutoff"
                  detail="$sample\tPASSED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	          echo -e "$detail" >> $QCfile
	          QCresult="PASSED"
               else
                  echo -e "$sample DID NOT pass second filter percent_mapped value $perc_mapped, minimum cutoff is $map_cutoff"
                  detail="$sample\tFAILED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	          echo -e "$detail"  >> $QCfile
	          QCresult="FAILED"	          
               fi
             else
               echo -e "$sample DID NOT pass first filter percent_duplicates value $perc_dup, minimum cutoff is $dup_cutoff"
               detail="$sample\tFAILED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	       echo -e "$detail"  >> $QCfile
	       QCresult="FAILED"
             fi

 	     #echo "#####################################################################################################"
             #echo "                     now using value of QCresult to decide whether to stop execution or not          "
 	     #echo "#####################################################################################################"
 	     #echo "#####################################################################################################"

             #if [ $QCresult != "PASSED" ]
             #then
             #    MSG="${bamprefix}.sorted.bam  failed QC test. Analysis will stop now"
	     #    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
             #    echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
             #    cp $qsubfile $AlignOutputLogs/FAILEDjobs/
             #    exit $exitcode;
             #fi

        fi  # end of QC for thissample
        
        echo `date`

 
        #######################################################################################################
        #######################################################################################################
        
echo -e "#######################################################################################################"
echo -e "########    Done with QC of aligned BAM, now we perform markduplicates                          #######"
echo -e "#######################################################################################################"
 
 
        $javadir/java -Xmx8g -Xms1024m -jar $picardir/MarkDuplicates.jar \
             INPUT=${bamprefix}.sorted.bam \
             OUTPUT=${bamprefix}.wdups \
             TMP_DIR=$outputdir \
             METRICS_FILE=${bamprefix}.dup.metric \
             ASSUME_SORTED=true \
             MAX_RECORDS_IN_RAM=null \
             CREATE_INDEX=true \
             VALIDATION_STRINGENCY=SILENT
             
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
                 MSG="picard-Markduplicates commands failed on ${bamprefix}.sorted.bam  exitcode=$exitcode. alignment failed"
	         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
                 echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
                 cp $qsubfile $AlignOutputLogs/FAILEDjobs/
                 exit $exitcode;
        fi
        if [ ! -s ${bamprefix}.wdups ]
        then
                 MSG="${bamprefix}.wdups aligned file not created. alignment failed"
	         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	         #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
                 echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
                 cp $qsubfile $AlignOutputLogs/FAILEDjobs/
                 exit 1;
        fi
fi # end of PICARD case

echo -e "#################      WRAP UP: sort by coordinate and index on the fly: required for creating an indexed bam, ###############"
echo -e "#################      which in turn is required for extracting alignments by chromosome                       ###############"

    $novodir/novosort --tmpdir $outputdir --threads $threads --index -m 16g --kt --compression 1 -o ${bamprefix}.wdups.sorted.bam ${bamprefix}.wdups
    exitcode=$?
    echo `date`
    if [ $exitcode -ne 0 ]
    then
         MSG="novosort command failed on ${bamprefix}.wdups.  exitcode=$exitcode alignment failed"
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" 
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
         cp $qsubfile $AlignOutputLogs/FAILEDjobs/
         exit 1;
    fi

    if [ ! -s ${bamprefix}.wdups.sorted.bam ]
    then
        MSG="${bamprefix}.wdups.sorted.bam aligned file not created. alignment failed"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
        echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
        cp $qsubfile $AlignOutputLogs/FAILEDjobs/
        exit 1;
    fi 

    $samdir/samtools view -H ${bamprefix}.wdups.sorted.bam > ${bamprefix}.wdups.sorted.bam.header
    exitcode=$?
    echo `date`
    if [ $exitcode -ne 0 ]
    then
         MSG="samtools view command failed on ${bamprefix}.wdups.sorted.bam.  exitcode=$exitcode. bwamem_pe_markduplicates stopped "
         echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
         #echo -e "program=$0 failed at line=$LINENO.\nReason=$MSG\n$LOGS" >> $AlignOutputLogs/FAILEDmessages
         cp $qsubfile $AlignOutputLogs/FAILEDjobs/
         exit $exitcode;
   fi

   
   truncate -s 0 ${bamprefix}.wdups.sorted.bam.RGline
   echo $RGparms > ${bamprefix}.wdups.sorted.bam.RGline
   echo `date`

