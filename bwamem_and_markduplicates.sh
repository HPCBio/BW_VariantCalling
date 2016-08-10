#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 11 ]
then
        MSG="parameter mismatch"
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
	
set -x
echo `date`
umask 0027
scriptfile=$0
outputdir=$1
inputfiles=$2
alignedbam=$3
outputbam=$4
RGparms=$5
runfile=$6
elog=$7
olog=$8
email=$9
qsubfile=${10}
failedlog=${11}

LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"



if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

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
alignerdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
refindex=$( cat $runfile | grep -w BWAMEMINDEX | cut -d '=' -f2 )
alignparms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 | tr " " "_" )
memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
dup_cutoff=$( cat $runfile | grep -w DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w MAP_CUTOFF | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
runqctest=$( cat $runfile | grep -w RUNDEDUPQC | cut -d '=' -f2 )
analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
QCfile=$rootdir/QC_Results.txt


if [ $runqctest == "YES" -a ! -s $qc_result ]
then
    MSG="$qc_result file for QC test results was not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi

if [ ! -d $alignerdir ]
then
    MSG="$alignerdir directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $refindexed ]
then
    MSG="$refindexed directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $samdir ]
then
    MSG="$samdir samtools directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $novodir ]
then
    MSG="$novodir novocraft directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $picardir ]
then
    MSG="$picardir picard directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi


if [ ! -d $samblasterdir ]
then
    MSG="$samblasterdir samblaster directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $sambambadir ]
then
    MSG="$sambambadir sambamba directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ ! -d $novodir ]
then
    MSG="$novodir novocraft directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ `expr ${#markduptool}` -lt 1 ]
then
    MSG="MARKDUPLICATESTOOL=$markduptool invalid value"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi

threads=`expr $thr "-" 1`

header=$( echo $RGparms  | tr ":" "\t" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${header}"  | tr "=" ":" )
threads=`expr $thr "-" 1`
fqfiles=$( echo $inputfiles | tr ":" " " )
prefix=${alignedbam%.bam}
samfile=${prefix}.sam
unsortedbam=${prefix}.temp.unsorted.bam
sortedbam=${prefix}.temp.sorted.bam
dedupbam=${prefix}.wdups.unsorted.bam
all_exitcodes=0

set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    Run bwa mem -> conversion to bam -> dedup                                           #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $outputdir
echo `date`

if [ $markduptool != "PICARD" ]
then
	set +x; echo -e "\n\n" >&2; 
	echo -e "#######################################   CASE1: MARKDUPLICATESTOOL != PICARD ##################################################" >&2
	echo -e "####################################### step 1: align, mark duplicates, convert to bam - on the fly ############################" >&2
	echo -e "\n\n" >&2; set -x;


        if [ $samprocessing == "SAMTOOLS" ]
        then
           echo `date`
           $alignerdir/bwa mem $alignparms -t $threads -R "${rgheader}" $refindex $fqfiles | $samblasterdir/samblaster |  $samdir/samtools view -bSu -@ $threads -> $dedupbam
           exitcode=$?
           echo `date`
           all_exitcodes=$(( $exitcode + $all_exitcodes ))
        elif [ $samprocessing == "SAMBAMBA" ]
        then 
           echo `date`
           $alignerdir/bwa mem $alignparms -t $threads  -R "${rgheader}" $refindex $fqfiles | $samblasterdir/samblaster -o $samfile
           exitcode=$?
           echo `date`
           all_exitcodes=$(( $exitcode + $all_exitcodes ))
           
           $memprof $sambambadir/sambamba view -t $thr -f bam -S $samfile -o $dedupbam
           exitcode=$?
           echo `date`
           all_exitcodes=$(( $exitcode + $all_exitcodes ))           
        fi


        set +x; echo -e "\n\n########################### checking to see if any of the commands in this block failed #################\n\n" >&2; set -x;
        
      
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem/samblaster one of more commands failed on $Rone.  exitcode=$exitcode. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit $exitcode;
        fi
        if [ ! -s $dedupbam ]
        then
            MSG="$dedupbam aligned-dedup file not created. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit 1;
        fi
        
        ### sometimes we may have a BAM file with NO alignmnets, just the header
 
        #numAlignments=$( $sambambadir/sambamba view -c -t $thr $dedupbam ) 
        numAlignments=$( $samdir/samtools view -c -@ $thr $dedupbam )
        echo `date`
        if [ `expr ${#numAlignments}` -lt 1 ]
        then
            MSG="bwa mem command did not produce alignments. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit 1;
        fi
        
        
        echo `date`
elif [ $markduptool == "PICARD" ]
then
        set +x; echo -e "\n\n" >&2; 
        echo -e "#######################################   CASE2: MARKDUPLICATESTOOL == PICARD   ##################################################" >&2
        echo -e "####################################### step 1: align, then convert sam to bam, then sort, then mark duplicates  #################" >&2
        echo -e "\n\n" >&2; set -x;
        
        cd $outputdir
        echo `date`
        $alignerdir/bwa mem $alignparms -t $threads -R "${rgheader}" $refindex $fqfiles |  $samdir/samtools view -bSu -@ $threads -> $unsortedbam
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="bwa mem command failed on $fqfiles   exitcode=$exitcode. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit $exitcode;
        fi        
        if [ ! -s $unsortedbam ]
        then
            MSG="$unsortedbam aligned bam file not created. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit 1;
        fi        

        ### sometimes we may have a BAM file with NO alignmnets, just the header
 
        numAlignments=$( $samdir/samtools view -c -@ $thr $unsortedbam ) 
        echo `date`
        if [ `expr ${#numAlignments}` -lt 1 ]
        then
            MSG="bwa mem command did not produce alignments. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit 1;
        fi

        # sorting the bam file
        $novodir/novosort --index --tmpdir $outputdir --threads $threads -m 16g --compression 1 -o $sortedbam $unsortedbam
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="novosort  failed on $unsortedbam  exitcode=$exitcode. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit $exitcode;

        fi        
        if [ ! -s $sortedbam ]
        then
            MSG="$sortedbam aligned file not created. alignment failed"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
            exit 1;
        fi

        # creating header file
        
        $samdir/samtools view -H $sortedbam > $sortedbam.header


        # running the deduplication command
        #$javadir/java -Xmx8g -Xms1024m -jar $picardir/MarkDuplicates.jar \

        java -Xmx8g -jar $picardir/picard.jar MarkDuplicates \
             INPUT=$sortedbam \
             OUTPUT=$dedupbam \
             TMP_DIR=$outputdir \
             METRICS_FILE=${sortedbam}.dup.metric \
             ASSUME_SORTED=true \
             MAX_RECORDS_IN_RAM=null \
             CREATE_INDEX=true \
             VALIDATION_STRINGENCY=SILENT
             
        exitcode=$?
        echo `date`
        if [ $exitcode -ne 0 ]
        then
                 MSG="picard-Markduplicates commands failed on $sortedbam  exitcode=$exitcode. alignment failed"
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
                 exit $exitcode;
        fi
        if [ ! -s $dedupbam ]
        then
                 MSG="$dedupbam aligned file not created. alignment failed"
                 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
                 exit 1;
        fi

fi

set +x; echo -e "\n\n" >&2;
echo -e "#######################################################################################################" >&2
echo -e "########    Sort by coordinate to make GATK happy                                           #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

$novodir/novosort --tmpdir $outputdir --threads $threads --index -m 16g --compression 1 -o $outputbam $dedupbam
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="novosort command failed on $dedupbam  exitcode=$exitcode alignment failed"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi

if [ ! -s $outputbam ]
then
	MSG="$outputbam file not created. alignment failed"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit 1;
fi 

$samdir/samtools view -H $outputbam > $outputbam.header
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
	MSG="samtools view command failed on $outputbam  exitcode=$exitcode. bwamem_pe_markduplicates stopped "
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
	exit $exitcode;
fi


truncate -s 0 $outputbam.RGline
echo $RGparms > $outputbam.RGline
echo `date`


set +x; echo -e "\n\n" >&2; 
echo -e "######################################################################################################" >&2
echo -e #######################################################################################################" >&2
echo -e ##      performing QC with percent mapped reads                       #################################" >&2
echo -e #######################################################################################################" >&2
echo -e #######################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

echo `date`

if [ $runqctest == "YES" ]
then
     set +x; echo -e "\n\n" >&2;         
     echo "#######################################################################################################" >&2
     echo "########    This is a with-QC sample. We need to perform QC and write results to a txt file     #######" >&2
     echo "#######################################################################################################" >&2
     echo "#######################################################################################################" >&2
     echo "########   QC rules: duplication cutoff <= $dup_cutoff AND mapped_reads cutoff >= $map_cutoff   #######" >&2
     echo "########                                                                                        #######" >&2
     echo "#######################################################################################################" >&2
     echo -e "\n\n" >&2; set -x;      

     
     ### Need to know the sample name           
     
     cd $outputdir/..
     sample=`basename $PWD`
     cd $outputdir

     ### generate the files with the pertinent stats
     
     bamstats=${outputbam}.flagstat

     $samdir/samtools flagstat $outputbam > $bamstats

     ####### making sure that stats were produced 
     
     if [ ! -s $bamstats ]
     then
	 MSG="samtools flagstat command produced an empty file for ${outputbam}"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>> $failedlog
	 exit $exitcode;
     fi

     set +x; echo -e "\n\n" >&2;
     echo "#####################################################################################################" >&2
     echo "#####################################################################################################" >&2
     echo "             now extracting the stats of interest from the  flagstat file                           " >&2
     echo "#####################################################################################################" >&2
     echo "#####################################################################################################" >&2
     echo -e "\n\n" >&2; set -x;            

     tot_mapped=$( cat $bamstats | grep "mapped (" | cut -d ' ' -f1 )
     tot_reads=$( cat $bamstats | grep "in total" | cut -d ' ' -f1 )
     tot_dups=$( cat $bamstats | grep "duplicates" | cut -d ' ' -f1 )

     #now testing if these variables have numbers

     if [ $tot_dups -eq $tot_dups 2>/dev/null ]
     then
       echo -e "ok val"
     else
       MSG="$bamstats samtools flagstat file parsed incorrectly"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>> $failedlog
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
       exit $exitcode;
     fi
     if [ $tot_reads -eq $tot_reads 2>/dev/null ]
     then
       echo -e "ok val"
     else
       MSG="$bamstats samtools flagstat file parsed incorrectly"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>> $failedlog
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
       exit $exitcode;
     fi

     if [ $tot_mapped -eq $tot_mapped 2>/dev/null ]
     then
       echo -e "ok val"
     else
       MSG="$bamstats samtools flagstat file parsed incorrectly"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>> $failedlog
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
       exit $exitcode;
     fi

     set +x; echo -e "\n\n" >&2;
     echo "#####################################################################################################" >&2
     echo "#####################################################################################################" >&2
     echo "      now calculating percentages to be used for the filtering step and making sure they have numbers" >&2
     echo "#####################################################################################################" >&2
     echo "#####################################################################################################" >&2
     echo -e "\n\n" >&2; set -x;            

     perc_dup=$(( tot_dups * 100 / tot_reads ))
     perc_mapped=$(( tot_mapped * 100 / tot_reads ))

     #now testing if these variables have numbers
     if [ $perc_dup -eq $perc_dup 2>/dev/null ]
     then
       echo -e "ok val"
     else
       MSG="$bamstats samtools flagstat file parsed incorrectly"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>> $failedlog
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
       exit $exitcode;
     fi

     if [ $perc_mapped -eq $perc_mapped 2>/dev/null ]
     then
       echo -e "ok val"
     else
       MSG="$bamstats samtools flagstat file parsed incorrectly"
       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>> $failedlog
       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
       exit $exitcode;
     fi

     set +x; echo -e "\n\n" >&2;
     echo "#####################################################################################################" >&2
     echo "#####################################################################################################" >&2
     echo "                     now applying the filters and populating the output reports  " >&2
     echo "                     WE ARE NOT HALTING EXECUTION IF SAMPLE FAILS THIS QC TEST   " >&2
     echo "#####################################################################################################" >&2
     echo "#####################################################################################################" >&2
     echo -e "\n\n" >&2; set -x;

     if [ $perc_dup -lt $dup_cutoff ]
     then
       echo -e "$sample passed first filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
       if [ $perc_mapped -gt $map_cutoff ]
       then
	  echo -e "$sample passed second filter percent_mapped with value $perc_mapped, minimum cutoff is $map_cutoff"
	  detail="$sample\tPASSED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	  echo -e "$detail" >> $QCfile
       else
	  echo -e "$sample DID NOT pass second filter percent_mapped value $perc_mapped, minimum cutoff is $map_cutoff"
	  detail="$sample\tFAILED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
	  echo -e "$detail"  >> $QCfile	          
       fi
     else
       echo -e "$sample DID NOT pass first filter percent_duplicates value $perc_dup, minimum cutoff is $dup_cutoff"
       detail="$sample\tFAILED\tpercent_duplication=$perc_dup\tduplication_cutoff=$dup_cutoff\tpercent_mapped=$perc_mapped\tmapping_cutoff=$map_cutoff"
       echo -e "$detail"  >> $QCfile
     fi

fi  # end of QC for thissample

echo `date`

set +x; echo -e "\n\n" >&2;        
echo -e "#######################################################################################################" >&2
echo -e "########    Send files to delivery if analysis ends here                                        #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x; 

if [ $analysis == "ALIGNMENT" -o $analysis == "ALIGN" -o $analysis == "ALIGN_ONLY" ]
then
	set +x; echo -e "\n ###### ANALYSIS = $analysis ends here. Wrapping up and quitting\n" >&2; set -x;

	deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
	rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )

	cd $outputdir/..
	sample=`basename $PWD`
	cd $outputdir

	if [ `expr ${#deliveryfolder}` -lt 2 ]
	then
	     delivery=$rootdir/delivery
	else
	     delivery=$rootdir/$deliveryfolder
	fi
	if [ ! -d $delivery/$sample ]
	then
	    mkdir -p $delivery/$sample
	fi
	
	cp $outputbam $delivery/$sample
	
	echo `date`

fi

set +x; echo -e "\n\n" >&2;        
echo -e "#######################################################################################################" >&2
echo -e "########    Done with QC of aligned BAM, now exiting                                            #######" >&2
echo -e "#######################################################################################################" >&2
echo -e "\n\n" >&2; set -x; 




