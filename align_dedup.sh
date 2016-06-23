#!/bin/bash
#
# align_dedup.sh <runfile> <SampleName> <read1> <read2> <log.in> <log.ou> <qsub>
# 
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 7 ]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <SampleName> <read1> <read2> <log.in> <log.ou> <qsub>"
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

echo -e "\n\n#####################################################################################"        
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
echo -e "#####################################################################################\n\n"        

echo -e "\n\n#####################################################################################"        
echo -e "#############             DECLARING VARIABLES                         ###############"
echo -e "#####################################################################################\n\n"        

umask 0027
set -x
echo `date`
scriptfile=$0
runfile=$1
SampleName=$2
R1=$3
R2=$4
elog=$5
olog=$6
qsubfile=$7
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
    exit 1;
fi

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 ) 
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
aligner=$( cat $runfile | grep -w BWADIR | cut -d '=' -f2  )
aligner_parms=$( cat $runfile | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
bwa_index=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
novoalign_index=$( cat $runfile | grep -w NOVOALIGNINDEX | cut -d '=' -f2 )
samtools_mod=$( cat $runfile | grep -w SAMTOOLSMODULE | cut -d '=' -f2 )
samblaster_mod=$( cat $runfile | grep -w SAMBLASTERMODULE | cut -d '=' -f2 )
sorttool_mod=$( cat $runfile | grep -w SORTMODULE | cut -d '=' -f2 )
markduplicates=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
gatk_mod=$( cat $runfile | grep -w GATKMODULE | cut -d '=' -f2 ) 
gatk_dir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
picard_mod=$( cat $runfile | grep -w PICARDMODULE | cut -d '=' -f2 )
java_mod=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
dup_cutoff=$( cat $runfile | grep -w  DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w  MAP_CUTOFF | cut -d '=' -f2 )
dbsnp_local=${refdir}/$dbSNP
outputdir=$rootdir/$SampleName

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          

module load $samtools_mod
module load $sorttool_mod


SampleDir=$outputdir
AlignDir=$outputdir/align
RealignDir=$outputdir/realign
VarcallDir=$outputdir/variant
DeliveryDir=$rootdir/$deliverydir/$SampleName

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
alignedbam=${SampleName}.nodups.bam                              # name of the aligned bam
alignedsortedbam=${SampleName}.nodups.sorted.bam                 # name of the aligned-sorted bam
dedupbam=${SampleName}.wdups.bam                                 # name of the deduplicated bam
dedupsortedbam=${SampleName}.wdups.sorted.bam                    # name of the dedup-sorted bam (output of this module)


echo -e "\n\n##################################################################################"        
echo -e "#############                       SANITY CHECK                   ###############"
echo -e "##################################################################################\n\n"        
 
if [ ! -d $tmpdir ]
then
    mkdir -p $tmpdir
fi

if [ ! -d $rootdir ]
then
    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $outputdir ]
then
    MSG="$outputdir outputdir not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $DeliveryDir ]
then
    mkdir -p $DeliveryDir
fi
if [ `expr ${#R1}` -lt 1]
then
    MSG="$R1 read one file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
elif [ ! -s $R1 ]
then
    MSG="$R1 read one file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1                
fi

if [ `expr ${#R2}` -lt 1]
then
    MSG="$R2 read two file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
elif [ ! -s $R2 ]
then
    MSG="$R2 read two  file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1                
fi

if [ `expr ${#SampleName}` -lt 1]
then
    MSG="$SampleName sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1     
else
    sID=$SampleName
    sPU=$SampleName
    sSM=$SampleName
fi
if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ] 
then
    MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${RGparms}"  | tr ":" "\t" | tr "=" ":" )


if [ `expr ${#markduplicates}` -lt 1 ]
then
    markduplicates="NOVOSORT"
fi

        

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          
echo -e "##################################################################################"        
echo -e "#############   ALIGN-DEDUPPLICATION  FOR SAMPLE $SampleName       ###############"
echo -e "##################################################################################" 
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          


echo `date` 

cd  $AlignDir

echo -e "\n\n##################################################################################"
echo -e "#############   select the dedup tool and then run the command        ############"
echo -e "#############   choices:  SAMBLASTER or NOVOSORT                      ############"
echo -e "##################################################################################\n\n"          



if [ $markduplicates == "SAMBLASTER" ]
then
	echo -e "\n\n##################################################################################"
	echo -e "##CASE1: dedup tool is $markduplciates we use a single command for align-deduplication ##"  
	echo -e "##################################################################################\n\n"

	module load $samblaster_mod

	echo -e "\n\n##################################################################################"
	echo -e "############# step one: alignment and deduplication                ############"	     
	echo -e "##################################################################################\n\n"
	
	$aligner/bwa mem $aligner_parms -t $thr -R "${rgheader}" $bwa_index $R1 $R2 | samblaster | samtools view -@ $thr -bSu -> $dedupbam 
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment-dedup step  failed for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit $exitcode;
	fi

	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step two: making sure that a file was produced with alignments    #####"
	echo -e "##################################################################################\n\n"

	if [ -s $AlignDir/$dedupbam ]
	then     
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
	    numAlignments=$( samtools view -c $AlignDir/$dedupbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for $AlignDir/$dedupbam\n" >> $qcfile	    
		MSG="bwa mem command did not produce alignments for $AlignDir/$dedupbam alignment failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    else
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########"
	    fi
	else 
	    MSG="bwa mem command did not produce a file $AlignDir/$dedupbam alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi       


	echo -e "\n\n##################################################################################"
	echo -e "#############  step three: sort                                      ############"
	echo -e "##################################################################################\n\n"

	novosort --index --tmpdir $tmpdir --threads $thr -m 16g --compression 1 -o $dedupsortedbam $dedupbam
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="align-sorting step failed for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"       
	    exit $exitcode;
	fi

	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step four: making sure that a file was produced with alignments    #####"
	echo -e "##################################################################################\n\n"

	if [ -s $AlignDir/$dedupsortedbam ]
	then     
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
	    numAlignments=$( samtools view -c $AlignDir/$dedupsortedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for $AlignDir/$dedupsortedbam\n" >> $qcfile	    
		MSG="novosort command did not produce a file for $AlignDir/$dedupsortedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		exit 1;
	    else
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########"
	    fi
	else 
	    MSG="novosort command did not produce a file $AlignDir/$dedupsortedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi       	
	
	module unload $samblaster_mod               

	echo -e "\n\n##################################################################################"
	echo -e "#############      END SAMBLASTER BLOCK                               ############"
	echo -e "##################################################################################\n\n"             

elif  [ $markduplicates == "NOVOSORT" ]
then
	echo -e "\n\n##################################################################################"
	echo -e "CASE2: dedup tool is NOVOSORT. one cmd for align and one for dedup-sort   ########"  
	echo -e "##################################################################################\n\n"


	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step one: alignment                                 ############"
	echo -e "##################################################################################\n\n"
        

        if [ $aligner == "BWA" ]
        then
	   bwa mem $aligner_parms -t $thr -R "${rgheader}" $ref_local $R1 $R2 | samtools view -@ $thr -bSu -> $alignedbam 
	   exitcode=$?
	   echo `date`
	   if [ $exitcode -ne 0 ]
	   then
	       MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	       exit $exitcode;
	   fi
        elif [ $aligner == "NOVOALIGN" ]
           novoalign $aligner_parms  -c $thr -d ${novoalign_index} -f $R1 $R2 | samtools view -@ $thr -bS - > $alignedbam
           exitcode=$?
           echo `date`
           if [ $exitcode -ne 0 ]
           then
               MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
               echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
               exit $exitcode;
           fi
        fi   

	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step two:making sure that a file was produced with alignments     #####"
	echo -e "##################################################################################\n\n"

	if [ -s $AlignDir/$alignedbam ]
	then            
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"

	    numAlignments=$( samtools view -c $AlignDir/$alignedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\taligner command did not produce alignments for $AlignDir/$alignedbam\n" >> $qcfile	    
		MSG="aligner command did not produce alignments for $AlignDir/$alignedbam alignment failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"      
		exit 1;
	    else
		echo -e "####### $AlignDir/$alignedbam seems to be in order ###########"
	    fi
	else 
	    MSG="aligner command did not produce a file $AlignDir/$alignedbam alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi       
	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step three: sort + dedup + indexing                       ############"
	echo -e "##################################################################################\n\n"

	novosort -markduplicates  -t $tmpdir -m 16G -c 8 -i -o $dedupsortedbam $alignedbam

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during sorting-deduplication for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit $exitcode
	fi 
	
	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step four: making sure that a file was produced with alignments #######"
	echo -e "##################################################################################\n\n"
	
	if [ -s $AlignDir/$dedupsortedbam ]
	then
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
	    
	    numAlignments=$( samtools view -c $AlignDir/$dedupsortedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for $AlignDir/$dedupsortedbam\n" >> $qcfile	    
		MSG="novosort command did not produce a file for $AlignDir/$dedupsortedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		exit 1;
	    else
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########"
	    fi
	else 
	    MSG="novosort command did not produce a file $AlignDir/$dedupsortedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi   
	
	echo -e "\n\n##################################################################################"
	echo -e "#############      END NOVOSRT  BLOCK                                 ############"
	echo -e "##################################################################################\n\n"             
	

elif  [ $markduplicates == "PICARD" ]
then
	echo -e "\n\n##################################################################################"
	echo -e "CASE2: dedup tool is PICARD. one cmd for align, one for sort, one for dedup   ########"  
	echo -e "##################################################################################\n\n"


        module unload java
        module load $java_mod       
        module load $picard_mod
        
 
        
	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step one: alignment                                 ############"
	echo -e "##################################################################################\n\n"

	bwa mem $aligner_parms -t $thr -R "${rgheader}" $ref_local $R1 $R2 | samtools view -@ $thr -bSu -> $alignedbam 
	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	    exit $exitcode;
	fi   

	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step two:making sure that a file was produced with alignments     #####"
	echo -e "##################################################################################\n\n"

	if [ -s $AlignDir/$alignedbam ]
	then            
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"

	    numAlignments=$( samtools view -c $AlignDir/$alignedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for $AlignDir/$alignedbam\n" >> $qcfile	    
		MSG="bwa mem command did not produce alignments for $AlignDir/$alignedbam alignment failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"      
		exit 1;
	    else
		echo -e "####### $AlignDir/$alignedbam seems to be in order ###########"
	    fi
	else 
	    MSG="bwa mem command did not produce a file $AlignDir/$alignedbam alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi       
	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step three: sort                                        ############"
	echo -e "##################################################################################\n\n"

	novosort -t $tmpdir -m 16G -c 8 -i -o $alignedsortedbam $alignedbam

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during sorting for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit $exitcode
	fi 


	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step four: dedup                                        ############"
	echo -e "##################################################################################\n\n"


java -Xmx8g -Djava.io.tmpdir=$tmpdir -jar $picardir/picard.jar  MarkDuplicates \
INPUT=$alignedsortedbam OUTPUT=$dedupsortedbam TMP_DIR=$tmpdir \
ASSUME_SORTED=true MAX_RECORDS_IN_RAM=null CREATE_INDEX=true \
METRICS_FILE=${SampleName}.picard.metrics \
VALIDATION_STRINGENCY=SILENT

	exitcode=$?
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during deduplication for sample $SampleName exitcode=$exitcode."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
            exit $exitcode
	fi 

	
	echo -e "\n\n##################################################################################"	     
	echo -e "#############  step five: making sure that a file was produced with alignments #######"
	echo -e "##################################################################################\n\n"
	
	if [ -s $AlignDir/$dedupsortedbam ]
	then
	    echo -e "### the file was created. But we are not done.     #############"
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
	    
	    numAlignments=$( samtools view -c $AlignDir/$dedupsortedbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tpicard command did not produce a file for $AlignDir/$dedupsortedbam\n" >> $qcfile	    
		MSG="novosort command did not produce a file for $AlignDir/$dedupsortedbam"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		exit 1;
	    else
		echo -e "####### $AlignDir/$dedupbam seems to be in order ###########"
	    fi
	else 
	    MSG="picard command did not produce a file $AlignDir/$dedupsortedbam"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"       
	    exit 1;          
	fi   


        `module unload $picard_mod`
        
	echo -e "\n\n##################################################################################"
	echo -e "#############      END PICARD  BLOCK                                 ############"
	echo -e "##################################################################################\n\n"             

	
else
	MSG="unrecognized deduplication tool $markduplicates"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;        

fi


        
echo -e "\n\n##################################################################################"
echo -e "#############     END ALIGNMENT-DEDUPLICATION BLOCK                   ############"
echo -e "##################################################################################\n\n"

echo `date`

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          
echo -e "##################################################################################"        
echo -e "########   ALIGNMENT QC TEST   FOR SAMPLE $SampleName              ###############"
echo -e "########   QC rule1: duplication cutoff <= $dup_cutoff             ###############"
echo -e "########   QC rule2: mapped_reads cutoff >= $map_cutoff            ###############"
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"
        
     

echo -e "\n\n##################################################################################"	     
echo -e "#############  step one: generating the relevant file with flagstat       ############"
echo -e "##################################################################################\n\n"

flagstats=${dedupsortedbam}.flagstats

echo `date`             
samtools flagstat $dedupsortedbam > $flagstats
echo `date`

echo -e "\n\n##################################################################################"	     
echo -e "#############  step two: sanity check                                 ############"
echo -e "##################################################################################\n\n" 

if [ ! -s $flagstats ]
then
	 echo -e "${SampleName}\tQCTEST\tFAIL\tsamtools flagstat command produced an empty file $flagstats\n" >> $qcfile
	 MSG="samtools flagstat command produced an empty file  $flagstats"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	 exit $exitcode;
fi

echo -e "\n\n##################################################################################"	     
echo -e "#############  step three: parsing the file and grabbing stats for the QC test     ###"
echo -e "##################################################################################\n\n"            


tot_mapped=$( cat $flagstats | grep "mapped (" | cut -d ' ' -f1 )
tot_reads=$( cat $flagstats | grep "in total" | cut -d ' ' -f1 )
tot_dups=$( cat $flagstats | grep "duplicates" | cut -d ' ' -f1 )

#now testing if these variables are numeric and have numbers

if [ $tot_dups -eq $tot_dups 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi
if [ $tot_reads -eq $tot_reads 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

if [ $tot_mapped -eq $tot_mapped 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

           
echo -e "\n\n##################################################################################"	     
echo -e "#############  step four: calculating stats according to QC rules                  ###"
echo -e "##################################################################################\n\n"            


perc_dup=$(( tot_dups * 100 / tot_reads ))
perc_mapped=$(( tot_mapped * 100 / tot_reads ))

#now testing if these variables have numbers

if [ $perc_dup -eq $perc_dup 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

if [ $perc_mapped -eq $perc_mapped 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstats samtools flagstat file parsed incorrectly"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit $exitcode;
fi

echo -e "\n\n##################################################################################"	     
echo -e "#############  step five: applying the  QC rules                                  ###"
echo -e "##################################################################################\n\n"            

if [ $perc_dup -lt $dup_cutoff ]
then
	echo -e "$#####  sample passed first filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
	
	if [ $perc_mapped -gt $map_cutoff ]
	then
	        echo -e "##### $sample passed second filter map_cutoff with value $perc_mapped, minimum cutoff is $map_cutoff"	
	        echo -e "${SampleName}\tQCTEST\tPASS\trule1 ok: percent_duplication=$perc_dup <? duplication_cutoff=$dup_cutoff\trule2 ok: percent_mapped=$perc_mapped >? mapping_cutoff=$map_cutoff" >> $qcfile
	else
	        echo -e "##### $sample DID NOT pass second filter map_cutoff with value $perc_mapped, minimum cutoff is $map_cutoff"	
	        echo -e "${SampleName}\tQCTEST\tFAIL\trule1 ok: percent_duplication=$perc_dup <? duplication_cutoff=$dup_cutoff\trule2 notok: percent_mapped=$perc_mapped >? mapping_cutoff=$map_cutoff" >> $qcfile	          
	fi
else
	echo -e "$#####  sample DID NOT pass first filter percent_duplicates with value $perc_dup, maximum cutoff is $dup_cutoff"
	echo -e "${SampleName}\tQCTEST\tFAIL\trule1 not ok: percent_duplication=$perc_dup <? duplication_cutoff=$dup_cutoff\trule2 not evaluated: percent_mapped=$perc_mapped >? mapping_cutoff=$map_cutoff" >> $qcfile
fi


echo -e "\n\n##################################################################################"
echo -e "#############       END QC TEST                                       ############"        
echo -e "##################################################################################\n\n"

echo `date`

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"		
echo -e "##################################################################################"        
echo -e "#############   WRAP UP                                               ############"        
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"	

module unload $samtools_mod
module unload $sorttol_mod
echo `date`

### perhaps this bam file is not necessary in the delivery folder           
### cp $AlignDir/${SampleName}.wdups.sorted.bam          $DeliveryDir   


MSG="ALIGNMENT-DEDUPLICATION for $SampleName finished successfully"
echo -e "program=$scriptfile at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"

echo `date`

echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"
