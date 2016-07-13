#!/bin/bash
#
# merge_bam.sh <runfile> <sample> <log.in> <log.ou> <qsub>
# 
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <sample> <log.in> <log.ou> <qsub> "
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
elog=$3
olog=$4
qsubfile=$5
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

if [ ! -s $runfile ]
then
    MSG="$runfile runfile not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
    exit 1;
fi

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 ) 
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
indeldir=$( cat $runfile | grep -w INDELDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
novocraftdir=$( cat $runfile | grep -w NOVOCRAFTDIR | cut -d '=' -f2  )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
ref_local=${refdir}/$refgenome
dbsnp_local=${refdir}/$dbSNP
indel_local=${refdir}/$indeldir
indices=$( cat $runfile | grep -w CHRNAMES | cut -d '=' -f2 | tr ':' ' ' )
outputdir=$rootdir/$SampleName

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          

SampleDir=$outputdir
RealignDir=$outputdir/realign
VarcallDir=$outputdir/variant
DeliveryDir=$rootdir/$deliverydir/$SampleName

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
tmpvariant=${SampleName}.raw.vcf                                 # name of raw variant file pre-sorting
rawvariant=${SampleName}.GATKCombineGVCF.raw.vcf                 # name of the raw variant file
outbam=${SampleName}.recalibrated.bam                            # name of the ready-for-analysis bam file

echo -e "\n\n##################################################################################" 
echo -e "##################################################################################"        
echo -e "#############                       SANITY CHECK                   ###############"
echo -e "##################################################################################"
echo -e "##################################################################################\n\n"

if [ ! -d $tmpdir ]
then
    mkdir -p $tmpdir
fi

if [ `expr ${#SampleName}` -lt 1 ]
then
    MSG="$SampleName sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1     
fi

if [ ! -d $rootdir ]
then
    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the runfile."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -d $VarcallDir ]
then
    MSG="$VarcallDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ ! -d $RealignDir ]
then
    MSG="$RealignDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi
if [ ! -d $DeliveryDir ]
then
    MSG="$DeliveryDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi



echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   MERGE BAMS BLOCK STARTS HERE  FOR              $SampleName        ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 

cd $RealignDir

echo -e "\n\n##################################################################################" 
echo -e "########### command one: skipping the merging for the bam files              #####"
echo -e "##################################################################################\n\n"

chr_bamList=$( ls -1 ${SampleName}.*.recalibrated.bam | tr "\n" " " )

if [ `expr ${#chr_bamList}` -lt 1 ]
then
    MSG="no bam files were found at $RealignDir"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

echo -e "\n\n##################################################################################" 
echo -e "########### command two: skipping the novosort part             #####"
echo -e "##################################################################################\n\n"

$novocraftdir/novosort --index --threads $thr --tmpdir $tmpdir -o $outbam  ${SampleName}.*.recalibrated.bam 

exitcode=$?

echo -e "\n\n##################################################################################" 
echo -e "########### command three: skipping the sanity check for novosort                        #####"
echo -e "##################################################################################\n\n"


echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="novosort command failed exitcode=$exitcode. merge for sample $SampleName stopped"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	 exit $exitcode;
fi
if [ -s "$outbam" ]
then     
    echo -e "### the file was created. But we are not done.     #############"
    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
    numAlignments=$( $samtoolsdir/samtools view -c $outbam ) 

    echo `date`
    if [ $numAlignments -eq 0 ]
    then
	echo -e "${SampleName}\tMERGE\tFAIL\tnovosort command did not produce alignments for $outbam\n" >> $qcfile	    
	MSG="novosort command did not produce alignments for $outbam"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
    else
	echo -e "####### $outbam seems to be in order ###########"
    fi
else 
    MSG="novosort command did not produce a file $outbam"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
    exit 1;          
fi	

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"		
echo -e "##################################################################################"        
echo -e "#############   COPY RESULTS TO DELIVERY and also to                              ############"        
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"	


echo `date`

cp $RealignDir/$outbam       $DeliveryDir

if [ ! -s $DeliveryDir/$outbam ]
then

    MSG="copy failed from $DeliveryDir/$outbam to $DeliveryDir/$outbam"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi


echo -e "${SampleName}\tVARCALLING\tPASS\tAll analyses completed successfully for this sample" >> $qcfile

echo `date`
echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"

