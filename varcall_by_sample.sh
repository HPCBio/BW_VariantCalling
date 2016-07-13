#!/bin/bash
#
# varcall_by_sample.sh <runfile> <sample> <log.in> <log.ou> <qsub>
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
elog=$4
olog=$5
qsubfile=$6
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
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 ) 
ref_local=${refdir}/$refgenome
dbsnp_local=${refdir}/$dbSNP
indel_local=${refdir}/$indeldir

outputdir=$rootdir/$SampleName

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          


SampleDir=$outputdir
AlignDir=$outputdir/align
RealignDir=$outputdir/realign
VarcallDir=$outputdir/variant
DeliveryDir=$rootdir/$deliverydir/$SampleName

qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
inputbam=$AlignDir/${SampleName}.wdups.sorted.bam                # name of the bam file that align-dedup produced
dedupsortedbam=${SampleName}.$chr.wdups.sorted.bam               # name of the aligned file -one chr out of inputbam
realignedbam=${SampleName}.$chr.realigned.bam                    # name of the realigned file
recalibratedbam=${SampleName}.$chr.recalibrated.bam              # name of the recalibrated file



qcfile=$rootdir/$deliverydir/docs/QC_test_results.txt            # name of the txt file with all QC test results
tmpvariant=${SampleName}.raw.vcf                                 # name of raw variant file pre-sorting
rawvariant=${SampleName}.GATKCombineGVCF.raw.vcf                 # name of the raw variant file (hc output could go to this program)
outbam=$RealignDir/${SampleName}.recalibrated.bam                # name of the ready-for-analysis bam file (input to hc)
rawvariant=${SampleName}.$chr.raw.g.vcf                          # name of the raw variant file (output from hc)

echo -e "\n\n##################################################################################" 
echo -e "##################################################################################"        
echo -e "#############                       SANITY CHECK                   ###############"
echo -e "##################################################################################"
echo -e "##################################################################################\n\n"

if [ ! -d $tmpdir ]
then
    mkdir -p $tmpdir
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
if [ ! -d $DeliveryDir ]
then
    MSG="$DeliveryDir directory not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ ! -s $inputbam ]
then
    MSG="$inputbam aligned-dedup file not found"
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

echo -e "\n\n##################################################################################"
echo -e "########### PREP WORK: split aligned.bam by chr and grab indels for chr            ###" 
echo -e "##################################################################################\n\n"

echo -e "\n### ploidy variable, its value is 2 for all chr except mitochondrial               ###\n"

if [ $chr == "M" ]
then
    ploidy=1
else
    ploidy=2
fi

echo -e "\n### split aligned.bam by region with  $inputbam and chr=$chr                       ###\n"     

echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"	
echo -e "##################################################################################"        
echo -e "#############    GATK VARIANT CALLING   FOR SAMPLE $SampleName     ###########"
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"

echo `date`        

cd $VarcallDir

echo -e "\n\n##################################################################################"            
echo -e "########### command one: executing GATK HaplotypeCaller command         ##########" 
echo -e "##################################################################################\n\n"



$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
	 -T HaplotypeCaller \
	 -R $ref_local \
	 --dbsnp $dbsnp_local \
	 -I $RealignDir/$recalibratedbam \
	 --emitRefConfidence GVCF \
	 -gt_mode DISCOVERY \
	 -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
	 -stand_call_conf 30 \
	 -stand_emit_conf 30 \
	 --sample_ploidy $ploidy \
	 -nt 1 -nct $thr \
	 -o $rawvariant


exitcode=$?
echo `date`

if [ $exitcode -ne 0 ]
then
	MSG="GATK HaplotypeCaller command failed exitcode=$exitcode for $rawvariant"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	exit $exitcode;
fi
if [ ! -s $rawvariant ]
then
	echo -e "${SampleName}\tVCALL\tFAIL\tHaplotypeCaller command did not produce results for $rawvariant\n" >> $qcfile	    
	MSG="GATK HaplotypeCaller command did not produce results for $rawvariant"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

echo -e "\n\n##################################################################################"
echo -e "#############       END VARIANT CALLING BLOCK                         ############"        
echo -e "##################################################################################\n\n"


echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"		
echo -e "##################################################################################"        
echo -e "#############   WRAP UP                                               ############"        
echo -e "##################################################################################"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"	

echo `date`
 
# we will merge all variant files for this sample and copy that file to delivery
#cp $VarcallDir/rawvariant=${SampleName}.raw.vcf $DeliveryDir  

echo `date`
echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"

