#!/bin/bash 
#	
#  script to perform variant calling with HaplotypeCaller ONLY
#  This module is called from within the realign module
######################################

DEBUG=0

redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 12 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi

set -x
echo `date`
ulimit -s unlimited
umask 0027
scriptfile=$0
sample=$1
outputdir=$2
inputdir=$3
inputfile=$4
chr=$5
target=$6
runfile=$7
elog=$8
olog=$9
email=${10}
qsubfile=${11}
failedlog=${12}
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;
       
if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
   exit 1;
fi

rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
ped=$( cat $runfile | grep -w PEDIGREE | cut -d '=' -f2 )
allsites=$( cat $runfile | grep -w EMIT_ALL_SITES | cut -d '=' -f2 )
snvcaller=$( cat $runfile | grep -w SNV_CALLER | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
dbsnp=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
indelfile=$( cat $runfile | grep -w INDELFILE | cut -d '=' -f2 )
genderinfo=$( cat $runfile | grep -w GENDERINFORMATION | cut -d '=' -f2 )
input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
variantcmd=$( cat $runfile | grep -w VARIANT_CMD | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
tabixdir=$( cat $runfile | grep -w TABIXDIR | cut -d '=' -f2 )

if [ `expr ${#deliveryfolder}` -lt 2 ]
then
    deliverydir=$rootdir/delivery/Vcfs
else
    deliverydir=$rootdir/$deliveryfolder/Vcfs
fi

if [ ! -d $deliverydir ]
then
    `mkdir -p $deliverydir`
fi

if [ $skipvcall != "1" -a $skipvcall != "0" -a $skipvcall != "YES" -a $skipvcall != "NO" ]
then
    MSG="Invalid value for SKIPVCALL=$skipvcall"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
else
    if [ $skipvcall == "1" -o $skipvcall == "YES" ]
    then
	echo "skipping the execution of this variant calling module"
	exit 0;
    fi
fi

if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
    input_type="WGS"
fi
if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
    input_type="WES"
fi

if [ $variantcmd == "HC" -o $variantcmd == "HAPLOTYPECALLER"  ]
then
    variantcmd="HAPLOTYPECALLER"
fi

if [ ! -d $inputdir ]
then
    MSG="$inputdir directory with realigned bams not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi
if [ ! -s ${inputdir}/${inputfile} ]
then
    MSG="${inputdir}/${inputfile} realigned bam file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi

if [ ! -d $outputdir ]
then
    mkdir -p $outputdir
fi

if [ ! -d $refdir ]
then
    MSG="$refdir reference genome directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi      
if [ ! -s $refdir/$ref ]
then
    MSG="$ref reference genome not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi
if [ ! -s $refdir/$dbsnp ]
then
    MSG="$refdir/$dbsnp dbSNP for reference genome not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi
if [ -z $snvcaller ]
then
    MSG="$snvcaller snvcaller tool was not specified in configuration file"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi

if [ $snvcaller == "GATK" -a $variantcmd != "HAPLOTYPECALLER" ]
then
    MSG="VARIANT_CMD=$variantcmd specified in configuration file. This script only runs HaplotypeCaller"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi

if [ ! -d $samdir ]
then
    MSG="$samdir samtools directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit 1;
fi
if [ -z $javadir ]
then
    MSG="A value must be specified for JAVADIR in configuration file"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi      
if [ ! -d $gatk ]
then
    MSG="$gatk GATK directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi

if [ ! -d $tabixdir ]
then
    MSG="$tabixdir tabix directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
    exit 1;
fi
echo `date`

#tabix needs to be on your path for GATK to produce *tbi files
export PATH=${PATH}:$tabixdir

set +x; echo "############################################################################################" >&2
echo "############################################################################################" >&2
echo "#############  extract $chr from  recalibrated bam                           ###############" >&2; 
echo "############################################################################################" >&2
set -x;       

cd $outputdir

infile=${chr}.${sample}.recalibrated.bam

$samdir/samtools view -bu -@ $thr -h $inputfile $chr > $infile  
  
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="split by chromosome samtools command failed exitcode=$exitcode  varcall stopped for $inputfile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi
if [ ! -s $infile ]
then
    MSG="$infile file not produced. split by chromosome samtools command failed. exitcode=$exitcode  varcall stopped for $inputfile"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog
    exit $exitcode;
fi

$sambambadir/sambamba index -t $thr $infile
exitcode=$?
echo `date`		
if [ $exitcode -ne 0 ]
then
   MSG="sambamba index command failed with $infile exitcode=$exitcode varcall stopped for $inputfile"
   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
   exit $exitcode;
fi


set +x; echo "############################################################################################" >&2
echo "############################################################################################" >&2
echo "#############  HaplotypeCaller cases: with or without gender information  ###############" >&2; 
echo "############################################################################################" >&2
set -x;       

if [ -s $genderinfo ]
then
       set +x; echo -e "\n\n" >&2;
       echo "############################################################################################" >&2       
       echo "############################################################################################" >&2
       echo "##########   HaplotypeCaller   with GENDERINFO file AND WGS input data           ###########" >&2
       echo "############################################################################################" >&2        
       echo -e "\n\n" >&2; set -x;
       
       echo `date`

       set +x; echo -e "\n\n#############    gathering gender and sample info  ###############\n\n" >&2; set -x;
	
       echo $genderinfo # Can take this line out later 
       sample_id=$sample
       gender=`cat $genderinfo | grep -v "#" | grep $sample_id | awk -F$'\t' '{print $2}'` # Can take this line out later 
       echo "Sample: $sample_id, gender: $gender" # Can take this line out later                                                                 
       echo $chr # Can take this line out later

       set +x; echo -e "\n\n#############    var-calling conditions start here  ###############\n\n" >&2; set -x;        

       ## Calling X 
       if [ $chr == "X" ];
       then
       
         set +x; echo -e "\n\n#############   calling chr=X with for TWO genders and THREE regions X_PAR1 X_PAR2 X_nonPAR ###############\n\n" >&2; set -x;        
      
	 # The X PAR coordinates can be found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt
         x_par1="X:60001-2699520"
         x_par2="X:154931044-155260560"

         ### Calling male X
         if [ $gender == "Male" ];
         then
 
		   set +x; echo -e "\n\n#############   calling chr=X gender=male X_PAR1  ###############\n\n" >&2; set -x;        

		   ploidy=2
		   site="-L "$x_par1
		   site_name="X_PAR1"
		   outfile=$infile.$site_name.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   $site \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		      eval $cmd
		   fi
		   echo `date`
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi           
 
		   set +x; echo -e "\n\n#############   calling chr=X gender=male X_PAR2  ###############\n\n" >&2; set -x;        


		   site="-L "$x_par2
		   site_name="X_PAR2"
		   outfile=$infile.$site_name.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   $site \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		     eval $cmd
		   fi
		   echo `date`
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi  


		   set +x; echo -e "\n\n#############   calling chr=X gender=male X_nonPAR  ###############\n\n" >&2; set -x;        


		   ploidy=1
		   site="-L X -XL "$x_par1" -XL "$x_par2
		   site_name="X_nonPAR"
		   outfile=$infile.$site_name.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   $site \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		     eval $cmd
		   fi
		   echo `date`  
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi  
		 fi
	 ### End calling male X  
	 echo `date`

	 ### Calling female X
         if [ $gender == "Female" ];
         then
 
		   set +x; echo -e "\n\n#############   calling chr=X gender=female   ###############\n\n" >&2; set -x;        

		   ploidy=2
		   outfile=$infile.$chr.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		     eval $cmd
		   fi
		   echo `date`
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi  
         fi
         ### End calling female X

 
       set +x; echo -e "\n\n#############   end calling chr=X   ###############\n\n" >&2; set -x;        


       echo `date`
       
       ## Calling Y
       elif [ $chr == "Y" ];
       then
       
         set +x; echo -e "\n\n#############   calling chr=Y with for TWO genders and THREE regions Y_PAR1 Y_PAR2 Y_nonPAR ###############\n\n" >&2; set -x;        
      
         # The X PAR coordinates can be found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt
         y_par1="Y:10001-2649520"
         y_par2="Y:59034050-59363566"

         ### Calling male Y
         if [ $gender == "Male" ];
         then
 
		   set +x; echo -e "\n\n#############   calling chr=Y gender=male Y_PAR1  ###############\n\n" >&2; set -x;        

		   ploidy=2
		   site="-L "$y_par1
		   site_name="Y_PAR1"
		   outfile=$infile.$site_name.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   $site \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		     eval $cmd
		   fi
		   echo `date`
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi  


		   set +x; echo -e "\n\n#############   calling chr=Y gender=male Y_PAR2  ###############\n\n" >&2; set -x;        

		   site="-L "$y_par2
		   site_name="Y_PAR2"
		   outfile=$infile.$site_name.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   $site \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		     eval $cmd
		   fi
		   echo `date`
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi  


		   set +x; echo -e "\n\n#############   calling chr=Y gender=male Y_nonPAR  ###############\n\n" >&2; set -x;        

		   ploidy=1
		   site="-L Y -XL "$y_par1" -XL "$y_par2
		   site_name="Y_nonPAR"
		   outfile=$infile.$site_name.raw.g.vcf.gz

		   cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
		   -T HaplotypeCaller \
		   -R $refdir/$ref \
		   -I ${inputdir}/${inputfile} \
		   --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
		   -gt_mode DISCOVERY \
		   -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
		   -stand_call_conf 30 \
		   -stand_emit_conf 30 \
		   --sample_ploidy $ploidy \
		   -nt 1 -nct 1 \
		   --dbsnp $refdir/$dbsnp  \
		   $site \
		   -o $outfile"

		   echo $cmd # Can take this line out later

		   if [ $DEBUG -eq 0 ]
		   then
		     eval $cmd
		   fi
		   echo `date`
		   if [ ! -s $outfile ] 
		   then
		    MSG="$outfile HaplotypeCaller file not created. vcall failed."
		    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
		    exit 1;
		   fi  
         fi
         ### End calling male Y


       set +x; echo -e "\n\n#############   end calling chr=X   ###############\n\n" >&2; set -x;        

       ## Calling MT
       elif [ $chr == "MT" ];
       then
       
         set +x; echo -e "\n\n#############   calling chr=MT   ###############\n\n" >&2; set -x;        

         ploidy=1
         site_name=$chr
         site="-L $chr"
         outfile=$infile.$site_name.raw.g.vcf.gz

         cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
         -T HaplotypeCaller \
         -R $refdir/$ref \
         -I ${inputdir}/${inputfile} \
         --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
         -gt_mode DISCOVERY \
         -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
         -stand_call_conf 30 \
         -stand_emit_conf 30 \
         --sample_ploidy $ploidy \
         -nt 1 -nct 1 \
         --dbsnp $refdir/$dbsnp  \
         $site \
         -o $outfile"

         echo $cmd # Can take this line out later

         if [ $DEBUG -eq 0 ]
         then
          eval $cmd
         fi
	 echo `date` 
         if [ ! -s $outfile ] 
         then
	    MSG="$outfile HaplotypeCaller file not created. vcall failed."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	    exit 1;
         fi  
	 
	 
         set +x; echo -e "\n\n#############   end calling chr=MT   ###############\n\n" >&2; set -x;        


       else
       
         set +x; echo -e "\n\n#############  Calling the autosomes, chr-1-chr22    ###############\n\n" >&2; set -x;        
       
         echo "Call 1->22" # Can take this line out later
         ploidy=2
         site_name=$chr
         outfile=$infile.$site_name.raw.g.vcf.gz
         
         cmd="$javadir/java -Xmx6g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
         -T HaplotypeCaller \
         -R $refdir/$ref \
         -I ${inputdir}/${inputfile} \
         --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
         -gt_mode DISCOVERY \
         -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
         -stand_call_conf 30 \
         -stand_emit_conf 30 \
         --sample_ploidy $ploidy \
         -nt 1 -nct 1 \
         --dbsnp $refdir/$dbsnp  \
         -L $chr \
         -o $outfile"

         echo $cmd # Can take this line out later

         if [ $DEBUG -eq 0 ]
         then
           eval $cmd
         fi
	 echo `date`
         if [ ! -s $outfile ] 
         then
	    MSG="$outfile HaplotypeCaller file not created. vcall failed."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	    exit 1;
         fi  
       fi
       ## End calling the autosomes

       set +x; echo -e "\n\n" >&2; 
       echo "############################################################################################" >&2
       echo "############################################################################################" >&2
       echo "####################     End HaplotypeCaller   with GENDERINFO file     ####################" >&2
       echo "############################################################################################" >&2
       echo -e "\n\n" >&2; set -x; 

else
       set +x; echo -e "\n\n" >&2; 
       echo "############################################################################################" >&2       
       echo "############################################################################################" >&2
       echo "##########   HaplotypeCaller   without GENDERINFO file                           ###########" >&2
       echo "############################################################################################" >&2        
       echo -e "\n\n" >&2; set -x;
       
       set +x; echo -e "\n\n#############  SELECT ACTIVE REGION                            ###############\n\n" >&2; set -x; 

       if [ $input_type == "WGS" -a $target != "NOTARGET" ]
       then
           set +x; echo -e "\n\n### for WES with available bed file, the target region is $target #########\n\n" >&2; set -x;         
           #target=$target
       else
           set +x; echo -e "\n\n### for any other case, the target region is $chr            ###############\n\n" >&2; set -x;         
           target=$chr
       fi
   
       outfile=$infile.$site_name.raw.g.vcf.gz
       
       cmd="$javadir/java -Xmx8g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
         -T HaplotypeCaller \
         -R $refdir/$ref \
         -I ${inputdir}/${inputfile} \
         --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
         -gt_mode DISCOVERY \
         -A Coverage -A FisherStrand -A StrandOddsRatio -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
         -stand_call_conf 30 \
         -stand_emit_conf 30 \
         -L $target \
         -nt 1 -nct 1 \
         --dbsnp $refdir/$dbsnp  \
         -o $outfile"

         echo $cmd # Can take this line out later

         if [ $DEBUG -eq 0 ]
         then
           eval $cmd
         fi
	 echo `date`
         if [ ! -s $outfile ] 
         then
	    MSG="$outfile HaplotypeCaller file not created. vcall failed."
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" >> $failedlog 
	    exit 1;
         fi 

       set +x; echo -e "\n\n" >&2; 
       echo "############################################################################################" >&2
       echo "############################################################################################" >&2
       echo "####################     End HaplotypeCaller   without GENDERINFO file     #################" >&2
       echo "############################################################################################" >&2
       echo -e "\n\n" >&2; set -x; 
         

fi

set +x; echo "############################################################################################" >&2
echo "############################################################################################" >&2
echo "#############  COPY OUTPUT FILES TO DELIVERY FOLDER                          ###############" >&2; 
echo "############################################################################################" >&2
set -x;

echo `date`

cp ${infile}.*.vcf.gz $deliverydir
cp ${infile}.*.tbi $deliverydir

echo `date`	 
