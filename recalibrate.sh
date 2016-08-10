#!/bin/bash
#	
#  script to 1) recalibrate realigned file(s) and 2) run QC with samtools calmd and VerifyBamId 
########################################################
redmine=hpcbio-redmine@igb.illinois.edu
set -x
if [ $# != 10 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi					

echo `date`
umask 0027
scriptfile=$0
realrecaldir=$1
sample=$2
inputfiles=$3 
outputfile=$4    
RGparms=$5
runfile=$6
elog=$7
olog=$8
email=$9
qsubfile=${10}
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
RealignOutputLogs=`dirname $elog`
infile=`basename $inputfile`

if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi

set +x; echo -e "\n\n" >&2; 
echo -e "####################################################################################################" >&2
echo -e "#####################################                       ########################################" >&2
echo -e "##################################### PARSING RUN INFO FILE ########################################" >&2
echo -e "##################################### AND SANITY CHECK      ########################################" >&2
echo -e "####################################################################################################" >&2
echo -e "\n\n" >&2; set -x;

input_type=$( cat $runfile | grep -w INPUTTYPE | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )        
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
threads=$( cat $runfile | grep -w PBSTHREADS | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
sambambadir=$( cat $runfile | grep -w SAMBAMBADIR | cut -d '=' -f2 )
gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
novodir=$( cat $runfile | grep -w NOVODIR | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
indelfile=$( cat $runfile | grep -w INDELFILE | cut -d '=' -f2 )
targetfile=$( cat $runfile | grep -w TARGETREGIONSFILE | cut -d '=' -f2 )
omnisites=$( cat $runfile | grep -w OMNISITES | cut -d '=' -f2 ) 
recalibrator=$( cat $runfile | grep -w RECALIBRATOR | cut -d '=' -f2 )
runverify=$( cat $runfile | grep -w RUNVERIFYBAM | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
outputrootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
thr=`expr $threads "-" 1`
RGparms=$( echo $RGparms | tr "::" ":" | sed "s/:/\tRG/g" | sed "s/ID\=/RGID\=/" )
freemix_cutoff=$( cat $runfile | grep -w FREEMIX_CUTOFF | cut -d '=' -f2 )
qc_result=$outputrootdir/QC_Results.txt

if [ ! -s $inputfiles ]
then
    MSG="$inputfiles input file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;

fi


if [ ! -d $realrecaldir ]
then
    MSG="$realrecaldir realign directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi

if [ ! -d $samdir ]
then
    MSG="$samdir samtools directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi
if [ ! -d $gatk ]
then
    MSG="$gatk GATK directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi

if [ -z $javadir ]
then
    MSG="Value for JAVADIR must be specified in configuration file"
    echo -e "program=$scriptfile WARNING at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    #exit 1;
fi

if [ ! -d $refdir ]
then
    MSG="$refdir reference genome directory not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi      
if [ ! -s $refdir/$ref ]
then
    MSG="$ref reference genome not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi
if [ ! -s $refdir/dbSNP ]
then
    MSG="$dbSNP DBSNP file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi
if [ ! -s $indelfile ]
then
    MSG="$indelfile INDELS file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_realign.${infile}.Anisimov.msg    
    exit 1;
fi


if [ $input_type == "GENOME" -o $input_type == "WHOLE_GENOME" -o $input_type == "WHOLEGENOME" -o $input_type == "WGS" ]
then
     input_type="WGS"
fi
if [ $input_type == "EXOME" -o $input_type == "WHOLE_EXOME" -o $input_type == "WHOLEEXOME" -o $input_type == "WES" ]
then
     input_type="WES"
fi

if [ $runverify == "NO" -o $runverify == "0" ]
then
     runverify="NO"
else
     runverify="YES"
fi

if [ $runverify == "YES" -a ! -s $omnisites ]
then
   MSG="OMNISITES=$omnisites file not found. A file must be specified if this parameter has been specified too RUNVERIFYBAM=$runverify"
   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
   exit 1;
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP1: define variables                                     ###" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

cd $realrecaldir

realignedfiles=$( cat $inputfiles | tr ":" " " )
realignedMergedFile=${sample}.realigned.merged.bam
recalreport=${sample}.recal_report.grp
covariatsFile=${sample}.recal_data.csv 
recalibratedFile=${sample}.recalibrated.bam
verifiedOut=${sample}.VerifyBamIdOutput

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP2: combine all realigned bams from regions into one file###" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

set +x; echo -e "\n\n###### Let us assume that all regions come from the same sample and that we do not need to redefine the RG line ########\n\n" >&2; set -x;
    
$novodir/novosort --index --threads $thr --tmpdir $realrecaldir -m 16g -o $realignedMergedFile  $realignedfiles 


exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="novosort command failed exitcode=$exitcode  recalibration stopped for $sample"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
    exit $exitcode;
fi
if [ ! -s $realignedMergedFile ]
then
    MSG="novosort command produce no output file exitcode=$exitcode  recalibration stopped for $sample"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
    exit $exitcode;
fi



set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP3: Run  Recalibration                   ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

if [ $recalibrator == "BQSR" ]
then
    set +x; echo -e "\n\n" >&2;       
    echo "#################################################################################" >&2
    echo "                            recalibrator == BQSR" >&2
    echo "#################################################################################" >&2
    echo -e "\n\n" >&2; set -x;         


    set +x; echo -e "\n\n###### run BaseRecalibrator ########\n\n" >&2; set -x;
        
    if [ $input_type == "WES" -a ! -s $targetfile ]
    then
        set +x; echo -e "\n\n###### WES sample but no target file ########\n\n" >&2; set -x;    

        $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	--knownSites $refdir/$dbSNP \
	--knownSites $ndelfile \
	-I $realignedMergedFile \
	-T BaseRecalibrator \
	--out $recalreport \
	-nct $thr

    else
        set +x; echo -e "\n\n###### WES sample with target file ########\n\n" >&2; set -x;    

        $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	--knownSites $refdir/$dbSNP \
	--knownSites $ndelfile \
        -L $targetfile \
	-I $realignedMergedFile \
	-T BaseRecalibrator \
	--out $recalreport \
	-nct $thr        
    fi
    
    exitcode=$?
    echo `date`
    if [ $exitcode -ne 0 ]
    then
	MSG="BQSR BaseRecalibrator command failed exitcode=$exitcode  recalibration stopped"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit $exitcode;
    fi


    if [ ! -s $recalreport ]
    then
	MSG="$recalreport file not created, BQSR BaseRecalibrator failed to produce file. recalibration stopped for for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit 1;
    fi
    echo `date`

    set +x; echo -e "\n\n###### run PrintReads with -BQSR option ########\n\n" >&2; set -x;

    $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	-I $realignedMergedFile \
	-T PrintReads \
	-BQSR $recalreport \
	--out $recalibratedFile \
	-nct $thr

    exitcode=$?
    echo `date`
    if [ $exitcode -ne 0 ]
    then
	MSG="BQSR PrintReads command failed exitcode=$exitcode  recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit $exitcode;
    fi


    if [ ! -s $recalibratedFile ]
    then
	MSG="recal.$chr.real.recal.bam recalibrated file not created, BQSR PrintReads failed to produce file. recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit 1;
    fi
else
    set +x; echo -e "\n\n" >&2;       
    echo "#################################################################################" >&2
    echo "                            recalibrator IS NOT BQSR" >&2
    echo "#################################################################################" >&2
    echo -e "\n\n" >&2; set -x;         

    set +x; echo -e "\n\n###### run CountCovariates ########\n\n" >&2; set -x;
        
    if [ $input_type == "WES" -a ! -s $targetfile ]
    then
        set +x; echo -e "\n\n###### WES sample but no target file ########\n\n" >&2; set -x;    

        $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	--knownSites $refdir/$dbSNP \
	--knownSites $ndelfile \
	-I $realignedMergedFile \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-recalFile $covariatsFile 
    else
        set +x; echo -e "\n\n###### WES sample with target file ########\n\n" >&2; set -x;    

        $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	--knownSites $refdir/$dbSNP \
	--knownSites $ndelfile \
        -L $targetfile \
	-I $realignedMergedFile \
	-T CountCovariates \
	-cov ReadGroupCovariate \
	-cov QualityScoreCovariate \
	-cov CycleCovariate \
	-cov DinucCovariate \
	-recalFile $covariatsFile 

    fi

    exitcode=$?
    echo `date`		
    if [ $exitcode -ne 0 ]
    then
	MSG="cCountCovariates command failed exitcode=$exitcode  recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit $exitcode;
    fi

    if [ ! -s $covariatsFile ]
    then
	MSG="$covariatsFile CountCovariates file not created, recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit 1;
    fi
    echo `date`

    set +x; echo -e "\n\n###### run TableRecalibration        ########\n\n" >&2; set -x;


    $javadir/java -Xmx8g -Xms1024m  -Djava.io.tmpdir=$realrecaldir -jar $gatk/GenomeAnalysisTK.jar \
	-R $refdir/$ref \
	-I $realignedMergedFile \
	-T TableRecalibration \
	--out $recalibratedFile \
	-recalFile $covariatsFile 

    exitcode=$?
    if [ $exitcode -ne 0 ]
    then
	MSG="tablerecalibration command failed exitcode=$exitcode  recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit $exitcode;
    fi
    echo `date`		


    if [ ! -s $recalibratedFile ]
    then
	MSG="$recalibratedFile  file not created recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit 1;
    fi

fi # end choosing between BQSR and CountCovariates/TableRecalibration

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  STEP4: Run  QC                              ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;

if [ $runverify != "YES" ]
then
    set +x; echo -e "\n\n###### SKIP QC: samtools-calmd and VerifyBamId steps   ########\n\n" >&2; set -x;

    mv $recalibratedFile $outputfile
else
    set +x; echo -e "\n\n###### PERFORM QC STEPS: first samtools calmd and then VerifyBamID   ########\n\n" >&2; set -x;
   
   
    set +x; echo -e "\n\n###### run  samtools calmd                         ########\n\n" >&2; set -x;

    $samdir/samtools calmd -Erb $recalibratedFile $refdir/$ref > $outputfile
    exitcode=$?
    echo `date`
    if [ $exitcode -ne 0 ]
    then
	MSG="samtools calmd command failed exitcode=$exitcode  recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit $exitcode;
    fi
    if [ ! -s $outputfile ]
    then
	MSG="samtools calmd command failed to produce output file exitcode=$exitcode  recalibration stopped for $sample"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
	exit $exitcode;
    fi

    set +x; echo -e "\n\n####### indexing the bam file to make VerifyBamID happy#######\n\n" >&2; set -x;


    $sambambadir/sambamba index -t $thr $outputfile
    exitcode=$?
    echo `date`		
    if [ $exitcode -ne 0 ]
    then
        MSG="sambamba index command failed with outputfile=$outputfile exitcode=$exitcode for $sample"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
        exit $exitcode;
    fi

    $samdir/samtools flagstat $outputfile > $outputfile.flagstat
    exitcode=$?
    echo `date`		
    if [ $exitcode -ne 0 ]
    then
        MSG="samtools flagstat command failed with outputfile=$outputfile exitcode=$exitcode for $sample"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
        exit $exitcode;
    fi

    set +x; echo -e "\n\n####### run VerifyBamID                                         #######\n\n" >&2; set -x;

    $verifydir/verifyBamID --vcf $omnisites --bam $outputfile --out $verifiedOut --verbose --ignoreRB 
    exitcode=$?
    `echo date`

    if [ $exitcode -ne 0 ]
    then
        MSG="verifyBamID command failed with outputfile=$outputfile exitcode=$exitcode for $sample"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
        exit $exitcode;
    fi
    
    if [ ! -s ${verifiedOut}* ]
    then
        MSG="verifyBamID command failed with outputfile=$outputfile exitcode=$exitcode for $sample"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #>>$RealignOutputLogs/FAILED_recalibrate.${sample}.Anisimov.msg
        exit $exitcode;
    fi
    
    set +x; echo -e "\n\n####### generate QC report with  results of file *selfSM        #######\n\n" >&2; set -x;

    verifyfile=`find ./ -name "*.selfSM"`
    freemix=$( cat $verifyfile | sed '2q;d' | cut -f 7 ) 

    set +x; echo -e "\n\n" >&2; 
    echo -e "##########################################################################" >&2
    echo -e "checking that we actually grabbed something that is a number" >&2
    echo -e "using a weird trick to do the comparison in bash between noninteger variables" >&2
    echo -e "##########################################################################" >&2
    echo -e "\n\n" >&2; set -x;    
 
    if [ $(echo "$freemix < $freemix_cutoff"|bc) -eq 1 ]
    then
        set +x; echo -e "\n\n" >&2;    
	echo -e "##########################################################################" >&2
        echo -e "sample-$sample PASSED verifyBamID filter" >&2
	echo -e "##########################################################################" >&2
        echo -e "\n\n" >&2; set -x;	
	
        filterflag="PASSED"	
        detailfilter=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
        echo -e "$detailfilter" >> $qc_result
    else
    
        set +x; echo -e "\n\n" >&2;  
	echo -e "##########################################################################" >&2
        echo -e "sample-$sample FAILED verifyBamID filter" >&2
	echo -e "##########################################################################" >&2
        echo -e "\n\n" >&2; set -x;	

        filterflag="FAILED"
        detailfilter=$( echo -e "$sample\t$filterflag\tfreemix_value=$freemix\tfreemix_cutoff=$freemix_cutoff" )
        echo -e "$detailfilter" >> $qc_result
    fi
    echo `date`
    # end else verifybamid
fi

set +x; echo -e "\n\n" >&2;
echo "#################################################################################" >&2
echo "################  DONE. Exiting now                           ###################" >&2
echo "#################################################################################" >&2
echo -e "\n\n" >&2; set -x;
   