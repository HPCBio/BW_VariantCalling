#!/bin/bash 
#	
#  script to perform variant calling with HaplotypeCaller ONLY
#  This module is called from within the realign module
######################################

DEBUG=0

redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 10 ];
then
	MSG="parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
else					
	set -x
	echo `date`
        ulimit -s unlimited
        umask 0037
	scriptfile=$0
        outputdir=$1
        inputdir=$2
        inputfile=$3
        chr=$4
        region=$5
        runfile=$6
	elog=$7
	olog=$8
	email=$9
        qsubfile=${10}
	LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


       echo "##############################################################################"
       echo "##############################################################################"
       echo "##########         sanity checks with input params                 ###########"
       echo "##############################################################################"
       echo "##############################################################################"
       
        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
        region=$( echo $region | sed 's/:knownSites:/ /' | tr ":" " " )
        refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
        ref=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
        picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
        samdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
        javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
        gatk=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
        ped=$( cat $runfile | grep -w PEDIGREE | cut -d '=' -f2 )
        allsites=$( cat $runfile | grep -w EMIT_ALL_SITES | cut -d '=' -f2 )
        snvcaller=$( cat $runfile | grep -w SNV_CALLER | cut -d '=' -f2 )
        snvmixdir=$( cat $runfile | grep -w SNVMIXDIR | cut -d '=' -f2 )
        snvmixparms=$( cat $runfile | grep -w SNVMIX2PARMS | cut -d '=' -f2 )
        snvmixfilter=$( cat $runfile | grep -w SNVMIX2FILTER | cut -d '=' -f2 )
        uparms=$( cat $runfile | grep -w UNIFIEDGENOTYPERPARMS | cut -d '=' -f2 )
        onlyontarget=$( cat $runfile | grep -w TARGETTED | cut -d '=' -f2 )
        #javamodule=$( cat $runfile | grep -w JAVAMODULE | cut -d '=' -f2 )
        skipvcall=$( cat $runfile | grep -w SKIPVCALL | cut -d '=' -f2 )
        memprof=$( cat $runfile | grep -w MEMPROFCOMMAND | cut -d '=' -f2 )
        dbsnp=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
        genderinfo=$( cat $runfile | grep -w GENDERINFORMATION | cut -d '=' -f2 )

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
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$email""
            #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$email""
            exit 1;
        else
            if [ $skipvcall == "1" -o $skipvcall == "YES" ]
            then
		echo "skipping the execution of this variant calling module"
		exit 0;
            fi
        fi

        if [ `expr ${#region}` -lt 1 ]
        then
	    MSG="$region interval for vcall was not specified"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $inputdir ]
        then
	    MSG="$inputdir directory with realigned bams not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s ${inputdir}/${inputfile} ]
        then
	    MSG="${inputdir}/${inputfile} realigned bam file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $genderinfo ]
        then
	    MSG="$genderinfo gender file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $outputdir ]
        then
	    mkdir -p $outputdir
        fi

        if [ ! -d $picardir ]
        then
	    MSG="$picardir picard directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $samdir ]
        then
	    MSG="$samdir samtools directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ -z $javadir ]
        then
	    MSG="A value must be specified for JAVADIR in configuration file"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        #else
            #`/usr/local/modules-3.2.9.iforge/Modules/bin/modulecmd bash load $javamodule`
        #    `module load $javamodule`
        fi      
        if [ ! -d $gatk ]
        then
	    MSG="$gatk GATK directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi

        if [ ! -d $refdir ]
        then
	    MSG="$refdir reference genome directory not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi      
        if [ ! -s $refdir/$ref ]
        then
	    MSG="$ref reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $refdir/$dbsnp ]
        then
	    MSG="$refdir/$dbsnp dbSNP for reference genome not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ -z $snvcaller ]
        then
	    MSG="$snvcaller snvcaller tool was not specified in configuration file"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
	echo `date`
	
       echo "##############################################################################"
       echo "##############################################################################"
       echo "##########      setting up filters and parameters for snvcaller    ###########"
       echo "##############################################################################"
       echo "##############################################################################"

       #infile=`basename $inputfile`
       infile=$inputfile
       cd $outputdir

       if [ $snvcaller == "GATK" ]
       then
	   echo "snvcaller is GATK"
           if [[ $allsites == "YES" && $input_type == "exome" ]]
           then
               pedfile=$infile.$chr.raw.all.pbt.vcf
	       outfile=$infile.$chr.raw.all.vcf # This line can be removed? The outputfile naming is specified a bit further down below based on the site on which calling is done.
	       umode="EMIT_ALL_SITES"
	       utype="BOTH"
           else
               pedfile=$infile.$chr.raw.pbt.vcf
	       outfile=$infile.$chr.raw.g.vcf.gz
	       umode="EMIT_VARIANTS_ONLY"
	       utype="BOTH"
           fi
       elif [ $snvcaller == "SNVMIX" -o $snvcaller == "SNVMix" ]
       then
	   echo "snvcaller is SNVMIX"
           if [ $allsites == "YES" -a $input_type == "exome" ]
           then
	       snvfile=$infile.$chr.raw.snv.all.vcf
	       outfile=$infile.$chr.raw.indel.all.vcf
	       combfile=$infile.$chr.raw.multi.vcf
	       combparms="-V $outputdir/$snvfile -V $outputdir/$outfile"
	       umode="EMIT_ALL_SITES"
	       utype="INDEL"
               smode="all"
           else
	       snvfile=$infile.$chr.raw.snv.vcf
	       outfile=$infile.$chr.raw.indel.vcf
	       combfile=$infile.$chr.raw.multi.vcf
	       combparms="$-V outputdir/$snvfile -V $outputdir/$outfile"
	       umode="EMIT_VARIANTS_ONLY"
	       utype="INDEL"
	       smode="target"
           fi
       elif [ $snvcaller == "BEAUTY_EXOME" ]
       then
	   echo "snvcaller is BEAUTY_EXOME"
	   snvfile=$infile.$chr.raw.snvmix.vcf
	   outfile=$infile.$chr.raw.gatk.vcf
	   combfile=$infile.$chr.raw.multi.vcf
	   combparms="-V:GATK $outputdir/$outfile -V:SNVMix $outputdir/$snvfile -priority GATK,SNVMix" 
	   umode="EMIT_VARIANTS_ONLY"
	   utype="BOTH"
	   smode="target"
       else
	   MSG="SNV_CALLER=$snvcaller. This case is not currently available at this site"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	   exit 1;
       fi


       echo "############################################################################################"
       echo "############################################################################################"
       echo "##########      calculating variant calling w HaplotypeCaller                    ###########"
       echo "##########We did not check, but we assume that variant caller is HaplotypeCaller ###########"
       echo "############################################################################################"        

       echo `date`
	
       echo $genderinfo # Can take this line out later 
       strip_right="${outputdir%/*}"
       sample_id="${strip_right##*/}"
       gender=`cat $genderinfo | grep -v "#" | grep $sample_id | awk -F$'\t' '{print $2}'` # Can take this line out later 
       echo "Sample: $sample_id, gender: $gender" # Can take this line out later                                                                 
       echo $chr # Can take this line out later

       # Start of calling conditions. 

       ## Calling X 
       if [ $chr == "X" ];
       then
         echo "call X" # Can take this line out later
	 # The X PAR coordinates can be found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt
         x_par1="X:60001-2699520"
         x_par2="X:154931044-155260560"

         ### Calling male X
         if [ $gender == "Male" ];
         then
           echo "Calling male X" # Can take this line out later
           #### Calling male X_PAR1
           ploidy=2
           site="-L "$x_par1
           site_name="X_PAR1"
           outfile=$infile.$site_name.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi           
	   cp $outfile $deliverydir
	   echo `date`	   
           #### Calling male X_PAR2
           site="-L "$x_par2
           site_name="X_PAR2"
           outfile=$infile.$site_name.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi  
           cp $outfile $deliverydir
	   echo `date`
	   
           #### Calling male X_nonPAR
           ploidy=1
           site="-L X -XL "$x_par1" -XL "$x_par2
           site_name="X_nonPAR"
           outfile=$infile.$site_name.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi  
	   cp $outfile $deliverydir
	   echo `date`	   
         fi
         ### End calling male X  
	 echo `date`
	 
         ### Calling female X
         if [ $gender == "Female" ];
         then
           echo "Calling female X" # Can take this line out later
           ploidy=2
           outfile=$infile.$chr.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi  
	   cp $outfile $deliverydir
	   echo `date`	   
         fi
         ### End calling female X

       ## End calling X
       echo `date`
       
       ## Calling Y
       elif [ $chr == "Y" ];
       then
         echo "call Y"; # Can take this line out later
         # The X PAR coordinates can be found here: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000001405.25.regions.txt
         y_par1="Y:10001-2649520"
         y_par2="Y:59034050-59363566"

         ### Calling male Y
         if [ $gender == "Male" ];
         then
           echo "Calling male Y"  # Can take this line out later
           #### Calling male Y_PAR1
           ploidy=2
           site="-L "$y_par1
           site_name="Y_PAR1"
           outfile=$infile.$site_name.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi  
	   cp $outfile $deliverydir
	   echo `date`
	   
           #### Calling male Y_PAR2
           site="-L "$y_par2
           site_name="Y_PAR2"
           outfile=$infile.$site_name.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi  
	   cp $outfile $deliverydir
	   echo `date`
	   
           #### Calling male Y_nonPAR
           ploidy=1
           site="-L Y -XL "$y_par1" -XL "$y_par2
           site_name="Y_nonPAR"
           outfile=$infile.$site_name.raw.g.vcf.gz

           cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
           fi  
	   cp $outfile $deliverydir
	   echo `date`	   
         fi
         ### End calling male Y

       ## End calling Y

       ## Calling MT
       elif [ $chr == "MT" ];
       then
         echo "call MT" # Can take this line out later
         ploidy=1
         site_name=$chr
         outfile=$infile.$site_name.raw.g.vcf.gz

         cmd="$javadir/java -Xmx4g -Xms1g -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
         fi  
	 cp $outfile $deliverydir
	 echo `date`	 
       ## End calling MT

       ## Calling the autosomes
       else
         echo "Call 1->22" # Can take this line out later
         ploidy=2
         site_name=$chr

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
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
         fi  
	 cp $outfile $deliverydir
	 echo `date`	 
       fi
       ## End calling the autosomes

       echo "############################################################################################"
       echo "############################################################################################"
       echo "####################     End of calling conditions         #################################"
       echo "############################################################################################"
 

        if [ $ped != "NA" -a $snvcaller == "GATK" ]
        then
	    echo "##############################################################################"
	    echo "##############################################################################"
	    echo "##########      calculating phasebytransmission                    ###########"
	    echo "##############################################################################"
	    echo "##############################################################################"        
	    echo `date`

            $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    -v $outfile \
	    -T PhaseByTransmission \
            --ped $ped \
	    --out $pedfile

            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="phasebytransmission command failed exitcode=$exitcode. vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

            if [ ! -s $pedfile ]
            then
		MSG="$pedfile phasebytransmission file not created. vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi


        elif [ $snvcaller != "GATK" ]
        then
	    echo "##############################################################################"
	    echo "##############################################################################"
	    echo "########## calculating variants by snvmix and then merging results ###########"
	    echo "##############################################################################"
	    echo "##############################################################################"
	    echo `date`

            pilefile=$outfile.pileup
            tmpfile=$infile.$chr.tmp.snv
            $memprof $samdir/samtools mpileup -f $refdir/$ref $inputdir/$infile > $pilefile 

            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="samtools mpileup command failed exitcode=$exitcode . vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

	    if [ ! -s $pilefile ]
            then
		MSG="$pilefile pileup file not created. snvmix failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi

            if [ $smode == "all" ]
            then
		## question: modefile Mu_pi.txt does not exist in that folder
		$memprof $snvmixdir/SNVMix2 -i $pilefile -f -m $snvmixdir/Mu_pi.txt -o $tmpfile $snvmixparms
	    else
		$memprof $snvmixdir/SNVMix2 -i $pilefile -m $snvmixdir/Mu_pi.txt -o $tmpfile $snvmixparms
	    fi

            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="snvmix2 command failed.  exitcode=$exitcode vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

	    echo `date`

	    if [ ! -s $tmpfile ]
            then
		MSG="$tmpfile snvmix file not created. snvmix failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi

	    $memprof perl $scriptdir/snvmix_to_vcf.pl -i $tmpfile -o $snvfile
            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="snvmix2 command failed.  exitcode=$exitcode vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

	    echo `date`

	    if [ ! -s $snvfile ]
            then
		MSG="snv to vcf conversion failed"
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
	    fi

            echo "combining VCF files"

            $memprof $javadir/java -Xmx8g -Xms1024m -Djava.io.tmpdir=$outputdir -jar $gatk/GenomeAnalysisTK.jar \
	    -R $refdir/$ref \
	    $combparms \
	    -T CombineVariants \
	    -o  $combfile

            exitcode=$?
	    echo `date`
            if [ $exitcode -ne 0 ]
            then
		MSG="combinevariants command failed.  exitcode=$exitcode vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit $exitcode;
            fi

	    echo `date`

            if [ ! -s $combfile ]
            then
		MSG="$combfile combineVariants file not created. vcall failed."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		#echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        else
	    echo "##############################################################################"
	    echo "##############################################################################"
	    echo "##########      Skipping PhaseByTransmission and SNVMix            ###########"
	    echo "##############################################################################"
	    echo "##############################################################################"
	fi
	

fi
