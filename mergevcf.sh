#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=lmainzer@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 8 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" 
        #echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        vardir=$1
        filenames=$2
        tabixdir=$3
        vcftoolsdir=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -d $vardir ]
	then
            MSG="$vardir directory not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            exit 1;	    
        fi

        # first, let's check that vcf files were produced for all chr        
	cd $vardir 
        if [ `expr ${#filenames}` -lt 1 ]
        then
            MSG="empty list of files to search for in variant calling folder"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            exit 1;	    
        fi

        listfiles=$( echo $filenames | tr ":" " " )
        prefix=$( echo $RANDOM )
        outfile=${prefix}".OUTPUT.vcf.gz"
        vcf_files="";

        for bam in $listfiles
        do
           if [ `find ./ -name ${bam}*.vcf | wc -l` -lt 1 ]
           then
               MSG="vcf file for $bam is missing in variant calling folder"
	       echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	       #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
               exit 1;
           else
               name=`find ./ -name ${bam}*.vcf`
               `$tabixdir/bgzip -c $name > ${name}.gz`
               exitcode=$?
               echo `date`
               if [ $exitcode -ne 0 ]
               then
		   MSG="indexing vcf file $name failed"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		   exit 1;
               fi
               `$tabixdir/tabix -p vcf ${name}.gz`
               exitcode=$?
               echo `date`
               if [ $exitcode -ne 0 ]
               then
		   MSG="indexing vcf file $name failed"
		   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
		   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		   exit 1;
               fi
               vcf_files=${vcf_files}" ${name}.gz"
           fi
        done        

        # next sep, concatenate the files
        `$vcftoolsdir/vcf-concat $vcf_files > ${outfile}.tmp`
        exitcode=$?
        echo `date`
        if [ $exitcode -eq 0 ]
        then
            `mv ${outfile}.tmp $outfile`
        else
            MSG="vcf-concat command failed. OUTPUT vcf file was not created."
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" 
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
            exit 1;	    
        fi          
        echo `date`
fi
