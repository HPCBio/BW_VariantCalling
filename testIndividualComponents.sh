#! /bin/bash

runfile=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/config/H3A_NextGen_assessment.Chr1_50X.set3.runfile

scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )

sample=H3A_NextGen_assessment_set3

outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )

TopOutputLogs=$outputdir/logs

sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )

deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )


rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )

echo $(find ${rootdir} -name "*.GATKCombineGVCF.raw.vcf" | sed "s/^/ --variant /g" | tr "\n" " ")

#$scriptdir/merge_vcf_and_bam.sh $runfile $sample $TopOutputLogs/log.mergeVcf.$sample.in $TopOutputLogs/log.merge.$sample.ou $TopOutputLogs/qsub.merge.$sample

while read sampleLine
do
    if [ `expr ${#sampleLine}` -lt 1 ]
    then 
        echo -e "\n\n########################################################################################" >&2
         echo -e "##############                 skipping empty line        ##############################" >&2
         echo -e "########################################################################################\n\n" >&2
     else
	echo -e "\n\n########################################################################################" >&2
        echo -e "#####         Processing next line $sampleLine                                ##########" >&2
        echo -e "##### col1=sample_name col2=read1 col3=read2  including full paths            ##########" >&2
        echo -e "##### sample_name will be used for directory namas and in RG line of BAM files##########" >&2
        echo -e "########################################################################################\n\n" >&2
      
        sample=$( echo "$sampleLine" | cut -d ' ' -f 1 )

        if [ `expr ${#sample}` -lt 1 ]
            then
                MSG="unable to parse line $sampleLine"
                echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1
        fi
       
      echo $sample >> $outputdir/$deliverydir/docs/samplesnames.txt

    fi

done <  $sampleinfo

