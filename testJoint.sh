#! /bin/bash

#PBS -S /bin/bash                       
#PBS -q default                         
#PBS -l nodes=1:ppn=23                  
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt          
#PBS -M aeahmed@illinois.edu
#PBS -m abe

runfile=/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/config/H3A_NextGen_assessment.Chr1_50X.set3.runfile

scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )

sample=H3A_NextGen_assessment_set3

outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )

TopOutputLogs=$outputdir/logs

sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )

deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )


rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )

echo $(find ${rootdir} -name "*.GATKCombineGVCF.raw.vcf" | sed "s/^/ --variant /g" | tr "\n" " ")


$scriptdir/joint_vcf.sh $runfile  $TopOutputLogs/log.mergeVcf.$sample.in $TopOutputLogs/log.merge.$sample.ou $TopOutputLogs/qsub.jointcall 


