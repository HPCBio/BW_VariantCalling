#!/bin/bash

################################################################################################ 
        #################### Torque preparation: PBS commands ###################
################################################################################################

#PBS -S /bin/bash                       
#PBS -q default                         
#PBS -l nodes=1:ppn=20                  
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt          
#PBS -M aeahmed@illinois.edu
#PBS -m abe

module load R/3.2.3 
cd /home/groups/hpcbio_shared/azza/GIAB/src


Rscript bqsrplotting.R

