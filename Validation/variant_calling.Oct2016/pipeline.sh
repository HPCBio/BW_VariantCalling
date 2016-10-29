#! /bin/bash

# This script runs the pipeline as is (in other words, it invokes the start.sh script with proper input, and stores the log in the config folder)
# Note that nohup is used to record everything, and also be able to log out of the system safely without interrupting the workflow

var=$(./autoincrement.sh)
sed -i  '/^OUTPUT/s/.$/'"$var"'/' sample.runfile


nohup /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/Workflow/BW_VariantCalling/start.sh /home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/src/variant_calling.Oct2016/sample.runfile > /home/groups/hpcbio_shared/azza/human_WES_multi_sample_synthetic/src/variant_calling.Oct2016/sample_variant_calling.nohup 
