#! /bin/bash

# This script runs the pipeline as is (in other words, it invokes the start.sh script with proper input, and stores the log in the config folder)
# Note that nohup is used to record everything, and also be able to log out of the system safely without interrupting the workflow

nohup /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/Workflow/BW_VariantCalling/start.sh /home/groups/hpcbio_shared/azza/GIAB/config/GIAB_NA12878_30X.runfile > /home/groups/hpcbio_shared/azza/GIAB/config/GIAB_NA12878_30X_variant_calling.nohup
