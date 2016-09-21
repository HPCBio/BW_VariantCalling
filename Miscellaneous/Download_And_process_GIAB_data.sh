#!/bin/bash

cd /home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw

wget -bqc https://s3-ap-southeast-2.amazonaws.com/kccg-x10-truseq-nano-v2.5-na12878/NA12878_V2.5_Robot_1_R1.fastq.gz

wget -bqc https://s3-ap-southeast-2.amazonaws.com/kccg-x10-truseq-nano-v2.5-na12878/NA12878_V2.5_Robot_1_R2.fastq.gz

wget -bqc https://s3-ap-southeast-2.amazonaws.com/kccg-x10-truseq-nano-v2.5-na12878/NA12878_V2.5_Robot_1.hc.vqsr.vep.vcf.gz

echo NA12878_V2.5_Robot_1 /home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw/NA12878_V2.5_Robot_1_R1.fastq.gz /home/groups/hpcbio_shared/azza/GIAB/reads/Gravan_raw/NA12878_V2.5_Robot_1_R2.fastq.gz >/home/groups/hpcbio_shared/azza/GIAB/config/GIAB_NA12878_30X.sampleinfo

at now + 5 hours -f /home/groups/hpcbio_shared/azza/GIAB/src/Invoking_variant_calling_pipeline.sh

at now + 72 hours -f /home/groups/hpcbio_shared/azza/GIAB/src/Validating_pipeline_variants.sh

