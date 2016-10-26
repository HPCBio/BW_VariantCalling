Variant Calling Workflow for HPC systems
========================================

This workflow originated from collaboration with Mayo Bioinformatics core group, and has now evolved mostly for testing the scalability of variant calling on many hundreds of genomes on large *HPC* systems.

  * [Variant Calling Workflow for HPC systems](#variant-calling-workflow-for-hpc-systems)
  * [Basics](#basics)
    * [Installation](#installation)
    * [Documentation](#documentation)
    * [Features](#Features)
  * [For users](#for-users)
    * [Usage examples](#usage-examples)
  * [For developers](#for-developers)
    * [Workflow file naming conventions](#workflow-file-naming-conventions)
    * [Error capture](#error-capture)
    * [Auto-archiving](#auto-archiving)
    * [Scalability testing](#scalability-testing)
    * [To-Dos](#to-dos)
    * [Example runfiles for commonly used configurations](#example-runfiles-for-commonly-used-configurations)


# Basics

## Introduction

This pipeline implements the GATK's best practices for germline variant calling in Whole Genome and Whole Exome Next Generation Sequencing datasets (https://software.broadinstitute.org/gatk/best-practices/), given a single sample or a cohort of samples. In its latest version, 3.6,  the best practices include the following stages:

1. Mapping to the reference genome
2. Marking duplicates
3. Base recalibration (BQSR)
4. Variant calling –----- (processing done per sample)
5. Joint genotyping –----- (processing done for all samples together)

These stages are implemented in our pipeline, with an optional  “Indel Realignment” step (which was recommended in previous GATK best practices < 3.6).



## Installation and dependencies

git clone https://github.com/HPCBio/BW_VariantCalling.git


The pipeline implements the stages of Figure [?] and [?], while allowing different software tools at some of the stages depending on user's preference. These are as shown in table [?] below, and it is assumed that the users would specify the path to each of them in theirhis runfile as shown in section ??.


Quality control
Fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ )

Illumina reads trimming
Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic )
read trimming implemented in the param sweep branch, not in main branch
TO-DO: work it inot main branch

Alignment
Bwa mem (https://github.com/lh3/bwa ), 
Novoalign (http://novocraft.com/ )

Marking duplicates
Samblaster (https://github.com/GregoryFaust/samblaster ), Novosort ( http://novocraft.com/ ),  
Picard (https://broadinstitute.github.io/picard/ ), 

Indel realignment
GATK (https://software.broadinstitute.org/gatk/download/ )

Base recalibration
GATK (https://software.broadinstitute.org/gatk/download/ )

Calling variants
GATK  (Haplotypecaller: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php )

Joint calling of variants
GATK (Genotypegvcf: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeGVCFs.php )

Miscelleaneous
Samtools (http://samtools.github.io/ )


Additionally, for the purposes of monitoring, parallelization and optimization, the workflow uses memprof, Anisimov launcher, YesWorkflow and parfu.

DECIDE WHETHER THESE NEED TO BE INSTALLED OR COULD BE OPTIONAL.


## Documentation

With an optional additional stage of checking the quality of input data, the pipeline can also be run as: Alignment stage only, Complete variant calling with realignment and Complete variant calling without realignment depending on the user’s ANALYSIS setting.

Under the hood, this pipeline splits and merges files at different stages to achieve optimal usage of resources. 

## Features


# Input files
Convention on input fastq names: must be in form of samplename_read?.fq(.gz) or samplename_read?.fastq(.gz)

For this pipeline to work, a number of standard files for calling variants are needed, namely the reference sequence, database of known variants and the adapter sequence to be trimmed. The full path to all these needs to be specified in the User’s runfile as specified in section ?.?

It is important to note that the reference sequence should be prepared first, following the GATK’s guideline (http://gatkforums.broadinstitute.org/wdl/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk). 

For working with human data, one can download most of the needed files from the GATK’s resource bundle: http://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it . Missing from the bundle are the index files for the aligner, which are specific to the tool that would be used for alignment (i.e., bwa or novoalign in this pipeline)

To achieve the parallelization of Figure [2] in the realignment/recalibration stages, the pipeline needs a separate vcf file for each chromosome/contig, and each should be named as: \*\${chr_name}.vcf. If working with the GATK bundle, the sample script (splitVCF-by-chromosome.sh) can be used to produce the needed files with some minor modifications (mainly, providing the right path to the referencedir, java and GenomeAnalysisTK.jar)


## Scripts
* qsub.whatever = qsub scripts
* log.whatever.in = output log for that qsub
* log.whatever.ou = pbs output message log for that qsub
* pbsWHATEVER = list of jobids for that step in the workflow 
* there are a few jobs that delineate major blocks of computation; their log files are in CAPS:
log.CONFIGURE.in

## User’s runfile and sample information files

To run a specific stage of the pipeline, in addition to specifying the needed script file, the user needs to supply 2 additional files, these are the runfile and the sampleinfo files.
The sampleinformation file contains the information about the samples to be processed by the pipeline. In its current implementation, it can analyze paired end WES/WGS data in fastq/fq/fastq.gz/fq.gz format only. These should be specified in tabular form of 3 columns separated by ‘space’, and according to the format below:
<sample name> <full path to read1 file> <full path to read2 file>
The runfile file contains all the details regarding a specific run of the pipeline, including the tools of section 2.1, resources of section 2.2, and the sampleinformation file path as well. It would change depending on the analysis type required.

## Results

The results from a typical run of the pipeline are organized according to the hierarchy shown in Figure [3] below. Overall, the DELIVERYFOLDER contains the key summarizing files of the run (the cleaned up bams, gvcfs and final vcf from joint calling; in addition to the summary reports regarding the quality of the data, and copies of the sampleinformation and runfile files). Each sample also has its own directory that contains the files generated after each stage. In Figure [3], a color coding schema is employed to differentiate the files that would be generated according to how the user specifies the ANALYSIS parameter in the runfile.


## Error capture

echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"

## Auto-archiving

The script autoArchive.sh is specific to Blue Waters and NCSA Globus endpoints.
It also makes use of Dart archiving utility, which is not yet available as opensource. 

The userids are now explicitly stated in the script, but will be moved to config file in the due course of time.


## Scalability testing

The scripts
 realrecal_NONsharedRefGenomeMode.sh
    and
 realrecal_sharedRefGenomeMode.sh    
    and
 schedule_realrecal_sharedRefGenomeMode.sh
are not a part of the workflow, but rather exist to enable scalability testing on the shared reference files.


## To-Dos

1. discuss mailing to RedMine in automatic fashion, and in a portable fashion
instead of commenting out the emails when debugging
2. remove #PBS -V from all scripts
3. incorporate checking the results of FastQC in automatic fashion; edit the AssessFastQCresults.sh script
4. incorporate checking vcf automatically?
5. make sure the parsing of jobids is still correct after the Jan 18 of 2016 maintenance:
   " During the upcoming Jan 18 maintenance (starting at 06:00 AM) Blue Waters scheduler will be upgraded and migrated to a new server. The current server name (nid11293) will change to bwsched and as a result the full jobid will change from <JOBID>.nid11293 to <JOBID>.bwsched (viewed from the qstat output).  If you have any scripts that use the "nid11293" string please use the new server name after the upgrade. Everything else regarding job submission will remain the same".
6. edit the comments statements so they look good in the logs
7. brush through and remove any mention of specific names or institutions or projects (found Gerrit mentioned somewhere)

## Example runfiles for commonly used configurations

1. Bundle all jobs by chromosome, do maximum QC, use qsub directly without aprun  
   1.1. Align-Realign-VariantCalling  
    ..1.1.1. analyze individual samples, sequenced 1 per lane  
    ..1.1.2. analyze multiplexed samples  
    ..1.1.3. perform a multisample analysis, where each sample was sequenced without multiplexing   
   1.2. AlignOnly - each individual fastq or pair of fastq  
2. Unbundle all jobs, do minimum QC, use qsub directly without aprun  
   2.1. Align-Realign-VariantCalling: analyze individual samples, sequenced 1 per lane  
   2.2. AlignOnly: analyze individual samples, sequenced 1 per lane   
