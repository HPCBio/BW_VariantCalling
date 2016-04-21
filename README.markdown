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

## Installation

git clone https://github.com/HPCBio/BW_VariantCalling.git

## Documentation

## Features

# For users

## Usage examples

# For developers

## Workflow file naming conventions

__Input files__
Convention on input fastq names: must be in form of samplename_read?.fq(.gz) or samplename_read?.fastq(.gz)

__Scripts__
* qsub.whatever = qsub scripts
* log.whatever.in = output log for that qsub
* log.whatever.ou = pbs output message log for that qsub
* pbsWHATEVER = list of jobids for that step in the workflow 
* there are a few jobs that delineate major blocks of computation; their log files are in CAPS:
log.CONFIGURE.in

## Error capture

echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"

## Auto-archiving

The script autoArchive.sh is specific to Blue Waters and NCSA Globus endpoints.
It also makes use of Dart archiving utility, which is not yet available as opensource. 


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
