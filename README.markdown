This workflow originated of collaboration with Mayo Bioinformatics core group, 
and has now evolved mostly for testing the scalability of variant calling on many hundreds of genomes on large *HPC* systems.

# Workflow file naming conventions

*Input files*

Convention on input fastq names: must be in form of samplename_read?.fq(.gz) or samplename_read?.fastq(.gz)

*Scripts*

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
are not a part of the workflow, but rather exist to enable scalability testing on the shared reference files.
