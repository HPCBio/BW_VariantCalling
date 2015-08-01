This workflow originated of collaboration with Mayo Bioinformatics core group, 
and has now evolved mostly for testing the scalability of variant calling on many hundreds of genomes on large *HPC* systems.

# Workflow file naming conventions

* qsub.whatever = qsub scripts
* log.whatever.in = output log for that qsub
* log.whatever.ou = pbs output message log for that qsub
* pbsWHATEVER = list of jobids for that step in the workflow 

* there are a few jobs that delineate major blocks of computation; their log files are in CAPS:
log.CONFIGURE.in

## Error capture

echo -e "$MSG\n\nDetails:\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
