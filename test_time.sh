#!/bin/sh
# written in collaboration with Mayo Bioinformatics core group
#tracejob -q -n 10 $963987 | grep -m 1 "resources_used.walltime=" | sed 's/^.*resources_used.walltime=//' |  awk 'BEGIN { FS=":"} { print $1*3600+$2*60+$3}'




cputime=`tracejob -q -n 10 963896 | grep -m 1 "resources_used.walltime=" | sed 's/^.*resources_used.walltime=//' |  awk 'BEGIN { FS=":"} { print $1*3600+$2*60+$3}'`
echo $cputime


