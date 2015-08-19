#!/bin/bash
######################################
#  script to convert bam files back to fastq as pre requisite to alignment
#
######################################
##redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 8 ];
then
	MSG= "parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine""
        exit 1;
else					
	set -x
	echo `date`
        scriptfile=$0
        inputdir=$1
        newbamfiles=$2
        runfile=$3
        samplefileinfo=$4
        elog=$5
        olog=$6
        email=$7
        qsubfile=$8
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )

        if [ ! -d $inputdir ]
        then
	    MSG="INPUTDIR=$inputdir directory not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -s $samplefileinfo ]
        then
	    MSG="SAMPLEFILENAMES=$samplefileinfo file not found"
            echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi


        # finding BAMs in input folder and updating config files

        cd $inputdir
        newsamples=""
        sep=" "
        newnames="\n"
        bam_samples=""
        bamfiles=$( echo $newbamfiles | tr ":" " " )

        for BAM in $bamfiles 
        do
            if [ -s $BAM ]
            then
                newsample=$( echo $BAM | sed 's/.bam//' )
		newnames="BAM:${newsample}=${inputdir}/${BAM}\n$newnames"
                newsamples=${newsample}${sep}${newsamples}
            else
		MSG="$inputdir/$BAM bam2newbam file not found. updateconfig stopped."
		echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        done

        # updating BOTH config files
        if [ `expr ${#newnames}` -gt 0 ]
        then
            directory=`dirname $samplefileinfo`
            oldfile=`basename $samplefileinfo`
            oldfilename=$directory/$oldfile.old
            newfilename=$directory/$oldfile
            mv $samplefileinfo $oldfilename
            echo -e $newnames >> $newfilename
            if [ ! -s $newfilename ]
            then
		MSG="Empty file samplefilenames after bam2fastq conversion. updateconfig failed"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        else
	    MSG="Empty file samplefilenames after bam2fastq conversion. updateconfig failed"
            echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        # updating runfile
        directory=`dirname $runfile`
        oldfile=`basename $runfile`
        oldrun=$directory/$oldfile.old
        newrun=$directory/$oldfile
        mv $runfile $oldrun
        `perl $scriptdir/lned.pl $oldrun $newrun SAMPLENAMES "$newsamples"`
        exitcode=$?
        if [ $exitcode -ne 0 ]
        then
 	    MSG="lned.pl failed.  exitcode=$exitcode. update config files after bam2fastq conversion failed"
	    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit $exitcode;
        else
            if [ ! -s $newrun ]
            then
		MSG="Empty configuration files after bam2fastq conversion"
		echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
		exit 1;
            fi
        fi
	echo `date`
fi
