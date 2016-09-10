#!/bin/sh

# This script checks the quality of reads in a folder, and trims known adapter sequences from reads using trimmomatic. It should be called as:
# trim_input.sh <runfile>

runfile=$1

if [ $# != 1 ]
then
     echo -e "program $0 stopped at line=$LINENO. \nREASON=Parameters mismatch"
     exit 1;
else
    set -x
    echo `date`
    rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
    tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
    email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
    sampleinfo=$( cat $runfile | grep SAMPLEINFORMATION | cut -d '=' -f2)
    adapters=$( cat $runfile | grep ADAPTERS | cut -d '=' -f2 )
    fastqcdir=$( cat $runfile | grep -w FASTQCDIR | cut -d '=' -f2)
    trimmomaticdir=$( cat $runfile | grep -w TRIMMOMATICDIR | cut -d '=' -f2)
    TopOutputLogs=$rootdir/logs
    thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
    nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
    queue=$( cat $runfile | grep -w PBSQUEUE | cut -d '=' -f2 )
    pbswalltime=$( cat $runfile | grep -w PBSWALLTIME | cut -d '=' -f2 )

    echo -e "\n\n############ parameter check                  #################\n\n"
    set -x 
    if [ ! -s $sampleinfo ]
    then
	echo -e "$0 stopped at line $LINENO. \nREASON=input file not found $sampleinfo"
	exit 1;
    fi
    if [ ! -s $adapters ]
    then
	echo -e "$0 stopped at line $LINENO. \nREASON=input file not found $adapters"
	exit 1;
    fi
    if [ ! -d $rootdir ]
    then
        echo -e "$0 stopped at line $LINENO. \sREASON=$rootdir directory not found"
        exit 1;
    fi
    if [ ! -d $TopOutputLogs ]
    then
	echo -e "creating Logs folder $TopOutputLogs";
	`mkdir $TopOutputLogs`
    fi
    if [ ! -d $tmpdir ]
    then
	echo -e "creating tmp folder $tmpdir";
        `mkdir $tmpdir`
    fi

    set +x
    echo -e "\n\n############ parameters ok                  #################\n\n"
    set -x 
    echo -e "\n\n############ trimming loop start here       #################\n\n"

    while read sampleLine
    do
       set +x	
       echo -e "\n############# processing next line in file...#################\n"
       set -x 
       if [ `expr ${#sampleLine}` -lt 1 ]
       then
	   set +x
           echo -e "\n############ skipping empty line            #################\n"
       else
           echo -e "\n############ processing: $sampleLine     #################\n"
           set -x 
           # we are assuming input is paired reads, three fields per row as follows

           samplename=$( echo "$sampleLine" | cut -d ' ' -f1 )
           R1=$( echo "$sampleLine" | cut -d ' ' -f2 )
           R2=$( echo "$sampleLine" | cut -d ' ' -f3 )
          
	   if [ `expr ${#samplename}` -lt 1 ]
 	   then
	       echo -e "$0 stopped at line $LINENO\nREASON=samplename string not found $samplename"
	       exit 1;
           fi

	   if [ ! -s $R1 ]
	   then
	       echo -e "$0 stopped at line $LINENO\nREASON=reads file not found $R1"
	       exit 1;
           fi

	   if [ ! -s $R2 ]
	   then
	       echo -e "$0 stopped at line $LINENO\nREASON=reads file not found $R2"
	       exit 1;
           fi
	   set +x
           echo -e "\n############ preparing output folders      #################\n"
	   set -x 
           outputdir=$rootdir/${samplename}/trimmed
           fqdir2=$rootdir/${samplename}/FastQC-trimmed
           fqdir1=$rootdir/${samplename}/FastQC-raw
           b1=${samplename}_R1           
           b2=${samplename}_R2


           mkdir -p $fqdir1
           mkdir -p $fqdir2
           mkdir -p $outputdir          
           set +x 
           echo -e "\n############ results will go here: outputdir=$outputdir    #################\n"           
	   set -x 
 

	   qsub1=$TopOutputLogs/qsub.trim.$samplename
	   echo "#PBS -S /bin/bash" > $qsub1
	   echo "#PBS -V" >> $qsub1
	   echo "#PBS -N trim.${samplename}" >> $qsub1
	   echo "#PBS -M $email" >> $qsub1
	   echo "#PBS -m ae" >> $qsub1
	   echo "#PBS -e $tmpdir/qsub.trim.${samplename}.er" >> $qsub1
	   echo "#PBS -o $tmpdir/qsub.trim.${samplename}.ou" >> $qsub1
	   echo "#PBS -l nodes=$nodes:ppn=$thr" >> $qsub1
	   echo "#PBS -l walltime=${pbswalltime}" >> $qsub1
	   echo "#PBS -l mem=$mem" >> $qsub1
	   echo "#PBS -q $queue" >> $qsub1
           echo "set -x" >> $qsub1
           echo "echo step1 fastqc on raw reads $R1 $R2" >> $qsub1           

           echo "$fastqcdir -o $fqdir1 -t $thr $R1" >> $qsub1
           echo "$fastqcdir -o $fqdir1 -t $thr $R2" >> $qsub1
           echo "echo `date`" >>  $qsub1           
           echo "echo step2 trim raw reads" >> $qsub1           
	   echo "java -jar $trimmomaticdir/trimmomatic-0.36.jar PE\
		   -threads $thr \
		   -trimlog $outputdir/${samplename}_trim.log \
		   $R1 $R2 \
		   $outputdir/${b1}.paired.fq.gz $outputdir/${b1}.unpaired.fq.gz \
		   $outputdir/${b2}.paired.fq.gz $outputdir/${b2}.unpaired.fq.gz \
		   ILLUMINACLIP:${adapters}:2:20:10 LEADING:5 TRAILING:5 MINLEN:25 " >> $qsub1
           echo "echo `date`" >>  $qsub1
           echo "echo step3 fastqc on trimmed reads" >> $qsub1           
           echo "$fastqcdir -o $fqdir2 -t $thr $outputdir/${b1}.paired.fq.gz" >> $qsub1
           echo "$fastqcdir -o $fqdir2 -t $thr $outputdir/${b2}.paired.fq.gz" >> $qsub1
           echo "echo `date`" >>  $qsub1
           echo "echo exiting now" >> $qsub1           
	   echo " echo "${samplename} $outputdir/${b1}.paired.fq.gz $outputdir/${b2}.paired.fq.gz" > ${rootdir}/sample.information " >> $qsub1
	   `chmod g+r $qsub1 `
	   jobid=`qsub $qsub1`
	   echo `date`
	  
    fi           
    done < $sampleinfo
    mv ${rootdir}/sample.information ${sampleinfo}
fi

