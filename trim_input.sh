#!/bin/sh

if [ $# != 1 ]
then
     echo -e "program $0 stopped at line $LINENO\nREASON=Parameters mismatch"
     exit 1;
else
    set -x
    echo `date`
    email=grendon@illinois.edu
    infiles=$1
    trimMod=trimmomatic/0.33
    fastqcMod=fastqc/0.11.4
    rootdir=/home/groups/hpcbio/projects/LifeSpan/exome-March2016
    nodes=1
    threads=12
    queue=default
    mem=4gb
    outdir=$rootdir/results
    tmpdir=$rootdir/src
    adapters=/home/apps/trimmomatic/trimmomatic-0.33/adapters/TruSeq3-PE-2.fa
    
    echo -e "\n\n############ parameter check                  #################\n\n"

    if [ ! -s $infile ]
    then
	echo -e "$0 stopped at line $LINENO\nREASON=input file not found $infile"
	exit 1;
    fi
    if [ ! -s $adapters ]
    then
	echo -e "$0 stopped at line $LINENO\nREASON=input file not found $adapters"
	exit 1;
    fi
    if [ ! -d $rootdir ]
    then
        echo -e "$0 stopped at line $LINENO\sREASON=$rootdir directory not found"
        exit 1;
    fi
    if [ ! -d $outdir ]
    then
	echo -e "creating output folder $outdir";
	`mkdir $outdir`
    fi
    if [ ! -d $tmpdir ]
    then
	echo -e "creating output folder $tmpdir";
        `mkdir $tmpdir`
    fi


    echo -e "\n\n############ parameters ok                  #################\n\n"
    
    echo -e "\n\n############ trimming loop start here       #################\n\n"

    while read sampledetail
    do
       echo -e "\n############# processing next line in file...#################\n"

       if [ `expr ${#sampledetail}` -lt 1 ]
       then
           echo -e "\n############ skipping empty line            #################\n"
       else
           echo -e "\n############ processing: $sampledetails     #################\n"

           # we are assuming input is paired reads, three fields per row as follows

           samplename=$( echo $sampledetail | cut -d ' ' -f1 )
           R1=$( echo $sampledetail | cut -d ' ' -f2 )
           R2=$( echo $sampledetail | cut -d ' ' -f3 )
          
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

           echo -e "\n############ preparing output folders      #################\n"
           outputdir=$outdir/${samplename}/trimmed
           fqdir2=$outdir/${samplename}/FastQC-trimmed
           fqdir1=$outdir/${samplename}/FastQC-raw
           b1=${samplename}_R1           
           b2=${samplename}_R2


           mkdir -p $fqdir1
#           mkdir -p $fqdir2
           mkdir -p $outputdir          
          
           echo -e "\n############ results will go here: outputdir=$outdir/${samplename}    #################\n"           

 

	   qsub1=$tmpdir/qsub.trim.$samplename
	   echo "#PBS -S /bin/bash" > $qsub1
	   echo "#PBS -V" >> $qsub1
	   echo "#PBS -N trim.${samplename}" >> $qsub1
	   echo "#PBS -M $email" >> $qsub1
	   echo "#PBS -m ae" >> $qsub1
	   echo "#PBS -e $tmpdir/qsub.trim.${samplename}.er" >> $qsub1
	   echo "#PBS -o $tmpdir/qsub.trim.${samplename}.ou" >> $qsub1
	   echo "#PBS -l nodes=$nodes:ppn=$threads" >> $qsub1
	   echo "#PBS -l mem=$mem" >> $qsub1
	   echo "#PBS -q $queue" >> $qsub1
           echo "set -x" >> $qsub1
           echo "echo step1 fastqc on raw reads $R1 $R2" >> $qsub1           
#           echo "module load $fastqcMod"  >> $qsub1
#           echo "fastqc -o $fqdir1 -t $threads $R1" >> $qsub1
#           echo "fastqc -o $fqdir1 -t $threads $R2" >> $qsub1
#           echo "module unload $fastqcMod"  >> $qsub1
#           echo "module unload java/1.8.0_65"  >> $qsub1
           echo "echo `date`" >>  $qsub1           
           echo "echo step2 trim raw reads" >> $qsub1           
#           echo "module load $trimMod"  >> $qsub1
#	   echo "java -classpath /home/apps/trimmomatic/trimmomatic-0.33/trimmomatic-0.33.jar \
#   org.usadellab.trimmomatic.TrimmomaticPE \
#   -threads $threads \
#   -trimlog $outputdir/${samplename}_trim.log \
#   $R1 $R2 \
#   $outputdir/${b1}.paired.fq.gz $outputdir/${b1}.unpaired.fq.gz \
#   $outputdir/${b2}.paired.fq.gz $outputdir/${b2}.unpaired.fq.gz \
#   ILLUMINACLIP:${adapters}:2:20:10 LEADING:5 TRAILING:5 MINLEN:25 " >> $qsub1
#           echo "module unload $trim"  >> $qsub1
           echo "ln -s $R1 $outputdir/${b1}.paired.fq" >>  $qsub1
           echo "ln -s $R2 $outputdir/${b2}.paired.fq" >>  $qsub1
           echo "echo `date`" >>  $qsub1
           echo "echo step3 fastqc on trimmed reads" >> $qsub1           
           echo "ln -s $fqdir1 $fqdir2" >>  $qsub1
#           echo "module load $fastqcMod"  >> $qsub1
#           echo "fastqc -o $fqdir2 -t $threads $outputdir/${b1}.paired.fq.gz" >> $qsub1
#           echo "fastqc -o $fqdir2 -t $threads $outputdir/${b2}.paired.fq.gz" >> $qsub1
#           echo "module unload $fastqcMod"  >> $qsub1
           echo "echo `date`" >>  $qsub1
           echo "echo exiting now" >> $qsub1           
	   `chmod g+r $qsub1 `
	   alnjobid=`qsub $qsub1`
	   echo `date`
    fi           
    done < $infiles

fi




