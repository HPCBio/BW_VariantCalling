#!/bin/bash

## This script sweeps the parameters for the aligner BWA MEM, and allows comaprison between the quality of mapping in each case. It should be called as: bwa_sweep.sh <runfile> <SampleName> <read1> <read2>

scriptfile=$0
runfile=$1
SampleName=$2
read1=$3
read2=$4
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 4 ]
then
	MSG="Parameter mismatch.\nRerun like this: $0 <runfile>\n"
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
	exit 1;
fi

if [ `expr ${#R1}` -lt 1 ]
then
	    MSG="$R1 read one file not found"
	        echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
		    exit 1
	    elif [ ! -s $R1 ]
	    then
		        MSG="$R1 read one file not found"
			    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
			        exit 1    
			fi


echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

set -x
echo `date`
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
bwa_index=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
bwamemdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )

results=$outputdir/bwa_sweep.$SampleName

######### Alignment: The default settings
cd $results
mkdir default
cd default

START=$(date +%s)
$bwamemdir/bwa mem -M -t 12 -R  '@RG\tID:H7.5.R1\tPL:illumina\tPU:H7LH3CCXX.5.0\tLB:R1\tPI:0\tDT:2016-7-1\tSM:NA12878-Garvan' $refdir/human $reads/read1.fastq $reads/read2.fastq > a.default.0.sam
# note that it is the indexed file(s) that was needed in this stage. It has a different name than the 	refdir (remember the -p argument), so I need to use its name (human)
END=$(date +%s)
[ -s a.default.0.sam ] && echo "Default alignment successeful!" || exit
alignments=$($samtoolsdir/samtools view -c a.default.0.sam)
if [ "$alignments" -eq 0]; then
	echo			
	echo " Unfortunately, I can NOT process the your request with default parameters" 
	echo
	exit
fi
$samtoolsdir/samtools view -bS a.default.0.sam > a.default.0.bam

DIFF=$(( $END - $START ))

echo 
echo "BWA Mem aligned :$alignments: using parameter :default: (*=0=)"  > a.default.0.summary.txt
echo 
echo "Execution time is :$DIFF: seconds" >> a.default.0.summary.txt
echo 
$samtoolsdir/samtools flagstat a.default.0.bam >> a.default.0.summary.txt  # Generating summary statistics


######### Alignment: The combinatorial settings: changing a variable at a time, with the objective of having a sense of how things work


#echo -e "\n\n########################################################################################"
#echo -e "#############                CHECKING PARAMETERS                         ###############"
#echo -e "########################################################################################\n\n"

declare -a parameters=(k r w d c D m W A B O E L U T)
declare -a min=(3 .5 20 20 300 .1 20 0 1 1 1 1 1 1 10)
declare -a step=(3 .5 20 20 300 .1 20 3 2 2 2 2 2 3 10)
declare -a max=(60 4 200 200 10000 1 200 30 20 20 20 20 20 40 80)

cd $results/p2_alignment/
mkdir ${parameters[@]}

echo The parameters being tested and their ranges are given below:
echo paramters: ${parameters[@]}, 
echo minimum  : ${min[@]} 
echo maximum  : ${max[@]}

pos=0
while [ $pos -lt ${#parameters[@]} ]; do
        par=${parameters[pos]}
	cd $results/p2_alignment/$par
        for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
		START=$(date +%s)	
		$bwamemdir/bwa mem -t 12 -$par "$i" -M -R  '@RG\tID:H7.5.R1\tPL:illumina\tPU:H7LH3CCXX.5.0\tLB:R1\tPI:0\tDT:2016-7-1\tSM:NA12878-Garvan'  $refdir/human $reads/read1.fastq $reads/read2.fastq  > "a.$par.$i.sam"
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		if [ -s "a.$par.$i.sam" ]; then
			echo "Alignment successeful! with -$par $i" 
		else 
			echo "BWA aligned :0: using default parameter (*=0=)"> "a.$par.$i.summary.txt"
			echo "Execution time is :0: seconds" >> "a.$par.$i.summary.txt"
			continue
		fi
		alignments=$($samtoolsdir/samtools view -c a.$par.$i.sam)
		if [ "$alignments" -eq 0]; then
			echo			
			echo " Unfortunately, I can NOT process the parameter $par = $i with bwa mem" > "a.$par.$i.summary.txt"
			echo
			echo "BWA aligned :0: using default parameter (*=0=)" > "a.$par.$i.summary.txt"
			echo "Execution time is :0: seconds" >> "a.$par.$i.summary.txt"
			continue
		fi

		echo 
		echo "BWA Mem aligned :$alignments: using parameter :$par: =$i" > "a.$par.$i.summary.txt"
		echo 
		echo "Execution time is :$DIFF: seconds" >> "a.$par.$i.summary.txt"
		echo 
		$samtoolsdir/samtools view -bS "a.$par.$i.sam" >"a.$par.$i.bam"
		$samtoolsdir/samtools flagstat a.$par.$i.bam >> "a.$par.$i.summary.txt"
	done
	let pos+=1
done


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
