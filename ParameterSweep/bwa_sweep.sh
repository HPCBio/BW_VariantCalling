#!/bin/bash
#PBS -S /bin/bash                       
#PBS -q default                         
#PBS -l nodes=1:ppn=23                  
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt          
#PBS -M aeahmed@illinois.edu
#PBS -m abe


## This script sweeps the parameters for the aligner BWA MEM, and allows comaprison between the quality of mapping in each case. It should be called as: bwa_sweep.sh <runfile> <SampleName> <read1> <read2>

scriptfile=$0
runfile=$1
SampleName=$2
read1=$3
read2=$4
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 4 ]; then
	MSG="Parameter mismatch.\nRerun like this: $0 <runfile>\n"
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
	exit 1;
fi

if [ `expr ${#read1}` -lt 1 ]; then
        MSG="$read1 read one file not found"
        echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        exit 1
elif [ ! -s $read1 ]; then
        MSG="$read1 read one file not found"
        echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        exit 1    
fi

if [ `expr ${#read2}` -lt 1 ]; then
    MSG="$read2 read two file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"        
    exit 1
elif [ ! -s $read2 ]; then
    MSG="$read2 read two  file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"        
    exit 1    
fi

if [ `expr ${#SampleName}` -lt 1 ] 
then
    MSG="$SampleName sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1
else
    sID=$SampleName
    sPU=$SampleName
    sSM=$SampleName
fi

if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ]
then
    MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${RGparms}"  | tr ":" "\t" | tr "=" ":" )

echo -e "\n\n########################################################################################"
echo -e "#############                Sweeping parameters starts here!              ###############"
echo -e "########################################################################################\n\n"

set -x
echo `date`
bwa_index=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
bwamemdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )

results=$outputdir/bwa_sweep.$SampleName

######### Alignment: The default settings
cd $results
mkdir default
cd default

START=$(date +%s)
$bwamemdir/bwa mem -M -t $thr -R "${rgheader}" $bwa_index $read1 $read2 > a.default.0.sam
END=$(date +%s)
[ -s a.default.0.sam ] && echo "Default alignment successeful!" || exit
alignments=$($samtoolsdir/samtools view -c a.default.0.sam)
if [ "$alignments" -eq 0]; then
	echo			
	echo " Unfortunately, I can NOT process the sample $SampleName with default parameters" 
	echo "Exiting now"
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

cd $results
mkdir ${parameters[@]}

echo The parameters being tested and their ranges are given below:
echo paramters: ${parameters[@]}, 
echo minimum  : ${min[@]} 
echo maximum  : ${max[@]}

pos=0
while [ $pos -lt ${#parameters[@]} ]; do
        par=${parameters[pos]}
	cd $results/$par
        for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
		START=$(date +%s)	
		$bwamemdir/bwa mem -t $thr -$par "$i" -M -R  "${rgheader}" $bwa_index $read1 $read2   > "a.$par.$i.sam"
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
