#! /bin/bash
################################################################################################ 
	#################### Torque preparation: PBS commands ###################
################################################################################################

#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=23			
#PBS -o localhost:$HOME/outputs.log.txt
#PBS -e localhost:$HOME/errors.log.txt		
#PBS -M aeahmed@illinois.edu
#PBS -m ae


echo -e "\n\n########################################################################################"
echo -e "#############               Summarizing the data :)              ###############"
echo -e "########################################################################################\n\n"


results="/home/groups/hpcbio_shared/azza/GIAB/results"
cd $results
echo parameter value Total_reads Total_aligned Mean_MAPQ Time > changing_parameters.txt

module load samtools/1.3.1
ls $results/p2_alignment > list
readarray parameters < list
rm list

for par in "${parameters[@]}"; do
	cd $results/p2_alignment/$par
	bams="a.$par.*.bam"
	bams=$(echo $bams | tr -d ' ')
	ls $bams|sed 's/.bam//g'|sed 's/a.//g'|cut -d'.' -f2- >range # need to account for cases when range is float
	readarray range < range
	i=0
	while [ $i -lt ${#range[@]} ] ; do
		file="a.$par.${range[i]}.summary.txt"
		summaryfile=$(echo $file | tr -d ' ')	#remove the space introduced from the array variable
		total=$(grep total $summaryfile | cut -d ' ' -f1)
		aligned=$(grep BWA $summaryfile|cut -d':' -f2)
		time=$(grep Execution $summaryfile|cut -d':' -f2)
		file="a.$par.${range[i]}.bam"
		bamfile=$(echo $file | tr -d ' ')
		mapq=$(samtools view $bamfile | awk '{sum+=$5} END { print sum/NR}')
		echo $par ${range[i]} $total $aligned $mapq $time >> $results/changing_parameters.txt
	        let i+=1
	done	
done


echo -e "\n\n########################################################################################"
echo -e "#############               Successful run :)              ###############"
echo -e "########################################################################################\n\n"
