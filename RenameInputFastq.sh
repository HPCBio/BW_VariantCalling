#!/bin/bash

#RenameInputFastq.sh
# This script is required because data coming from different facilities tends to have differently formatted filenames
# thus it is difficult to construct code which will work on all kinds of file names
# and we decided to adopt a convention:
# our workflow looks for _read?.fq or _read?.fastq or _read?.fq.gz or _read?.fastq.gz

# this script will rename files to gave _read?.fq or _read?.fq.gz

# first variable = input folder
# second variable = read name ending pattern; 
# third variable = left read (1) or right read (2) 
#for example: H7FCKCCXX-6-IDX706_S8_L006_R1_001_1_sequence.txt.gz has ending pattern _R1_001_1_sequence.txt.gz and $3=1
#for example: H7FCKCCXX-6-IDX706_S8_L006_R2_001_2_sequence.txt.gz has ending pattern _R2_001_2_sequence.txt.gz and $3=2
#


umask 0027

InputFolder=$1
EndingPattern=$2
ReadSide=$3


for fastq in ${InputFolder}/*${EndingPattern} 
do
   #find the file name base
   file_name_base=`basename $fastq $EndingPattern`

   # determine if this is fastq or fastq.gz
   if [[ ( $fastq == *.fq ) || ( $fastq == *.fastq ) ]] 
   then
      file_extension="fq"
   elif [[ ( $fastq == *.gz ) ]] 
   then
      file_extension="gz"
   fi

   OldName=`basename $fastq`
   NewName=${file_name_base}_read${ReadSide}.fq.${file_extension}
   echo "renaming $OldName into $NewName"
   mv $fastq ${InputFolder}/${NewName}

done

