#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
if [ $# != 4 ]
then
        MSG="Parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"        
        exit 1;
fi
	set -x
	echo `date`	
        scriptfile=$0
        runfile=$1
        email=$2
        elog=$3
        olog=$4
        umask 0027

        outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
        deliveryfolder=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
        TopOutputLogs=$outputdir/logs
        pbsprj=$( cat $runfile | grep -w PBSPROJECTID | cut -d '=' -f2 )
        pbsqueue=$( cat $runfile | grep -w PBSQUEUEWGEN| cut -d '=' -f2 )
        pipeid=$( cat $TopOutputLogs/CONFIGUREpbs )        
        archivefolder=$outputdir/archive
        dartcmd=/projects/sciteam/jti/builds/dart_0_8_5
        skipmd5="YES"   #this could be changed later on
        
        echo "##############################################################################"
        echo "##############################################################################"
        echo "##########   sanity checks with input params and transfer end points  ########"
        echo "##############################################################################"
        echo "##############################################################################"

        if [ ! -s $runfile ]
        then
	    MSG="$runfile configuration file not found"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	   #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ `expr ${#outputdir}` -lt 2 ]
        then
	    MSG="$outputdir folder not specified in configuration file"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi
        if [ ! -d $outputdir ]
        then
  	    MSG="$outputdir folder not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	elif [ ! "$(ls -A $outputdir)" ]
        then
	    MSG="$outputdir Is empty"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	    
        fi

        if [ ! -d $TopOutputLogs ]
        then
  	    MSG="$TopOutputLogs folder not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	fi
	
        if [ `expr ${#deliveryfolder}` -lt 2 ]
        then
             delivery=$outputdir/delivery
        else
             delivery=$outputdir/$deliveryfolder
        fi
        if [ ! -d $delivery ]
        then
  	    MSG="$delivery folder not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	elif [ ! "$(ls -A $delivery)" ]
        then
	    MSG="$delivery folder is empty. Execution of pipeline did not produce results. Nothing to archive. Exiting now"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
	    
        fi
       
        # names of output folders for each sample should be specified in SAMPLENAMES.list. Checking that the file exists
        if [ ! -s $outputdir/SAMPLENAMES.list ]
        then
	    MSG="$outputdir/SAMPLENAMES.list file not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
	    exit 1;
        fi        

        # this is where the dart files will go
        if [ ! -d $archivefolder ]
        then
            echo -e "$archivefolder directory does not exist, creating it"
            mkdir -p $archivefolder
        else
            echo -e "$archivefolder directory does  exist, resetting it" 
            rm -r $archivefolder
            mkdir -p $archivefolder
        fi
        chmod -R 770 $archivefolder

        if [ ! -s $dartcmd ]
        then
     	    MSG="$dartcmd dart executable  not found"
   	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
   	    #echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""
   	    exit 1;
        fi
        

        echo "##########################################################################################"
	echo "##########                 declaring END POINTS                                   ########"
        echo "##########################################################################################"

         ###############################################################################################################
         #this is how end points are normally specified
         #Destination="ncsa#Nearline/projects/sciteam/jti/TestingAutoTransfer_baylor_barch_1" 
         #Source="ncsa#BlueWaters/scratch/sciteam/lmainzer/H3A_DiversityGenotyping_Space/DartPackages/Striped_Minus1"
         #we need to declare them here for this particular project/batch
         ###############################################################################################################

        thisbatch=`basename $outputdir`
        thisproject=/projects/sciteam/jti/Results
        #thisproject=/projects/sciteam/jti/TestingAutoTransfer_batch2
        source_endpoint="ncsa#BlueWaters"
        destination_endpoint="ncsa#Nearline"       
        source=${source_endpoint}${archivefolder}
        destination=${destination_endpoint}${thisproject}/$thisbatch
        #user=$( echo $email | cut -d '@' -f1 )


        echo "##########################################################################################"
	echo "########## TESTING END POINTS: 1) that they exist 2) that we have access to them  ########"
        echo "##########################################################################################"
        #We assume that the userid is identical on both end points AND on globusonline

        user=$( echo $USER )
        if [ `expr ${#user}` -lt 2 ]
        then
            MSG="$user value is missing." 
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""           
            exit 1
        elif [ $user == "rendong" -o $user == "grendon" ]
        then
             user="gloriarendon"     #because  userid on BW is different from  userid on GO
        fi
        
        echo "checking ssh cli.globusonline.org on $source" 
        ssh $user@cli.globusonline.org "ls $source" 

        if [ $? -ne 0 ]
        then
            MSG="ssh cli.globusonline.org command failed on SOURCE $source." 
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""           
            exit 1
        else
            echo -e "$source is ok"
        fi

        echo "checking ssh cli.globusonline.org on DESTINATION $destination" 
        ssh $user@cli.globusonline.org "ls $destination" 
        if [ $? -ne 0 ]
        then
            echo "$destination folder does not exist, creating it now" 
            ssh $user@cli.globusonline.org "mkdir $destination" 
        else
            echo -e "$destination is ok"
        fi
        
	echo `date`

        echo "##########################################################################################"
        echo "##########    Now lets create the README file and send a copy to archive          ########"
        echo "##########################################################################################"

        ListGOFiles=$( cat $outputdir/SAMPLENAMES.list | sed "s/$/\.drt\n/g" | tr " " "\n" )
        readme=$delivery/docs/README.txt
        chmod 770 $readme
        truncate -s 0 $readme
        thisDate=$( echo `date` )         
        echo -e "On $thisDate $user archived the sample folders \nfrom: $outputdir \nto: $destination\n\n" >> $readme
        echo -e "The command used to archive was \n$dartcmd -C -e 20 -B 1000 -F archivefolder/{SampleFolder}.drt -r outputdir/SampleFolder \n\n"   >> $readme
        echo -e "The command to extract one or more files is \n$dartcmd -X -F $destination/{SampleFolder}.drt  -r yourCurrentPath \n\n" >> $readme
        echo -e "Notice that the extract command MUST be issued from the folder where the file will be EXTARCTED TO\n\n"  >> $readme
        echo -e "List of files in $destination \n\n" >> $readme
        echo -e "$ListGOFiles"  >> $readme
        echo -e "\n\n\nNote: Each drt file is a packaged sample folder\n\n\n"   >> $readme

        cp $readme $archivefolder/README.txt

        ssh $user@cli.globusonline.org "transfer --label README --verify-checksum -- $source/README.txt $destination/README.txt "
        if [ $? -ne 0 ]
        then
            MSG="transfer of $readme to $destination failed" 
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" #| ssh iforge "mailx -s '[Support #200] variant identification pipeline' "$redmine,$email""           
            exit 1
        else
            echo -e "transfer of $readme to $destination was ok"
        fi
        
        echo `date`

	
        echo "##########################################################################################"
        echo "##########              Create al the qsub jobs  here for md5sum                  ########"
        echo "##########          They will be executed from within a launcher                  ########"
        echo "##########################################################################################"

        
        if [ $skipmd5 == "YES" ]
        then 
             echo "##########################################################################################"
	     echo "##########       SKIP md5sum STEP because GO does it for us                       ########"
             echo "##########################################################################################"
        else     
            echo "##########################################################################################"
            echo "##########       Let's create the qsubs for the md5sum launcher                   ########"
            echo "##########################################################################################"
        
            echo -e "resetting the joblists and qsub files ..."
        
            truncate -s 0 $TopOutputLogs/md5sumBAMS.AnisimovJoblist
            truncate -s 0 $TopOutputLogs/md5sumVCFS.AnisimovJoblist
            truncate -s 0 $TopOutputLogs/md5sum.pbs
 
            qsub_md5sumBAMS=$TopOutputLogs/qsub_md5sumBAMS.AnisimovJoblist
            qsub_md5sumVCFS=$TopOutputLogs/qsub_md5sumVCFS.AnisimovJoblist


            echo -e "walltimes for these jobs..."        
            md5sumBAMS_cputime="06:00:00"
            md5sumVCFS_cputime="02:00:00"  
        
            echo -e "counters for the anisimov launchers"
            md5sumBAMS_files=1
            md5sumVCFS_files=1

            echo -e "these are the lists of files that MAY need md5sum"
        
            bamfiles1=`find $delivery -name *.bam -type f`   
            vcfFiles=`find  $delivery -name *.vcf.gz -type f`
        
            echo -e "populating the anisimov launcher for the md5sum of BAMS"
            for bam in $bamfiles1
            do
                echo -e "nohup md5sum $bam > ${bam}.md5\n" >> $TopOutputLogs/md5sumBAMS.AnisimovJoblist
                (( md5sumBAMS_files++ ))
            done

            echo -e "populating the anisimov launcher for the md5sum of VCFS" 
            for vcf in $vcfFiles
            do
                echo -e "nohup md5sum $vcf > ${vcf}.md5\n" >> $TopOutputLogs/md5sumVCFS.AnisimovJoblist
                (( md5sumVCFS_files++ ))
            done
        
            echo -e "calculating number of nodes needed for each launcher"  

            NumberOfProcPerNode=16  #there are 32 cpus per node, each job uses 2 cpus
        
            if [ $md5sumBAMS_files -lt 17 ]
            then
                 md5sumBAMS_nodes=1
                 NumberOfProcPerNode=$md5sumBAMS_files
            else
                 md5sumBAMS_nodes=$(( md5sumBAMS_files/NumberOfProcPerNode ))
                 if [ `expr $md5sumBAMS_nodes % $NumberOfProcPerNode` -gt 0 ]
                 then
                     (( md5sumBAMS_nodes++ )) # there is a remainder in that division, and we give those remaining jobs an extra node
                 fi
            fi        

            if [ $md5sumVCFS_files -lt 17 ]
            then
                 md5sumVCFS_nodes=1
                 NumberOfProcPerNode=$md5sumVCFS_files
            else
                 md5sumVCFS_nodes=$(( md5sumVCFS_files/NumberOfProcPerNode ))
                 if [ `expr $md5sumVCFS_nodes % $NumberOfProcPerNode` -gt 0 ]
                 then
                     (( md5sumVCFS_nodes++ )) # there is a remainder in that division, and we give those remaining jobs an extra node
                 fi
            fi        
       
            # we need to add one more node for the launcher
        
            (( md5sumBAMS_nodes++ ))
            (( md5sumVCFS_nodes++ ))
        
             echo -e "NOW we have all the pieces together to form the qsub files for the TWO anisimov launchers"
        
             cat $outputdir/qsubGenericHeader > $qsub_md5sumBAMS
             echo "#PBS -N ${pipeid}_md5sumBAMS_Anisimov" >> $qsub_md5sumBAMS
             echo "#PBS -o $TopOutputLogs/log.md5sumBAMS_Anisimov.ou" >> $qsub_md5sumBAMS
             echo "#PBS -e $TopOutputLogs/log.md5sumBAMS_Anisimov.in" >> $qsub_md5sumBAMS
             echo "#PBS -l nodes=${md5sumBAMS_nodes}:ppn=32:xe" >> $qsub_md5sumBAMS       
             echo "#PBS -l walltime=$md5sumBAMS_cputime" >> $qsub_md5sumBAMS
             echo "aprun -n $md5sumBAMS_files -N $NumberOfProcPerNode -d 2 ~anisimov/scheduler/scheduler.x $TopOutputLogs/md5sumBAMS.AnisimovJoblist /bin/bash > $TopOutputLogs/log.md5sumBAMS_Anisimov.ou " >> $qsub_md5sumBAMS
             md5sumBAMS_job=`qsub $qsub_md5sumBAMS`  
             `qhold -h u $md5sumBAMS_job`
             echo $md5sumBAMS_job >> $TopOutputLogs/md5sumBAMS.AnisimovJoblist
             echo $md5sumBAMS_job >> $TopOutputLogs/md5sum.pbs
        
             cat $outputdir/qsubGenericHeader > $qsub_md5sumVCFS
             echo "#PBS -N ${pipeid}_md5sumVCFS_Anisimov" >> $qsub_md5sumVCFS
             echo "#PBS -o $TopOutputLogs/log.md5sumVCFS_Anisimov.ou" >> $qsub_md5sumVCFS
             echo "#PBS -e $TopOutputLogs/log.md5sumVCFS_Anisimov.in" >> $qsub_md5sumVCFS
             echo "#PBS -l nodes=${md5sumVCFS_nodes}:ppn=32:xe" >> $qsub_md5sumVCFS       
             echo "#PBS -l walltime=$md5sumVCFS_cputime" >> $qsub_md5sumVCFS
             echo "aprun -n $md5sumVCFS_files -N $NumberOfProcPerNode -d 2 ~anisimov/scheduler/scheduler.x $TopOutputLogs/md5sumVCFS.AnisimovJoblist /bin/bash > $TopOutputLogs/log.md5sumVCFS_Anisimov.ou " >> $qsub_md5sumVCFS
             md5sumVCFS_job=`qsub $qsub_md5sumVCFS`  
             `qhold -h u $md5sumVCF_job`
             echo $md5sumVCFS_job >> $TopOutputLogs/md5sumVCFS.AnisimovJoblist
             echo $md5sumVCFS_job >> $TopOutputLogs/md5sum.pbs
        fi

        
        echo "##########################################################################################"
        echo "######       prep work for the dart and transfer operations                            ###"
        echo "##########################################################################################"
        
        echo -e  "this variable will be empty if we skip the md5sum block. However, the code will still work"

        if [ $skipmd5 != "YES" ]
        then
             md5sum_all=$( cat $TopOutputLogs/md5sum.pbs | sed "s/\..*//" | tr "\n" ":" ) 
        fi

        echo -e "####################       resetting lists and resources                   #############"
        GO_blockSize=32                # number of GO jobs that can be packaged together
        GO_count=0                     # number of GO packages
        i=1                            # for keeping track of GO_blocks
        GO_cputime="01:00:00"
        DART_cputime="02:00:00"
        DART_nodes=10                  # this is the total PER job
        numGOfiles=$( cat $outputdir/SAMPLENAMES.list | wc -l )
        GO_joblist=$TopOutputLogs/GO_Joblist_${i}        
        truncate -s 0 $TopOutputLogs/DART2GO_Joblist_${i}.pbs        
        truncate -s 0 $GO_joblist
        truncate -s 0 $TopOutputLogs/DART.pbs
        truncate -s 0 $TopOutputLogs/GOTransfer.pbs
        echo `date`
        
        echo -e "##########################################################################################"
        echo -e "######  Next loop for creating qsub jobs for DART and populating lists of GO transfers  ##"
        echo -e "##########################################################################################"
        
        while read SampleFolder 
        do
           if [ `expr ${#SampleFolder}` -lt 1 ]
           then
               echo -e "##########################################################################################"           
               echo -e "###########                skipping empty line                                      ######"
               echo -e "##########################################################################################"               
           else
               echo -e "##########################################################################################"
               echo -e "#####           qsub job for DART packaging of $SampleFolder                        ######"
               echo -e "##########################################################################################"
               
               qsub_dart=$TopOutputLogs/qsub_DART_$SampleFolder
               cat $outputdir/qsubGenericHeader > $qsub_dart
               echo "#PBS -N ${pipeid}_DART_$SampleFolder" >> $qsub_dart
               echo "#PBS -o $TopOutputLogs/log.DART_$SampleFolder.ou" >> $qsub_dart
               echo "#PBS -e $TopOutputLogs/log.DART_$SampleFolder.er" >> $qsub_dart              
               echo "#PBS -l nodes=${DART_nodes}:ppn=32" >> $qsub_dart   
               echo "#PBS -l walltime=$DART_cputime" >> $qsub_dart
               if [ $skipmd5 != "YES" ]
               then
                    echo "#PBS -W depend=afterok:$md5sum_jobs " >>  $qsub_dart
               fi
               echo "set -x" >>  $qsub_dart
               echo "newgrp ILL_jti " >>  $qsub_dart
               echo "chmod -R g+rX $outputdir/$SampleFolder " >>  $qsub_dart
               echo "aprun -n 80 -d 4 $dartcmd -C -e 20 -B 1000 -F $archivefolder/${SampleFolder}.drt -r $outputdir/$SampleFolder " >> $qsub_dart
               echo "chmod g+rX $archivefolder/${SampleFolder}.drt " >>  $qsub_dart
               chmod 770 $qsub_dart
               dart_job=`qsub $qsub_dart`
               `qhold  -h u $dart_job`
               echo $dart_job >> $TopOutputLogs/DART.pbs

               echo -e "##########################################################################################"
               echo -e "##### adding job to anisimov launcher for GO archival of $SampleFolder              ######"
               echo -e "##########################################################################################"

               cmd="\"transfer --label $SampleFolder --verify-checksum --encryp -- $source/${SampleFolder}.drt $destination/$SampleFolder/${SampleFolder}.drt \""

               if [ $GO_count -lt $GO_blockSize ]
               then
                    echo -e "#####            GO tranfer for $SampleFolder goes in this package $GO_joblist ######"               
                    echo -e "ssh $user@cli.globusonline.org $cmd" >> $GO_joblist
                    (( GO_count++ ))
                    echo $dart_job >> $TopOutputLogs/DART2GO_Joblist_${i}.pbs                    
               else
                    (( i++ ))
                    GO_count=0
                    GO_joblist=$TopOutputLogs/GO_Joblist_${i}       
                    truncate -s 0 $GO_joblist
                    truncate -s 0 $TopOutputLogs/GO_Joblist_${i}.pbs
                    echo -e "ssh $user@cli.globusonline.org $cmd" >> $GO_joblist
                    echo $dart_job >> $TopOutputLogs/DART2GO_Joblist_${i}.pbs      
                    echo -e "#####            GO tranfer for $SampleFolder goes in this package $GO_joblist  ######"                   
               fi
           fi  
           echo -e "#####           NEXT SAMPLE PLEASE....        #####"
        done < $outputdir/SAMPLENAMES.list

        echo "##########################################################################################"
        echo "##########                  OUTSIDE  THE MAIN LOOP                                ########"
        echo "##########       Now generating the qsubs for all GO joblists                     ########"
        echo "##########################################################################################"
        echo `date`     
        numGOLists=$i   #this is the number of GO job lists that were generated inside the loop. Each one has up to 32 jobs inside

        for inx in $numGOLists
        do
            GO_joblist=$TopOutputLogs/GO_Joblist_${inx}
            DART2GO_pbs=$( cat $TopOutputLogs/DART2GO_Joblist_${inx}.pbs | sed "s/\..*//" | tr "\n" ":" )
	    qsubGO=$TopOutputLogs/qsub_GO_Joblist_$inx
	    
            echo "##########################################################################################" 	    
            echo -e "#####            creating qsub for this GO package $GO_joblist                   ######"
            echo "##########################################################################################" 	    
            
	    cat $outputdir/qsubGenericHeader > $qsubGO
            echo "#PBS -N ${pipeid}_GO_Joblist_$inx" >> $qsubGO
            echo "#PBS -o $TopOutputLogs/log_GO_Joblist_${inx}.ou" >> $qsubGO
            echo "#PBS -e $TopOutputLogs/log_GO_Joblist_${inx}.in" >> $qsubGO
            echo "#PBS -l nodes=1:ppn=32:xe" >> $qsubGO
            echo "#PBS -l walltime=$GO_cputime" >> $qsubGO
            echo "#PBS -W depend=afterok:${DART2GO_pbs} "  >> $qsubGO
            echo "set -x " >> $qsubGO
            cat "$GO_joblist"   >> $qsubGO
            chmod 770 $qsubGO
            GO_job=`qsub $qsubGO`  
            `qhold -h u $GO_job`
            echo $GO_job >> $TopOutputLogs/GOTransfer.pbs
        done

        echo "##########################################################################################"
        echo "##########    Now we release all the jobs for execution                           ########"
        echo "##########################################################################################"
        md5sum_all=$( cat $TopOutputLogs/md5sum.pbs | sed "s/\..*//" | tr "\n" " " )
        globus_all=$( cat $TopOutputLogs/GOTransfer.pbs | sed "s/\..*//" | tr "\n" " " )
        dart_all=$( cat $TopOutputLogs/DART.pbs | sed "s/\..*//" | tr "\n" " " )


        `qrls -h u $md5sum_all`
        `qrls -h u $dart_all`   
        `qrls -h u $globus_all`        

        echo `date`
        echo -e "exiting now..."
