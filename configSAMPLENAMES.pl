#!/usr/bin/perl
# first argument: full path to the output folder of the workflow
# second argument: full path to, and name of the information sheet provided by user
# third argument: output file name where sample names will be recorded
# fourth argument: output filename where the processed information about samples and lanes etc will be recorded in a way that the workflow can ingest
# fifth argument: output filename where all this information is grouped across lanes, to facilitate merging of data per flowcell.


use strict;


if ($#ARGV < 5) { 
    print "parameter mismatch\nTo run type this command:\nperl $0 dirname infile filetype outfile1 outfile2 outfile3\n\n"; 
    print "where:\ndirname is the directory where the infile is located\n";
    print "infile is a tab-separated infosheet of the input files\n";
    print "filetype is the file format of the input file, FASTQ or BAM\n";
    print "outfile1 is usually set to SAMPLENAMES.list\n";
    print "outfile2 is usually set to SAMPLENAMES_multiplexed.list\n";
    print "outfile3 is usually set to SAMPLEGROUPS.list\n";
    exit 1; 
}

my $dirname=shift;
my $infile=shift;
my $filetype=shift;
my $outfile1=shift;
my $outfile2=shift;
my $outfile3=shift;
my $date=localtime();
my $intype=uc $filetype;
print "starting $0 with $dirname $infile $intype $outfile1 $outfile2 $outfile3\n$date\n";

chdir $dirname;
open(IN,"<$infile") || die("cannot open Baylor info sheet file $infile\n");
open(RAW,">$outfile1") || die("cannot open $outfile1\n");
open(MERGED,">$outfile2") || die("cannot open $outfile2\n");
open(GROUPS,">$outfile3") || die("cannot open $outfile3\n");

if ( $intype ne "BAM" && $intype ne "FASTQ" ) { 
     print "$intype unknown value. It must be FASTQ or BAM\n"; exit 1 
}

my $prevsampleid="";
my %samplexlane=();

# infile is a tab delimited text file with 3 or 5 columns. 3 columns for non-multiplexed samples or 5 columns for multiplexed samples
# col1: sample_id
# col2: filename of read1
# col3: filename of read2
# col4: flowcell     (optional for multiplexed samples only)
# col5: library name (optional for multiplexed samples only)


my %sampledet=();
my $counter=0;
my $prevSize=0;
while (my $line = <IN>) {
    #skip empty lines
    if ($line !~ /^\s*$/ ) {
        #skipt header
	if ($counter == 0) {
            # initialize counters and skip header line
	    $counter++;
	    chomp $line;
	    my @det = split(/\t/,$line);
	    my $arraySize = @det;
	    $prevSize = $arraySize;
	} 
        else {
            # removing nasty hidden codes left by excel
            chomp $line;
	    $line =~ s/\r//g;
	    $line =~ s/\n//g;
	    $counter++;
            # extracting sample details, one column in the input file per array element
	    my @det = split(/\t/,$line);

            # sanity check; quit or skip incomplete record?  we'll quit for now 
            my $arraySize = @det;
            if ( $arraySize != $prevSize ) {
                        # infile has lines with different number of fields
	    	        print "Unable to process line=$line\ninfile has lines with different number of fields\n\n";
	    	        exit 1;
	    } else {
	    		$prevSize = $arraySize;
            }
            if ( $arraySize == 5 && $intype eq "FASTQ" ) {
		       # we have a multiplexed sample, paired reads, fastq format

		       # parsing the row begins
		       my $sampleid= $det[0];               
		       my $read1   = $det[1];
		       my $read2   = $det[2];
		       my $flowcell= $det[3];
		       my $library = $det[4];

		       # remove whitespace 
		       $sampleid =~ s/\s*//g;              
		       $read1 =~ s/\s*//g;
		       $read2 =~ s/\s*//g;
		       $flowcell =~ s/\s*//g;
		       $library =~ s/\s*//g;


		       # populating hashes

		       $samplexlane{$flowcell}=1;

		       $sampledet{$sampleid}{$flowcell}{LB}=$library;

		       $sampledet{$sampleid}{$flowcell}{source}=$read1."\t".$read2;

		  
	    } elsif (  $arraySize == 3 && $intype eq "FASTQ"  )  {
		       # we have a non-multiplexed sample, paired reads, fastq format
		       
		       # parsing the row begins
		       my $sampleid= $det[0];               
		       my $read1   = $det[1];
		       my $read2   = $det[2];

		       # remove whitespace 
		       $sampleid =~ s/\s*//g;              
		       $read1 =~ s/\s*//g;
		       $read2 =~ s/\s*//g;

		       # populating hashes

		       $sampledet{$sampleid}=$read1."\t".$read2;		       
		       
	    } elsif (  $arraySize == 2 && $intype eq "BAM"  )  {
	    
		       # we have a non-multiplexed sample, bam format
		       
		       # parsing the row begins
		       my $sampleid= $det[0];               
		       my $bamfile = $det[1];


		       # remove whitespace 
		       $sampleid =~ s/\s*//g;              
		       $bamfile =~ s/\s*//g;

		       # populating hashes

		       $sampledet{$sampleid}=$bamfile;
		       	    
	    } else {
	    	       # none of the above. exit with error
	    	       print "Unable to process line=$line\n";
	    	       exit 1;
	    }
	}
        # next line please
    }
}
print "\nDone parsing input\t$counter lines read.\nProducing $outfile1: A list of input files.\nNOTE: Only one line for paired reads\n";
$counter=0;
# lets find out which file we just processed, 5-column or 3-column or 2-column
if ( $prevSize == 5 ) {
        # infile has 5 columns, it is multiplexed,  two hashes were populated
        
        foreach my $key (sort keys %samplexlane) {  print RAW "$key\n"; $counter++; }
} else {
        # infile has fewer than 5 columns, it is non-multiplexed,  one hash was populated
        
        foreach my $key (sort keys %sampledet) {  print RAW "$key\n"; $counter++; }
}
close(RAW);
print "done writing $outfile1\t$counter lines written.\n\n\nProducing $outfile2 and $outfile3\n";


if ( $prevSize == 5 ) {
        # infile has 5 columns, it is multiplexed,  two hashes were populated

	print "SAMPLEID\tREAD1\tREAD2\tLANE_FLOWCELL_NAME\tLIBRARY\n";
	# sort hash by sample id
	foreach my $SID (sort keys %sampledet) {
	    my $grouplanes="";
	    # sort hash by flowcell-lane-ientifier
	    foreach my $FL (sort keys %{ $sampledet{$SID} }) {
		   # group all fastq that were run in all flowcells across all lanes, belonging to the same sample - so we could merge them all per sample after realrecal
		   $grouplanes = $FL." ".$grouplanes;
		   print MERGED "$SID\t$sampledet{$SID}{$FL}{source}\t$FL\t$sampledet{$SID}{$FL}{LB}\n";
		   #$counter++;
	    }
	    print GROUPS "$SID\t$grouplanes\n";
	}
} else  {
        # infile has fewer than 5 columns, it is non-multiplexed,  one hash was populated
	print "SAMPLEID\tREAD1\tREAD2\n" if $prevSize == 3;
	print "SAMPLEID\tBAMFILE\n" if $prevSize == 2;	
	# sort hash by sample id
	foreach my $SID (sort keys %sampledet) {
            print MERGED "$SID\t$sampledet{$SID}\n";
	    print GROUPS "$SID\t$SID\n";
	} 	
 	
}        
	
close(MERGED);
close(GROUPS);
$date=localtime();
print "done. Exiting now\n$date\n\n";
exit 0;
