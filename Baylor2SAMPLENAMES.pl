#!/usr/bin/perl
# first argument: full path to the output folder of the workflow
# second argument: full path to, and name of the information sheet provided by Baylor
# third argument: output file name where sample names will be recorded
# fourth argument: output filename where the processed information about samples and lanes etc will be recorded in a way that the workflow can ingest
# fifth argument: output filename where all this information is grouped across lanes, to facilitate merging of data per flowcell.


use strict;


if ($#ARGV < 4) { 
    print "parameter mismatch\nTo run type this command:\nperl $0 dirname infile outfile1 outfile2 outfile3\n\n"; 
    print "where:\ndirname is the directory where the infile is located\n";
    print "infile is a Baylor-style info sheet of the input files\n";
    print "outfile1 is usually set to SAMPLENAMES.list\n";
    print "outfile2 is usually set to SAMPLENAMES_multiplexed.list\n";
    print "outfile3 is usually set to SAMPLEGROUPS.list\n";
    exit 1; 
}

my $dirname=shift;
my $infile=shift;
my $outfile1=shift;
my $outfile2=shift;
my $outfile3=shift;
my $date=localtime();
print "starting $0 with $dirname $infile $outfile1 $outfile2 $outfile3\n$date\n";

chdir $dirname;
open(IN,"<$infile") || die("cannot open Baylor info sheet file $infile\n");
open(RAW,">$outfile1") || die("cannot open $outfile1\n");
open(MERGED,">$outfile2") || die("cannot open $outfile2\n");
open(GROUPS,">$outfile3") || die("cannot open $outfile3\n");

my $prevsampleid="";
my %samplexlane=();

# infile is a tab delimited text file with these columns
# col1: sample_id
# col2: flowcell
# col3: lane number
# col4: library type
# col5: library name
# last column: filename of read file
#
# we will construct a hash : {sample id} {name_of_flowcell_and_lane} {library, read1 or read2}

my %sampledet=();
my $counter=0;
while (my $line = <IN>) {
    #skip empty lines
    if ($line !~ /^\s*$/ ) {
        #skipt header
	if ($counter == 0) {
	    $counter++;
	} 
        else {
            # removing nasty hidden codes left by excel
            chomp $line;
	    $line =~ s/\r//g;
	    $line =~ s/\n//g;
	    $counter++;
            # extracting sample details, one column in the input file per array element
	    my @det = split(/\t/,$line);
            
            # parsing filename of read, assumed to be the last element of the array
            my $read=$det[$#det];
            # remove whitespace 
            $read =~ s/\s*//g;


            # build flowcell information
            my $flowcell="";
            $flowcell=$det[1];
            $flowcell =~ s/\s*//g;
            $samplexlane{$flowcell}=1;

            #build library information
            $sampledet{$det[0]}{$flowcell}{LB}=$det[4];
           

            # determine if this is a left read or a right read
            if ( $read =~ /(.*)_S\d_(.*)_R1_(.*)\.gz/ ) {
               $sampledet{$det[0]}{$flowcell}{R1}=$read;
            }
            else {
               if ( $read =~ /(.*)_S\d_(.*)_R2_(.*)\.gz/ ) {
                  $sampledet{$det[0]}{$flowcell}{R2}=$read;
               }  
               else {
	          print "failed to parse line $counter\tread=$read\n";
		  exit 1;
               }
	    }
	}
    }
}
print "\nDone parsing input\t$counter lines read.\nProducing $outfile1: A list of input files.\nNOTE: Only one line for paired reads\n";
$counter=0;
foreach my $key (sort keys %samplexlane) {  print RAW "$key\n"; $counter++; }
close(RAW);
print "done writing $outfile1\t$counter lines written.\n\n\nProducing $outfile2: A tab file with these columns\n";


$counter=0;
my $groupcounter=0;
print "SAMPLEID\tLANE_FLOWCELL_NAME\tLIBRARY\tREAD1\tREAD2\n";
# sort hash by sample id
foreach my $SID (sort keys %sampledet) {
    my $grouplanes="";
    # sort hash by flowcell-lane-ientifier
    foreach my $FL (sort keys %{ $sampledet{$SID} }) {
           # group all fastq that were run in all flowcells across all lanes, belonging to the same sample - so we could merge them all per sample after realrecal
           $grouplanes = $FL." ".$grouplanes;
           print MERGED "$SID\t$FL\t$sampledet{$SID}{$FL}{LB}\t$sampledet{$SID}{$FL}{R1}\t$sampledet{$SID}{$FL}{R2}\n";
           #$counter++;
    }
    print GROUPS "$SID\t$grouplanes\n";
    $groupcounter++;
}
close(MERGED);
close(GROUPS);
$date=localtime();
print "done writing $outfile2\t$counter lines written.\n\n\n";
print "Producing $outfile3: A tab file with these columns\n";
print "SAMPLEID\tlist of lanes separated by <SPACE>\n";
print "done writing $outfile3\t$groupcounter lines written.\n$date\n\n";
exit 0;
