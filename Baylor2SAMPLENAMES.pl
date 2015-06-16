#!/usr/bin/perl


if ($#ARGV < 4) { print "param mismatch\nTo run type this command:\nperl $0 dirname infile outfile1 outfile2 outfile3\n\n"; exit 1; }
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
my $prevlane="";
my %samplexlane=();
my %sampledet=();
my $counter=0;
while ($line = <IN>) {
    #skip empty lines
    if ($line !~ /^\s*$/ ) {
        #skipt header
	if ($counter == 0) {
	    $counter++;
	} else {
            chomp $line;
	    $line =~ s/\r\n//g;
	    $counter++;
	    my @det = split(/\t/,$line);
            
	    if ($det[0] ne "") {
		# parsing filename of read1
		my $read1=$det[$#det];
                my $laneinfo="";
		if ( $read1 =~ /(.*)_S0_(.*)_1.sequence.txt.gz/ ) {
                    $laneinfo=$1;
                    $laneinfo =~ s/\s*//g;
		    $samplexlane{$laneinfo}=1;
		}  else {
		    print "failed to parse line $counter\tread1=$read1\n";
		    exit 1;
		}
		$sampledet{$det[0]}{$det[2]}{LB}=$det[4];
		$sampledet{$det[0]}{$det[2]}{IN}=$laneinfo;
		$sampledet{$det[0]}{$det[2]}{R1}=$read1;
		$prevsampleid=$det[0];
		$prevlane=$det[2];
	    } else {
		# parsing filename of read2
		my $read2=$det[$#det];
                $read2 =~ s/\s*//g;
		if ( $read2 !~ /(.*)_S0_(.*)_2.sequence.txt.gz/ ) {
		    print "failed to parse line $counter\tread2=$read2\n";
		    exit 1;
		}	
		$sampledet{$prevsampleid}{$prevlane}{R2}=$read2;
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
$groupcounter=0;
print "SAMPLEID\tLANE_NUMBER\tLANE_NAME\tLIBRARY\tREAD1\tREAD2\n";
foreach my $SID (sort keys %sampledet) {
    my $grouplanes="";
    foreach my $LN (sort keys %{ $sampledet{$SID} }) {
        $grouplanes =$sampledet{$SID}{$LN}{IN}." ".$grouplanes;
        print MERGED "$SID\t$LN\t$sampledet{$SID}{$LN}{IN}\t$sampledet{$SID}{$LN}{LB}\t$sampledet{$SID}{$LN}{R1}\t$sampledet{$SID}{$LN}{R2}\n";
	$counter++;
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
