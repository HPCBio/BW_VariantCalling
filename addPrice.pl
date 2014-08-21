#!/usr/bin/perl
# written in collaboration with Mayo Bioinformatics core group

#  program: addPricing.pl
#  parses a Summary report produced by GGPS v2.2 and adds pricing info
#  Author: Gloria Rendon
#  Date: Aug 2013

#  getting input 
#use lib "/projects/mayo/builds/Date-Calc-6.3/lib/";
use lib "/projects/mayo/builds/Carp-Clan-6.04/lib";
#use Date::Calc qw(Delta_Days);
use Switch;

my $infile="";
my $outfile="";
my $rate="";
my $sucorrection=16/12;
my $redmine="hpcbio-redmine\@igb.illinois.edu";
#my $redmine="grendon\@illinois.edu";
my $email="";

if ($#ARGV<2) {
    print "Program that parses a Summary.Report file produced during a run of the GGPS v2.2 pipeline and adds SUs and COST information.\nIt produces a new updated report that goes in the same folder as that of the Summary.Report. It also send an email notification to the redmine system\n\n";
    print "EXAMPLE:$0 file=fullPathToSummaryReport rate=Rate email=someEmailAccount\n\nNote: the \@ in the email value must be typed as a  \\\@ like this  email=youremail\\\@yourdomain \n";
    exit 1;
}

foreach my $argval (@ARGV) {
    switch ($argval) {
        case /file/ { $infile = substr($argval,5); }
        case /rate/ { $rate = substr($argval,5); }
        case /email/ { $email = substr($argval,6); }
        else { die("$argval Unkown input parameter\n"); }
   }
}


die("A filename must be specified\n") if $infile eq "";
die("A Rate must be specified\n") if $rate eq "";
die("An SUcorrection must be specified\n") if $sucorrection eq "";


open(SUMMARY,"<$infile") || die("Cannot open $infile\n");
$outfile= $infile.".new.wpricing";
open(REPORT,">$outfile") || die("Cannot open $outfile\n");

#here goees the code for figuring out recency of report
#my $filedateflat="";
#my $numdays="";

#$filedateflat=`ls --full-time $infile | cut -d ' ' -f 6`;
#my @filedate=split(/-/,$filedateflat) if $filedateflat ne "";
#my (undef,undef,undef,$day,$month,$year)=localtime();
#my @today= ( $year+1900, $month+1, $day );
#$numdays=(Delta_Days(@filedate,@today))+2;
#$numdays=60 if $numdays eq "";
#print "file=$infile recency of file in days=$numdays applicable rate=$rate";
#sleep(5);

# here goes the code for reading and parsing SummaryReport

$/='This jobid:';  # set record delimiter
my $totprice=0;
my $totSUs=0;
my $MSG="REVISED SUMMARY REPORT\nEXPANDED  WITH ESTIMATED COSTS\n\n";
my $jobids="";
while (<SUMMARY>) {
    if ( /GGPS pipeline/ )  {
        # top portion
	$MSG = $MSG.$_;
        $MSG =~ s/This jobid://;
        $MSG = $MSG."\nDetails of each job:\n\n";
    }
    if ( /JOBIDS((.|\n)*)/ ) {
        # table portion
	$jobids=$1;
        $jobids =~ s/This jobid://;
        $jobids =~ s/\n\n/ /g;
        $jobids =~ s/  / /g;
    }
}

# here goes the code for gathering stats on the jobs

if ( $jobids ne "" ) {
   $qhist= `time qhist -f jobid,jobname,queue,status,nds,ppn,usedwall,su --with-tabs $jobids`;
   my $exitcode=$?;
   if ( $exitcode != 0 ) {
    $MSG=$MSG."\nqhist command failed.\nFile:  $infile \ncost calcutating program has been stopped.\n";
    `echo -e "$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email,$redmine""`;
    exit 1;
   }
   # the result here does not look like what I get from the command line
   # the last two TOTAL lines do not show up
   # so, we need to do post procession here
   @qhistlines=split(/\n/,$qhist);
   foreach my $qline (@qhistlines) {
       if ( $qline =~ /^JobId/ ) {
	   $MSG = $MSG.$qline."\n";
       }
       elsif ( $qline =~ /\t/) {
           @values=split(/\t/,$qline);
           $MSG=$MSG."@values\n";
           $SUs= pop @values;
	   $totSUs+=$SUs;
       }
   }
} else {
    $MSG=$MSG."\nparsing of Summary report failed.\nFile:  $infile \ncost calcutating program has been stopped.\n";
    `echo -e "$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email,$redmine""`;
    exit 1;
}

# totals go here and geneating reports

$totprice= $sucorrection * $totSUs * $rate;

$MSG=$MSG."\nTotal Service Units (SUs) for this run: $totSUs\nRate per Service Unit: \$ $rate\nCorrection per Service Unit: $sucorrection\nTotal Cost for this run:  \$ $totprice\n\n\nNote: The formula for calculating total price: Total SUs * correction * Rate per SU.\n\nThe correction is based on the priority used when launching the job:\nCorrection for Normal queue with default Priority is 16/12\nCorrection for Himem queue is 32/64\n\n\nDisclaimer: This is NOT a bill for the jobs included in this analysis. All statistics and estimated costs reported in this analysis are subject to revisions and updates.\n";
#print "$MSG";
print REPORT "$MSG";
close(SUMMARY);
close(REPORT);

# code for sending email to redmine goes here
`echo -e "$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$email,$redmine""`;
exit 0;
