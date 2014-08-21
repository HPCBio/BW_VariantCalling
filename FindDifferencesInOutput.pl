#!/usr/bin/perl
#
# FindDifferencesInOutput.pl
# written by LSY
# in the lab of Matt Hudson
# University of Illinois
# March 30, 2011
#
# segregate reads that are found in file 2, not file 1, and thpse found in both files.
# 

use strict;
use IO::File;

my $usage = "./FindDifferencesInOutput.pl *TopHits1 *TopHits2 *NotIn1 *InBoth1and2
                  \n*=required ";

#ARGS--------------------------------
die "Usage: $usage\n" unless ( @ARGV > 0 );


open(TopHits1,      '<', $ARGV[0]) || die("Could not open file!");
open(TopHits2,      '<', $ARGV[1]) || die("Could not open file!");
open(NotIn1,        '>', $ARGV[2]) || die("Could not open file!");
open(InBoth1and2,   '>', $ARGV[3]) || die("Could not open file!");





################ BEGIN CONSTRUCTING HASHES 1


# initialize the hash of reads names for end 1
my %reads_names1 = ();
my $count = 0;
print "COLLECTING reads names 1 \n";
while (<TopHits1>) {
   $count++;
   $_ =~ s/\n|\r//;

   my @line = split('\t', $_);
   $line[0] =~ s/>//; # read name, stripped of the > sign
   #print "$line[0]\n";

   #$reads_names1{ $line[0] } = 1; # just setting the flag to say the read is there
   $reads_names1{ $line[0] } = $_; 

   # to keep track of progress
   #unless ($count%10000) {
      #print "read 1 number $count \n";
   #}
}



################ print out to files


my $count = 0;
print "COLLECTING reads names 2 \n";
while (<TopHits2>) {
   $count++;
   #$_ =~ s/\n|\r//;

   my @line = split('\t', $_);
   #my @line = split('\s+', $_);



   $line[0] =~ s/>//; # read name, stripped of the > sign
   # in case we need to chop off the reading frame info
   #$line[0] =~ s/_orf_\d+_frame=(|-)\d//;

   #print "$line[0]\n";

   unless ( $line[0] eq 'read' ) {

      if ( exists $reads_names1{ $line[0] } ) {
         # if the read exists in the first file too, then it goes into InBoth1and2:
         #print InBoth1and2 "$line[0]\n";
         print InBoth1and2 "$reads_names1{ $line[0] }\t$_";
         #print InBoth1and2 "$_";
      }
      else {
         # otherwise it goes into NotIn1:
         #print NotIn1 "$line[0]\n";
         print NotIn1 "$_";
      }
   }

   # to keep track of progress
   #unless ($count%10000) {
      #print "read 1 number $count \n";
   #}
}

