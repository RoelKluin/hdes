#!/usr/bin/env perl

use warnings;
use strict;

my $bn = shift;

my %h;
my $end = 0;
my $ch = "new\n";
$ch .= "genome ${bn}.fa\n";
$ch .= "load ${bn}.bwa.bam\n" if -f "${bn}.bwa.bam";
$ch .= "load ${bn}.uqct_prev.bam\n" if -f "${bn}.uqct_prev.bam";
$ch .= "load ${bn}.uqct.bam\nsnapshotDirectory snaps\n";
print $ch;
my $start;
my $i = 0;

while (<>) {
 next unless /^[+-][^+-]/;
 chomp;
 my @e = split /\t/;
 my $p = $e[3];
 if (($p > $end) || ($ch ne $e[2])) {
     if (defined $start) {
         print "goto $ch:$start-$end\nsort base\ncollapse\nsnapshot ".$i++.".jpg\n";
     }
     $ch = $e[2];
     $start = $e[3];
     $end = $start + 51;
 } else {
     $end = $e[3] + 51;
 }
}
die "No changes\n" if not defined $start;
print "goto $ch:$start-$end\nsnapshot ".$i++.".jpg\n";
