#!/usr/bin/env perl

use warnings;
use strict;

my %h;
my $end = 0;
my $ch = "new\nload hg19_GL.1.bwa.bam\nload hg19_GL.1.uqct_prev.bam\nload hg19_GL.1.uqct.bam\nsnapshotDirectory snaps\n";
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
     } else {
         print $ch;
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
