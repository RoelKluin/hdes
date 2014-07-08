#!/usr/bin/perl
# (c) Roel Kluin, 2012, NKI
# Usage `samtools view $bam | nacolor | less -R'

use warnings;
use strict;

my %h;
my @seq = qw/T A C G U N/;
my $X = "\x1b[38;5;";

$h{$seq[$_]} = "\x1b[".(31+$_).";1m".$seq[$_] for 0..$#seq;

while (<>) {
	$_ = <>;
	chomp;
	print join('', map { exists $h{$_} ? $h{$_} : $_ } split(//))."\x1b[0m\n";
}

