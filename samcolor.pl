#!/usr/bin/perl
# Usage `samtools view $bam | csam | less -R'

use warnings;
use strict;

my %h;
my @seq = qw/T A C G U N/;
my $X = "\x1b[38;5;";

$h{$seq[$_]} = "\x1b[".(31+$_).";1m".$seq[$_] for 0..$#seq;

my $l;
while ($l = <>) { # print header, if present
	last if $l !~ /^@/;
	print $l;
}

while ($l) {
	my @e = split(/\t/, $l);
	$e[9] = join('', map { $h{$_} } split(//, $e[9]))."\x1b[0m";
	$e[10] = $X.join($X, (map { int(232+(ord($_)-33)/3)."m$_" } split(//, $e[10])))."\x1b[0m";
	print join("\t", @e);
} continue {
	$l = <>;
}

