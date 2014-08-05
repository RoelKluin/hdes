#!/usr/bin/perl
# Usage `samtools view $bam | csam | less -R'

use warnings;
use strict;

my (%h, %g);
my @seq = qw/T A C G U N/;
my $X = "\e[38;5;";

$h{$seq[$_]} = "\e[".(31+$_).";1m".$seq[$_] for 0..$#seq;
$g{$seq[$_]} = "\e[".(91+$_).";1m".$seq[$_] for 0..$#seq;

while (<>) {
        print and next if /^(==|$)/;
        #if (/^([^ ]+) (.*[ACTGN]{4,})/) {
        #	print $1.join('', map { $_ ? $h{$_} // $_ : () } split(//, substr ($2, 0)))."\x1b[0m\n";
        #} else {
                my $re = /\s([ACTG]+(\|[ACTG]+)*)\s/ ? qr/($1)?/ : qr//;
		print;
        #}
        print and next while (($_ = <>) =~ /^(==|$)/);
	chomp;
	print join('', map { $_ ? (exists $h{$_} ? $h{$_} : "\e[7m".join('', map { exists $g{$_}? $g{$_} : warn "\n$_\n" } split(//, $_))."\e[0m") : () } split($re))."\e[0m\n";
        print and next while (($_ = <>) =~ /^(==|$)/);
	print;
        print and next while (($_ = <>) =~ /^(==|$)/);
	chomp;
	print $X.join($X, (map { int(232+(ord($_)-33)/3)."m$_" } split(//)))."\e[0m\n";
}

