#!/usr/bin/perl
# Usage `samtools view $bam | csam | less -R'

use warnings;
use strict;

my %h;
my @seq = qw/T A C G U N/;
my $X = "\x1b[38;5;";

$h{$seq[$_]} = "\x1b[".(31+$_).";1m".$seq[$_] for 0..$#seq;

while (<>) {
        print and next if /^==/;
        #if (/^([^ ]+) (.*[ACTGN]{4,})/) {
        #	print $1.join('', map { $_ ? $h{$_} // $_ : () } split(//, substr ($2, 0)))."\x1b[0m\n";
        #} else {
		print;
        #}
        print and next while (($_ = <>) =~ /^==/);
	chomp;
	print join('', map { $h{$_} } split(//))."\x1b[0m\n";
        print and next while (($_ = <>) =~ /^==/);
	print;
        print and next while (($_ = <>) =~ /^==/);
	chomp;
	print $X.join($X, (map { int(232+(ord($_)-33)/3)."m$_" } split(//)))."\x1b[0m\n";
}

