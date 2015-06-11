#!/usr/bin env perl
use warnings;
use strict;

my (%h, %sz);
while (<>) {
    chomp;
    my ($ct, $f) = split /\t/, $_;
    next if -d $f or $f =~ /^Makefile|main\.cpp$/;
    my ($bn, @etc) = split /\./, $f;
    $sz{$bn} = 0 if not exists $h{$bn};
    $sz{$bn} += $ct;
    push @{$h{$bn}}, $f;
}

my @keys = ("__init__", sort { scalar(@{$h{$b}}) <=> scalar(@{$h{$a}}) ||
    (join("\t", @{$h{$b}}) =~ /\.cpp/) <=> (join("\t", @{$h{$a}}) =~ /\.cpp/) ||
    $sz{$b} <=> $sz{$a} || $a cmp $b } keys %h);

$h{__init__} = ['Makefile', 'main.cpp'];

print "open ".join("\ntabnew ", map { join("\nvs ", reverse sort @{$h{$_}})} @keys)."\ntabr\n";
