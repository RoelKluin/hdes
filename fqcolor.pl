#!/usr/bin/perl
# Usage `samtools view $bam | csam | less -R'
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

my (%h, %g);
my @seq = qw/T A C G U N/;
my $X = "\e[38;5;";

my $phred = 33;
my ($help, $man, $regex);

Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-)");
GetOptions('help|?|h' => \$help, man => \$man, 'phred=i' => \$phred, 're=s' => \$regex);

pod2usage(-verbose => 2) if defined $man;
pod2usage(0) if defined $help;


$h{$seq[$_]} = "\e[".(31+$_).";1m".$seq[$_] for 0..$#seq;
$g{$seq[$_]} = "\e[".(91+$_).";1m".$seq[$_] for 0..$#seq;

while (<>) {
        print and next if /^(==|$)/;
        #if (/^([^ ]+) (.*[ACTGN]{4,})/) {
        #	print $1.join('', map { $_ ? $h{$_} // $_ : () } split(//, substr ($2, 0)))."\x1b[0m\n";
        #} else {
                my $re = $regex ? $regex : (/\s([ACTG]+(\|[ACTG]+)*)\s/ ? qr/($1)?/ : qr//);
		print;
        #}
        print and next while (($_ = <>) =~ /^(==|$)/);
	chomp;
	print join('', map { $_ ? (exists $h{$_} ? $h{$_} : "\e[7m".join('', map { exists $g{$_}? $g{$_} : warn "\n$_\n" } split(//, $_))."\e[0m") : () } split($re))."\e[0m\n";
        print and next while (($_ = <>) =~ /^(==|$)/);
	print;
        print and next while (($_ = <>) =~ /^(==|$)/);
	chomp;
	print $X.join($X, (map { int(232+(ord($_)-$phred)/3)."m$_" } split(//)))."\e[0m\n";
}
__END__

=head1 NAME

fqcolor.pl - print fastq with coloring

=head1 SYNOPSIS

fqcolor.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-h|-?|--help>

Print options message

=item B<--man>

Print manual page

=item B<--phred> = i

Use this phred offset

=item B<--re> = s

highlight this regular expression

=back

=head1 DESCRIPTION

fqcolor.pl - print fastq with coloring

=head1 EXAMPLE


=cut


