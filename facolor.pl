#!/usr/bin/perl
# Usage `samtools view $bam | csam | less -R'
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

my (%h, %g);
my @seq = qw/T A C G U N/;

my ($help, $man, $re);

Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-)");
GetOptions('help|?|h' => \$help, man => \$man, 're=s' => \$re);

pod2usage(-verbose => 2) if defined $man;
pod2usage(0) if defined $help;
$re = defined($re) ? ($re =~ /^\(/ ? qr/$re/ : qr/($re)/) : qr//;

$h{$seq[$_]} = "\e[".(31+$_).";1m".$seq[$_] for 0..$#seq;
$g{$seq[$_]} = "\e[".(91+$_).";1m".$seq[$_] for 0..$#seq;

while (<>) {
	chomp;
        print "\e[1m$_\e[0m\n" and next if /^(?:==|\s*$|>)/;
	print join("", map { $_ =~ $re ? "\e[7m".join('',
                    map { $g{$_} // $_ } split //)."\e[0m" :
                "\e[1m".join('', map { $h{$_} // $_ } split //)."\e[0m"
                } split $re)."\e[0m\n";
}
__END__

=head1 NAME

facolor.pl - print fasta with colors

=head1 SYNOPSIS

facolor.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-h|-?|--help>

Print options message

=item B<--man>

Print manual page

=item B<--re> = s

highlight this regular expression

=back

=head1 DESCRIPTION

facolor.pl - print fastq with colors

=head1 EXAMPLE


=cut


