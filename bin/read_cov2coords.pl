#!/usr/bin/env perl
use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

use Fasta::Parser;
use Fastq::Parser;

our $VERSION = '0.02';

my %opt;

=head1 NAME 

reads_cov2coords.pl

=cut

=head1 DESCRIPTION

=cut

=head1 USAGE

  $ pv <READS.cov> |  read_cov2coords.pl -u 1000 -l 20 -s 50 -k 19  > READS.coords

=cut

=head1 CHANGELOG

=cut

=head2 0.02

=over

=item [PodFix] USAGE, obsolete parameter ...

=item [Refracture] Output only one tab delimited line foreach read of format
 ID OFFSET,LENGTH OFFSET,LENGTH ... Also output ID only lines for reads w/o
 annotation.

=back

=head2 0.01

=over

=item [Initial]

=back

=cut

=head1 OPTIONS

=over 25

=cut

=item -c|--cov=<file>

Path to reads coverage file, default STDIN. 

=cut

$opt{'c|cov=s'} = \(my $opt_cov);

#=item --reads=<file>
#
#the read file in fastq/fasta format, default is STDIN, set explicitly with '-'.
#
#=cut
#
#$opt{'reads=s'} = \(my $opt_reads = undef);

=item -u|--upper=<INT>

Upper cutoff for kmer frequency.

=cut

$opt{'u|upper=s'} = \(my $opt_upper = 100);

=item -l|--lower=<INT>

Lower cutoff for kmer frequency.

=cut

$opt{'l|lower=s'} = \(my $opt_lower = 10);

=item -s|--min-stretch-length=<INT>

Minimum length of a stretch to be reported.

=cut

$opt{'s|min-stretch-length=s'} = \(my $opt_min_stretch_length = 100);

=item -k|--kmer-size=<INT>

kmer size.

=cut

$opt{'k|kmer-size=s'} = \(my $opt_kmer_size = 19);

=item [--[no]verbose] 

verbose is default.

=cut

$opt{'verbose!'} = \(my $opt_verbose = 1);

=item [--help] 

show help

=cut

$opt{'help|?'} = \(my $opt_help);

=item [--man] 

show man page

=cut

$opt{'man'} = \(my $opt_man);

=back

Output: one line per read, fields tab separated, offset and length comma separated, id only lines are allowed

  read1   8,130   342,789
  read2
  read3   3,343
  ...

=head1 CODE

=cut


GetOptions(%opt) or pod2usage(1);

pod2usage(1) if($opt_help);
pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION|USAGE|OPTIONS|LIMITATIONS|AUTHORS") if($opt_man);

my $covf;
if($opt_cov){
	open ($covf, '<', $opt_cov) or die "$!: $opt_cov";
}else{
	$covf = \*STDIN;
}

while (<$covf>){
	my ($id, $med, $cov) = split(/\t/, $_);
	my @cov = split(/ /, $cov);
	my @coords;

	my $BM = 2; # 0000_0010
	if( $cov[0] > $opt_upper ){
		$BM = $BM << 1; 
	}elsif( $cov[0] < $opt_lower ){
		$BM = $BM >> 1; 
	}

	push @coords, 0 if $BM == 2;

	my $bm;
	for( my $i=1; $i<@cov; $i++ ){
		
		$bm = $cov[$i] > $opt_upper ? 4
			: ($cov[$i] < $opt_lower	? 1	: 2);
		
		# state switch
		if ($BM != $bm){
			if( $bm == 2 ){ # start of range
				push @coords, $i
			}elsif($BM == 2){ # end of range
				push @coords, $i-1
			}
			$BM=$bm;
		}
	}
	
	# add last pos if ended in range
	push @coords, $#cov if ($bm == 2 && $coords[$#coords] < @cov);
	
	# merge overlapping
	my $j=0;
	while($j < @coords-3){
		if((my $ovl=$coords[$j]+$coords[$j+1]-$coords[$j]+$opt_kmer_size-$coords[$j+2]) >= 0){
			$coords[$j+1]=$coords[$j+3];
			splice(@coords, $j+2, 2);  
		}else{
			$j+=2;
		}
	}	
	
	my $s = '';
	for (my $j=0; $j<@coords; $j+=2){
		if((my $len = $coords[$j+1]-$coords[$j]+$opt_kmer_size) >= $opt_min_stretch_length){
			$s.= "\t$coords[$j],$len";
		}
	}
	print "$id$s\n";
}

close $covf;



















