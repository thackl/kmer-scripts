#!/usr/bin/env perl
use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;
use Log::Log4perl qw(:easy :no_extra_logdie_message);

use Data::Dumper;

use Kmer;
use Jellyfish;
use Fasta::Parser;
use Fastq::Parser;

our $VERSION = '0.02';


my %opt;

=head1 NAME 

mask_reads.pl

=cut

=head1 DESCRIPTION

=cut

=head1 USAGE

  $ pv <READS.f[aq]> | read_mask.pl --coords <READS.coords>  > READS.masked.f[aq]

=cut

=head1 CHANGELOG

=cut

=head2 0.02

=over

=item [PodFix] USAGE, obsolete parameter ...

=item [BugFix] Read first --read record before looping.

=item [Feature] Log::Log4perl, --quiet

=item [Refacture] Digest one-line-per-read output of read_cov2coords.pl-0.02.

=back

=head2 0.01

=over

=item [Initial]

=back

=cut

=head1 OPTIONS

=over 25

=cut

=item --coords=<file>

path to coordinates file for reads to mask. 

=cut

$opt{'coords=s'} = \(my $opt_coords);

=item --reads=<file>

the read file in fastq/fasta format, default is STDIN, set explicitly with '-'.

=cut

$opt{'reads=s'} = \(my $opt_reads = undef);

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



=head1 CODE

=cut


GetOptions(%opt) or pod2usage(1);

pod2usage(1) if($opt_help);
pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION|USAGE|OPTIONS|LIMITATIONS|AUTHORS") if($opt_man);

pod2usage(-msg => "--coords required", -verbose => 0) unless ($opt_coords);
$opt_reads = undef if (defined $opt_reads && $opt_reads eq '-');

Log::Log4perl->init(\<<'CFG');
	log4perl.logger.main				= DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{MM-dd HH:mm:ss}] [rkf] %m%n

CFG

my $L = Log::Log4perl::get_logger();
$L->level($opt_verbose ? $INFO : $WARN);

unless(defined $opt_reads){
	$opt_verbose && $L->debug("Reading from STDIN\n");
}

my $fp;
my $is_fastq;
# try fastq;
if($fp = Fastq::Parser->new(file => $opt_reads)->check_format){
	$is_fastq = 1;
}elsif($fp = Fasta::Parser->new(file => $opt_reads)->check_format){
	$is_fastq = 0;
}else{
	$L->logcroak("--reads $opt_reads neither FASTQ nor FASTA");
}

open(COORDS, '<', $opt_coords) or $L->logcroak("$!: $opt_coords");

while(
	defined(my $c = <COORDS>) && 
	(my $fs = $fp->next_seq)
){
    chomp $c;
    my ($c_id, @cs) = split(/[\t,]/, $c);
	
	if($fs->id ne $c_id){
		$L->logcroak("--reads and --coords ids differ: ".$fs->id." !=  $c_id");
	}
	
	my $seq = lc($fs->seq);
	for(my $i=0; $i<@cs; $i+=2){
		substr($seq, $cs[$i], $cs[$i+1]) = uc(substr($seq, $cs[$i], $cs[$i+1]));
	}
	$fs->seq($seq);
	print "$fs";
}

close COORDS;














