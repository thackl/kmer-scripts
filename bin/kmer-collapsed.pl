#!/usr/bin/env perl
use warnings;
no warnings 'qw';
use strict;

use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

use Log::Log4perl qw(:no_extra_logdie_message);
use Log::Log4perl::Level;

use FindBin qw($RealBin);
use lib "$RealBin/../lib/";

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

# additional modules

use Jellyfish 0.04; # get_kmer_size
use Kmer;
use Verbose::ProgressBar;

#------------------------------------------------------------------------------#
# Globals

our $VERSION = '0.03';

# get a logger
my $L = Log::Log4perl::get_logger();
Log::Log4perl->init(\<<'CFG');
	log4perl.rootLogger                 = DEBUG, Screen
	log4perl.appender.Screen			= Log::Log4perl::Appender::Screen
	log4perl.appender.Screen.stderr		= 0
	log4perl.appender.Screen.layout		= PatternLayout
	log4perl.appender.Screen.layout.ConversionPattern = [%d{MM-dd HH:mm:ss}] [kcd] %m%n
CFG
$L->level($INFO);


my %opt = (
    c1 => undef,
    c2 => undef,
    x1 => 1,
    x2 => 1,
    "window-size" => 1,
);


#------------------------------------------------------------------------------#
# GetOptions
Getopt::Long::Configure("no_ignore_case");
GetOptions(\%opt, qw(
	c1|1=s
	c2|2=s
	x1|x=i
	x2|y=i
        window-size=i
	quiet
	debug
	help|h
	version
)) or $L->logcroak($!);

pod2usage(1) if $opt{help};

# version
if($opt{version}){
	print "$VERSION\n"; 
	exit 0;
}

$opt{quiet} && $L->level($WARN);
$opt{debug} && $L->level($DEBUG);

$L->debug("GetOptions:\n", Dumper(\%opt));



##------------------------------------------------------------------------##	
# required	
for(qw(c1 c2)){
    if(ref $opt{$_} eq 'ARRAY'){
	pod2usage("required: --$_") unless @{$opt{$_}}
    }else{
	pod2usage("required: --$_") unless defined ($opt{$_})
    }
};


# check fq files -e -s
my $c1 = $opt{c1};
my $c2 = $opt{c2};
my $x1 = $opt{x1};
my $x2 = $opt{x2};
my $ws = $opt{"window-size"};

foreach my $file($c1, $c2){
	$L->logcroak("Cannot find file: $file ") unless -e $file && -f $file;
}

open(my $cf1, '<', $c1) or $L->logcroak("$c1: $!");
open(my $cf2, '<', $c2) or $L->logcroak("$c2: $!");

##------------------------------------------------------------------------##	
# read both files simultaneously
my ($s1, $s2);
while( 
    defined ($s1 = <$cf1>) && 
    defined ($s2 = <$cf2>)
){

    chomp($s1, $s2);

    my($id1, $med1, @cov1) = split(/\s/,$s1);
    my($id2, $med2, @cov2) = split(/\s/,$s2);

    # lets play it safe
    $L->logcroak("Inconsistent coverage files (line $.)") if ($id1 ne $id2) or (@cov1 != @cov2);
    
    my @covr;
    if($ws < 2){
	for(my $i=0; $i<@cov1; $i++){
	    $covr[$i] = sprintf("%.02f", ($cov1[$i]/$x1) / ($cov2[$i]/$x2));
	}
    }else{
	next unless @cov1>=$ws; # ignore short seqs
	my $wn=$ws-1;
	my $wh=int($ws/2);
	my ($wcov1, $wcov2);
	for(my $i=0; $i<@cov1-$ws; $i+=$ws){
	    $wcov1 = (sort{$a<=>$b}(@cov1[$i..$i+$wn]))[$wh];
	    $wcov2 = (sort{$a<=>$b}(@cov2[$i..$i+$wn]))[$wh];
    	    push @covr, sprintf("%.02f", ($wcov1/$x1) / ($wcov2/$x2) );
	}
    }
    print $id1, "\t", "MED", "\t", join(" ", @covr), "\n";
}


