#!/usr/bin/perl

use strict;
use File::Copy;
use Getopt::Long;

my ($inputfile, $prefix, $suffix, $subdir, $phimin, $phimax, $dphi, $fuel);
#my $inputfile = "base-input.txt";
#my $prefix = "syngas-x0.50-phi";
#my $suffix = ".txt";
#my $subdir = "run";

GetOptions( "phimin=s" => \$phimin,
			"phimax=s" => \$phimax,
			"dphi=s" => \$dphi,
			"if=s" => \$inputfile,
			"prefix=s" => \$prefix,
			"suffix=s" => \$suffix,
			"subdir=s" => \$subdir,
			"fuel=s" => \$fuel);

if (not defined $inputfile || 
	not defined $phimin || 
	not defined $phimax || 
	not defined $dphi || 
	not defined $subdir ||
	not defined $suffix ||
	not defined $prefix) 
{
	print "Usage: All of the following arguments must be supplied:\n" .
		"\t--if: A valid 1dflameV2 input file\n" .
		"\t--prefix: A string to use as the prefix for the generated input files\n" .
		"\t--suffix: A string to use as the suffix for the generated input files\n" .
		"\t--phimin: The minimum equivalence ratio\n" .
		"\t--phimax: the maximum equivalence ratio\n" .
		"\t--dphi: the interval between equivalence ratios\n" .
		"\t--fuel: specification of the fuel composition\n\n";
	
	exit;
}

# read the original input file;

open(INFILE,"<$inputfile");
my @lines = <INFILE>;

my $phi = $phimin;
while ($phi <= $phimax) {
	my $phistr = sprintf  "%4.2f", $phi;
    my $newname = $prefix . $phistr . $suffix;
	my $dirname = $subdir . "/" . $prefix . $phistr;

	print "Creating file:\"$newname\"\n";
	open(OUTFILE,">$newname");
	foreach my $line (@lines) {
		if ($line =~ /outputDir/) {
			$line =~ s/\".*\"/\"$dirname\"/;	
		} elsif ($line =~ /equivalenceRatio/) {
			$line =~ s/=.*;/= $phistr;/;
		} elsif ($line =~ /fuel/) {
			$line =~ s/=.*;/= \"$fuel\";/;
		}
		
		print OUTFILE $line;
	}
	close(OUTFILE);
	$phi += $dphi;
}
