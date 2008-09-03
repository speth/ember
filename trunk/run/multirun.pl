#!/usr/bin/perl

use strict;
use Getopt::Long;

my $nProcs = 4;
my $commandfile;
GetOptions( "commands=s" => \$commandfile,
			"n=i" => \$nProcs);

open(INFILE,"<$commandfile");
my @commands = <INFILE>;
close(INFILE);

my @pids;
my $runningProcs = 0;

print "Total tasks to run: " . ($#commands+1) . "\n";
print "Simultaneous tasks: $nProcs\n";

while ($#commands >= 0 && $runningProcs < $nProcs) {
    my $cmd = shift @commands;
	chomp $cmd;
    $runningProcs++;
    sleep(1);
    my $pid = fork();
    push @pids, $pid;
    if (not defined $pid) {
		print "Failure!\n";
    } elsif ($pid == 0) {
		print "executing command \"" . $cmd . "\"\n";
	exec($cmd);
		print "completed command \"" . $cmd . "\"\n";
    }
}


foreach my $oldpid (@pids) {
	wait();
	sleep(1);
    if ($#commands >= 0) {
		my $cmd = shift @commands;
		chomp $cmd;
		my $pid = fork();
		push @pids, $pid;
		if (not defined $pid) {
			print "Failure!\n";
		} elsif ($pid == 0) {
			print "executing command \"" . $cmd . "\"\n";
			exec($cmd);
			print "completed command \"" . $cmd . "\"\n";
		}    
	}
}
