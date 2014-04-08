#!/usr/bin/env perl

# this script takes in the file produced by extract_deletions.sh
# usage: cat deletions.input | perl generate_deletion_profile.pl > deletions.positions

use strict;
use warnings;

while (<>) {
    chomp;
    my @line = split(/\s+/, $_);
    for (my $i=1; $i<@line; $i++) {
	if ($line[$i] eq "D") {
	    my $k=0;
	    for (my $j=1; $j<$i; $j++) {
		if ($line[$j] ne "D") {   
		    $k=$k+$line[$j]
		}
	    }
	    for (my $l=1; $l<=$line[$i-1]; $l++) {
		print $line[0], "\t", $k+$l, "\n";
	    }
	}
    }
}
