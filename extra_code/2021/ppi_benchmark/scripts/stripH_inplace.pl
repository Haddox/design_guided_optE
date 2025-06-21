#!/usr/bin/perl
##
##
###############################################################################

use strict;

if ($#ARGV < 0) {
	print STDERR "usage: $0 <pdbfile>\n";
	exit -1;
}
my $pdbfile = shift @ARGV;

my @file;
open (PDB, $pdbfile);
while (my $line = <PDB>) {
	push @file, $line;
}
close (PDB);

open (PDBR, ">$pdbfile");
foreach my $line (@file) {
	if ($line !~ /^ATOM/ && $line !~ /^HETATOM/) {
		print PDBR $line;
		#print $line;
		next;
	}
	if ( substr ($line, 13, 1) ne "H" &&  substr ($line, 12, 1) ne "H" ) {
		print PDBR $line;
		#print $line;
	}
}
#close (PDBR);
