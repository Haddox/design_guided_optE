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

my $KEEPWAT=0;

## amino acid map
my %one_to_three = (
	'G' => 'GLY', 'A' => 'ALA', 'V' => 'VAL', 'L' => 'LEU',
	'I' => 'ILE', 'P' => 'PRO', 'C' => 'CYS', 'M' => 'MET',
	'H' => 'HIS', 'F' => 'PHE', 'Y' => 'TYR', 'W' => 'TRP',
	'N' => 'ASN', 'Q' => 'GLN', 'S' => 'SER', 'T' => 'THR',
	'K' => 'LYS', 'R' => 'ARG', 'D' => 'ASP', 'E' => 'GLU' );

my %three_to_one = (
	'GLY' => 'G', 'ALA' => 'A', 'VAL' => 'V', 'LEU' => 'L', 'ILE' => 'I',
	'PRO' => 'P', 'CYS' => 'C', 'MET' => 'M', 'HIS' => 'H', 'PHE' => 'F',
	'TYR' => 'Y', 'TRP' => 'W', 'ASN' => 'N', 'GLN' => 'Q', 'SER' => 'S',
	'THR' => 'T', 'LYS' => 'K', 'ARG' => 'R', 'ASP' => 'D', 'GLU' => 'E',

	# nonstd
	'5HP' => 'Q', 'ABA' => 'C', 'AGM' => 'R', 'CEA' => 'C', 'CGU' => 'E',
	'CME' => 'C', 'CSB' => 'C', 'CSE' => 'C', 'CSD' => 'C', 'CSO' => 'C',
	'CSP' => 'C', 'CSS' => 'C', 'CSW' => 'C', 'CSX' => 'C', 'CXM' => 'M',
	'CYM' => 'C', 'CYG' => 'C', 'DOH' => 'D', 'FME' => 'M', 'GL3' => 'G',
	'HYP' => 'P', 'KCX' => 'K', 'LLP' => 'K', 'LYZ' => 'K', 'MEN' => 'N',
	'MGN' => 'Q', 'MHS' => 'H', 'MIS' => 'S', 'MLY' => 'K', 'MSE' => 'M',
	'NEP' => 'H', 'OCS' => 'C', 'PCA' => 'Q', 'PTR' => 'Y', 'SAC' => 'S',
	'SEP' => 'S', 'SMC' => 'C', 'STY' => 'Y', 'SVA' => 'S', 'TPO' => 'T',
	'TPQ' => 'Y', 'TRN' => 'W', 'TRO' => 'W', 'YOF' => 'Y',

	# fpd
	'MIS' => 'S'
);

my @file;
open (PDB, $pdbfile);
while (my $line = <PDB>) {
	push @file, $line;
}
close (PDB);

my $pdbout = $pdbfile;
if ($KEEPWAT==0) {
	$pdbout =~ s/\.pdb/_clean.pdb/;
} else {
	$pdbout =~ s/\.pdb/_cleanwat.pdb/;
}

open (PDBOUT, ">$pdbout");
my $linecount = 0;
foreach my $line (@file) {
	last if ($line =~ /^ENDMDL/ && $linecount>0);
	if ($line =~ /^CRYST1/) {
		print PDBOUT $line;
	}

	if ($line =~/^ATOM  / || $line =~/^HETATM/) {
		my $resname = substr($line,17,3);

		# remove multiple confs
		my $conf = substr($line,16,1);
		next if ($conf ne " " && $conf ne "A");

		if ($KEEPWAT == 1 && ($resname eq "WAT" || $resname eq "HOH")) {
			substr($line,17,3) = "WAT";
		} else {
			next if (!defined $three_to_one{ $resname });
			# sanitization
			substr($line,0,6) = "ATOM  ";
			substr($line,17,3) = $one_to_three{ $three_to_one{ $resname } };
			substr($line,56,4) = "1.00";

			# MSE
			if ($resname eq "MSE") {
				my $atomname = substr ($line,12,4);
				if ($atomname eq "SE  ") {
					substr ($line,12,4) = " SD ";
				}
			}
			# all other non-standard ... throw away sidechain
			elsif (substr($line,17,3) ne $resname) {
				my $atomname = substr ($line,12,4);
				next if ($atomname ne " C  " && $atomname ne " CA "  && $atomname ne " O  "  && $atomname ne " N  "  && $atomname ne " CB ");
			}
		}
		print PDBOUT $line;
		$linecount++;
	}
}
close (PDBR);
