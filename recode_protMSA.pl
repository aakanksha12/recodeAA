#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Initialize variables
############################################################################

my($progname) = $0;
my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($zter);
my($tvar1);
my($tvar2);

my @tmparr1;
my @tmparr2;

if ( @ARGV < 2 || @ARGV > 3 ) {
	print "Usage:\n  \$ $progname <prefix> <alphabet> <authority>\n";
	print "  prefix    = prefix for input and output files\n";
	print "      use '.phy' as the extension for the input file\n";
	print "  alphabet  = codon dictionary or reduced alphabet\n";
	print "      univRY  -- universal code, binary nucleotides\n";
	print "      univNT  -- universal code, nucleotides (A, C, G, T, and ?)\n";
	print "      Dayhoff -- six-state Dayhoff recoding\n";
	print "      SR6     -- six-state S&R recoding\n";
	print "      KGB6pam -- six-state KGB recoding from PAM matrix\n";
	print "      KGB6wag -- six-state KGB recoding from WAG matrix\n";
	print "      Hanada  -- four-state Hanada recoding\n";
	print "  authority = authority file name (optional, default to input file)\n";
	print "exiting...\n";
	exit;
}

my($prefix) = $ARGV[0];
my($suffix) = $ARGV[1];

my($seqfname)  = "$prefix" . ".phy";
my($authfname) = "$prefix" . ".phy";
if ( @ARGV == 3 ) {
	$authfname = $ARGV[2];
}
# print "$authfname\n";

############################################################################
# Initialize arrays for codons and alternative alphabets
############################################################################

# Flag for whether the output will be codons or recoded amino acids
my($tripletflag) = 3; 
# amino acid alphabet = 1
# codons              = 3

# Amino acid numbers and universal code
my @aminoacid = ("F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G");
# 0.  F -- YYY (UUU, UUC)
# 1.  L -- YYN (UUA, UUG, CUN)
# 2.  I -- RYN (AUU, AUC, AUA)
# 3.  M -- RYR (ATG)
# 4.  V -- RYN (GUN)
		
# 5.  S -- NNN (UCN, AGU, AGC)
# 6.  P -- YYN (CCN)
# 7.  T -- RYN (ACN)
# 8.  A -- RYN (GCN)
		
# 9.  Y -- YRY (UAU, UAC)
# 10. H -- YRY (CAU, CAC)
# 11. Q -- YRR (CAA, CAG)
# 12. N -- RRY (AAU, AAC)
# 13. K -- RRR (AAA, AAG)
# 14. D -- RRY (GAU, GAC)
# 15. E -- RRR (GAA, GAG)
		
# 16. C -- YRY (UGU, UGC)
# 17. W -- YRR (UGG)
# 18. R -- NRN (CGN, AGA, AGG)
# 19. G -- RRN (GGN)

# Universal code, with  R = 0; Y = 1
#  AA order:   F      L      I      M      V      S      P      T      A      Y      H      Q      N      K      D      E      C      W      R      G
my @univRY = ("111", "11?", "01?", "010", "01?", "???", "11?", "01?", "01?", "101", "101", "100", "001", "000", "001", "000", "101", "100", "?0?", "00?");
#  RY codon:   YYY    YYN    RYN    RYR    RYN    NNN    YYN    RYN    RYN    YRY    YRY    YRR    RRY    RRR    RRY    RRR    YRY    YRR    MRN    RRN

# Universal code, nucleotides with ambiguity
#  AA order:   F      L      I      M      V      S      P      T      A      Y      H      Q      N      K      D      E      C      W      R      G
my @univNT = ("TTY", "YT?", "AT?", "ATG", "GT?", "???", "CC?", "AC?", "GC?", "TAY", "CAY", "CAR", "AAY", "AAR", "GAY", "GAR", "TGY", "TGG", "?G?", "GG?");
#  RY codon:   YYY    YYN    RYN    RYR    RYN    NNN    YYN    RYN    RYN    YRY    YRY    YRR    RRY    RRR    RRY    RRR    YRY    YRR    MRN    RRN

# Dayhoff codes (six state)
#  AA order:    F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @dayhoff = ("5", "4", "4", "4", "4", "1", "1", "1", "1", "5", "3", "2", "2", "3", "2", "2", "0", "5", "3", "1");
#   0 = C
#   1 = AGPST
#   2 = NDEQ
#   3 = RHK
#   4 = ILMV
#   5 = FWY

# S&R6 saturation codes (six state) 
# from Susko, E., and Roger, A.J. (2007). On reduced amino acid alphabets for phylogenetic 
# inference. Mol. Biol. Evol. 24, 2139–2150.
#  AA order:    F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @sr6sat  = ("5", "3", "3", "3", "3", "0", "0", "0", "0", "5", "0", "2", "1", "2", "1", "1", "4", "4", "2", "1");
#   0 = APST
#   1 = DENG
#   2 = QKR
#   3 = ILMV
#   4 = WC
#   5 = FYH

# KGB PAM codes (six state) 
# from Kosiol, C., Goldman, N., and Buttimore, N.H. (2004). A new criterion and method for 
# amino acid classification. J. Theor. Biol. 228, 97–106.
#  AA order:     F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @kgb6pam  = ("4", "2", "2", "2", "5", "0", "0", "1", "0", "4", "1", "1", "1", "1", "1", "1", "5", "3", "1", "0");
#   0 = AGPS
#   1 = DENQHKRT
#   2 = MIL
#   3 = W
#   4 = FY
#   5 = CV

# KGB WAG codes (six state) 
# from Kosiol, C., Goldman, N., and Buttimore, N.H. (2004). A new criterion and method for 
# amino acid classification. J. Theor. Biol. 228, 97–106.
#  AA order:     F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @kgb6wag  = ("4", "3", "3", "3", "2", "0", "0", "0", "0", "4", "0", "0", "0", "0", "0", "0", "2", "5", "0", "1");
#   0 = HRKQNEDSTPA
#   1 = G
#   2 = CV
#   3 = IML
#   4 = FY
#   5 = W

# Hanada codes (four state, recoded as nucleotides)
#  AA order:   F    L    I    M    V    S    P    T    A    Y    H    Q    N    K    D    E    C    W    R    G
my @hanada = ("G", "C", "C", "C", "C", "A", "A", "A", "A", "G", "G", "G", "A", "G", "T", "T", "A", "G", "G", "A");
#   A = ANCGPST
#   C = ILMV
#   G = RQHKFWY
#   T = DE

my @alphabet;
for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $univNT[$iter]; } # default is univNT
if ( lc($suffix) eq "univry" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $univRY[$iter]; } }
elsif ( lc($suffix) eq "dayhoff" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $dayhoff[$iter]; } $tripletflag = 1; }
elsif ( lc($suffix) eq "sr6" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $sr6sat[$iter]; } $tripletflag = 1; }
elsif ( lc($suffix) eq "kgb6pam" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $kgb6pam[$iter]; } $tripletflag = 1; }
elsif ( lc($suffix) eq "kgb6wag" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $kgb6wag[$iter]; } $tripletflag = 1; }
elsif ( lc($suffix) eq "hanada" ) { for ($iter=0; $iter<20; $iter++) { $alphabet[$iter] = $hanada[$iter]; } $tripletflag = 1; }

# for ( $iter=0; $iter<20; $iter++) {
#	print "  $iter";
#	print ". $aminoacid[$iter] -- $alphabet[$iter]\n";
# }

############################################################################
# Read the authority file
############################################################################

my @authlist;
open (my $AUTHF, "$authfname") or die "Could not open file < $authfname > for input.\n";
@authlist = <$AUTHF>; # Read the authority file
close($AUTHF) or die "Could not close file < $authfname >\n";
my($authlistnum) = $#authlist + 1;

my($ntax);
my @taxname;
my($maxnamelen) = 10;

chomp($authlist[0]);
($ntax) = split(/\s+/, $authlist[0]);
# print "$ntax\n";
for ($iter=0; $iter<$ntax; $iter++) {
	$jter = $iter + 1;
	chomp($authlist[$jter]);
	($taxname[$iter]) = split(/\s+/, $authlist[$jter]);
	if ( length($taxname[$iter]) > $maxnamelen ) { $maxnamelen = length($taxname[$iter]); }
#	print "  $taxname[$iter] -- $maxnamelen\n";
}

############################################################################
# Read the sequence file
############################################################################
my @seqlist;
open (my $SEQF, "$seqfname") or die "Could not open file < $seqfname > for input.\n";
@seqlist = <$SEQF>; # Read the relaxed phylip format input file
close($SEQF) or die "Could not close file < $seqfname >\n";
my($seqlistnum) = $#seqlist + 1;
for ($iter=0; $iter<$seqlistnum; $iter++) { chomp($seqlist[$iter]); }

my($nchar);

($tvar1,$nchar) = split(/\s+/, $seqlist[0]);
my($totalchar) = $nchar * $tripletflag;

############################################################################
# Output the recoded data
############################################################################
my($taxonfound);	my($padding);
my($seqname);		my($sequence);
my($recodedseq);

my($outfname) = "$prefix" . "." . "$suffix" . ".phy";
open (my $OUTF, ">$outfname") or die "Could not open file < $outfname > for output.\n";
print $OUTF "$ntax $totalchar\n";

# iterate through taxon names from authority file
for ($iter=0; $iter<$ntax; $iter++) {
	
	# output the name and pad with spaces
	print $OUTF "$taxname[$iter]";
	$padding = $maxnamelen - length($taxname[$iter]) + 2;
	for ($jter=0; $jter<$padding; $jter++) { print $OUTF " "; }
	
	# search for the taxon name
	$taxonfound = 0;
	for ($kter=1; $kter<$seqlistnum; $kter++) {
		($seqname,$sequence) = split(/\s+/, $seqlist[$kter]);
		if ( $seqname eq $taxname[$iter] ) {
			$recodedseq = recodeprotein($sequence);
			print $OUTF "$recodedseq\n";
			$taxonfound = 1;
			$kter = $seqlistnum;
		}
	}
	
	# write missing data for taxa that are not found
	if ( $taxonfound == 0 ) {
		for ($lter=0; $lter<$totalchar; $lter++) { print $OUTF "?"; }
		print $OUTF "\n";
	}
}

close($OUTF) or die "Could not close file < $outfname >\n";
print "Recoded file exported as $outfname\n";

exit;

############################################################################
# Subroutine to recode the protein sequence, either by reverse translating
# or using a reduced alphabet
############################################################################
sub recodeprotein {
	my($localinseq) = @_;
	my($aalen) = length($localinseq);

	my $localoutseq;
	my($aa);	my($codon);
	my($locfirstaa) = 1;
	
	my($allmissing) = "???";
	my($allgap) = "---";
	if ( $tripletflag == 1 ) { $allmissing = "?"; $allgap = "-"; }

	for ($zter=0; $zter<$aalen; $zter++) {

		$aa = substr($localinseq,$zter,1);
		
		# then back translate the amino acid
		$codon = "$allmissing";
		if ( $aa eq "-" ) { $codon = "$allgap"; }
		elsif ( uc($aa) eq "F" ) { $codon = "$alphabet[0]";  } # amino acid 0
		elsif ( uc($aa) eq "L" ) { $codon = "$alphabet[1]";  } # amino acid 1
		elsif ( uc($aa) eq "I" ) { $codon = "$alphabet[2]";  } # amino acid 2
		elsif ( uc($aa) eq "M" ) { $codon = "$alphabet[3]";  } # amino acid 3
		elsif ( uc($aa) eq "V" ) { $codon = "$alphabet[4]";  } # amino acid 4
		
		elsif ( uc($aa) eq "S" ) { $codon = "$alphabet[5]";  } # amino acid 5
		elsif ( uc($aa) eq "P" ) { $codon = "$alphabet[6]";  } # amino acid 6
		elsif ( uc($aa) eq "T" ) { $codon = "$alphabet[7]";  } # amino acid 7
		elsif ( uc($aa) eq "A" ) { $codon = "$alphabet[8]";  } # amino acid 8
		
		elsif ( uc($aa) eq "Y" ) { $codon = "$alphabet[9]";  } # amino acid 9
		elsif ( uc($aa) eq "H" ) { $codon = "$alphabet[10]"; } # amino acid 10
		elsif ( uc($aa) eq "Q" ) { $codon = "$alphabet[11]"; } # amino acid 11
		elsif ( uc($aa) eq "N" ) { $codon = "$alphabet[12]"; } # amino acid 12
		elsif ( uc($aa) eq "K" ) { $codon = "$alphabet[13]"; } # amino acid 13
		elsif ( uc($aa) eq "D" ) { $codon = "$alphabet[14]"; } # amino acid 14
		elsif ( uc($aa) eq "E" ) { $codon = "$alphabet[15]"; } # amino acid 15
		
		elsif ( uc($aa) eq "C" ) { $codon = "$alphabet[16]"; } # amino acid 16
		elsif ( uc($aa) eq "W" ) { $codon = "$alphabet[17]"; } # amino acid 17
		elsif ( uc($aa) eq "R" ) { $codon = "$alphabet[18]"; } # amino acid 18
		elsif ( uc($aa) eq "G" ) { $codon = "$alphabet[19]"; } # amino acid 19
		
		if ( $locfirstaa == 1 ) { 
			$localoutseq = "$codon";
			$locfirstaa = 0;
		}
		else { 
			$localoutseq = "$localoutseq" . "$codon";
		}
	}
	
	return $localoutseq;
}




