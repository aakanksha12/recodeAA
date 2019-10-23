# recodeAA

### Simple amino acid recoding program for phylogenetics
### - reduced amino acid alphabets and RY-coding
-----------------------------------------
Program developed for: Pandey, A., & Braun, E. L. (submitted). Phylogenetic analyses of 
sites in different protein structural environments result in distinct placements of the 
metazoan root. Biology

### Usage:
-----------------------------------------
"recode_protMSA.pl" takes a relaxed phylip format (non-interleaved) protein sequence
alignment as input and generates recoded files. The recoding can involve simple changes
to the alphabet or back-translation, either with or without RY-coding.

Simply invoke the program as follows:

./recode_protMSA.pl prefix alphabet authority

prefix = first part of the name for the input and output file
  - prefix.phy = input file, assumed to be relaxed phylip format amino acid alignment
  - prefix.alphabet.phy = output file, relaxed phylip format recoded dataset
  - authority = optional authority file to organize output file

Alternative alphabets that are currently implemented are stored as arrays in the section
of the program called "Initialize arrays for codons and alternative alphabets". They are:

  - RY universal code (F=111, L=11?, I=01?, etc) (note: purine = 0; pyrimidine = 1)
  - Universal code with ambiguities (F=TTY, L=YT?, I=AT?, etc)
  - Dayhoff codes (six-state, 0-5)
  - Susko-Roger "saturation codes" (six-state, 0-5)
  - Kosiol et al. PAM codes (six-state, 0-5)
  - Kosiol et al. WAG codes (six-state, 0-5)
  - Hanada codes (four states, coded as nucleotides)

It is straightforward to add new reduced alphabets or genetic codes as arrays. Simply
inspect the source code in this section and modify the if...elsif... commands at the end
of this section. The default is to recode to the universal code with ambiguities.

The authority file is simply the number of taxa followed by a list of taxon names. The
recoded data will be output in the order they are listed in the authority file. If the
authority file name is omitted the input file (i.e., "prefix.phy") will be treated as the
authority file (i.e., the output file will have the same taxa as the input file and they
will be output in the same order).


-----------------------------------------
References:

Dayhoff, M. O., Schwartz, R. M., & Orcutt, B. C. (1978). A model of evolutionary change in 
proteins. In Atlas of protein sequence and structure (Vol. 5, pp. 345-352). Dayhoff, M. O.
(ed.) National Biomedical Research Foundation, Silver Spring MD.

Hanada, K., Shiu, S. H., & Li, W. H. (2007). The nonsynonymous/synonymous substitution rate 
ratio versus the radical/conservative replacement rate ratio in the evolution of mammalian 
genes. Mol. Biol. Evol. , 24(10), 2235-2241.

Kosiol, C., Goldman, N., & Buttimore, N.H. (2004). A new criterion and method for amino 
acid classification. J. Theor. Biol. 228, 97–106.

Susko, E., & Roger, A.J. (2007). On reduced amino acid alphabets for phylogenetic inference. 
Mol. Biol. Evol. 24, 2139–2150.
