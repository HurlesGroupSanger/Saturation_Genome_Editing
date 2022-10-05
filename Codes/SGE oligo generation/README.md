## This code has been replaced by [VaLiAnT](https://github.com/cancerit/VaLiAnT)

### This page is for record purpose only. Please use [VaLiAnT](https://github.com/cancerit/VaLiAnT) for the oligo design.

The scripts needed are all in this folder except for the transeq from EMBOSS. You have to install transeq (EMBOSS:6.5.7.0) by yourself. You have to ensure all of these scripts are in your path. Also, you need the parallel from gnu. Your machine may not have the gnu version but you can conda it.

1) saturated_mutagenesis_real.sh
2) saturated_mutagenesis_subreal.sh
3) rm_same_peptideseq.sh
4) translate.sh
5) fasta2tab.pl
6) RetainInList, RemoveInList, InsertInList (C++ compiled code. It should work on all Linux 64 platform. If not you can get the source code from me. These three codes are the subcode for a package that we wrote for "Widespread interaction between ADAR1 and transcriptional byproducts, Wu et al, 2019 bioRxiv" few years back in my PhD time.)


I have included a input table for the DDX3X all exon as an example, and one of its output file for ready-to-order oligo. 


**Main script: saturated_mutagenesis_real.sh**

It does: 
1) Generating working file (testchip)
2) Generating folder based on the name in column 1
3) Run the saturation mutagenesis script, saturated_mutagenesis_subreal.sh.
4) Add the constant sequence and adaptor sequence on the output of saturated_mutagenesis_subreal.sh.

**Running:**

saturated_mutagenesis_real.sh [full path for input tab file] [core]

**Dependencies:**

1) saturated_mutagenesis_subreal.sh
	Dependencies:
	1) rm_same_peptideseq.sh
		Dependencies: translate.sh -> transeq (EMBOSS), fasta2tab.pl	
	2) RetainInList, RemoveInList, InsertInList
2) gnu-parallel
3) gnu-(awk, sed, uniq, sort)

**Input tab file format:**

Column 1: Name (do not contain special character and white space and comma "," or underscore "_".)</br>
Column 2: 5intron_seq (5' intronic sequence for saturation)</br>
Column 3: exon_seq (exonic sequence for saturation, should have mutated PAM or mutated protospacer)</br>
Column 4: 3intron_seq (3' intronic sequence for saturation)</br>
Column 5: 5intron_coord (5' intronic start coordinate, for saturation)</br>
Column 6: exon_coord (exon start coordinate, for saturation)</br>
Column 7: 3intron_coord (3' intronic start coordinate, for saturation)</br>
Column 8: 5intron_action (snv/1del/2del/3del/NA)</br>
Column 9: exon_action (snv/1del/2del/3del/inframe/snvfull/snvre/NA)</br>
Column 10: 3intron_action (snv/1del/2del/3del/NA)</br>
Column 11: inframe (inframe CDS sequence)</br>
Column 12: 5base.no (additinal base at 5' end for adjusting inframe CDS, can be 0,1,2)</br>
Column 13: 3base.no (additinal base at 3' end for adjusting inframe CDS, can be 0,1,2)</br>
Column 14: 5constant (5' constant sequence, non-SGE and primer binding sequence)</br>
Column 15: 3constant (3' constant sequence, non-SGE and primer binding sequence)</br>
Column 16: 5adaptor (adaptor sequence for oligo pool amplification, eg: P5/P7 sequence. Cannot be empty, remove the last base of 5constant and input here if adaptor not needed)</br>
Coulmn 17: 3adaptor (adaptor sequence for oligo pool amplification, eg: P5/P7 sequence. Cannot be empty, remove the last base of 3constant and input here if adaptor not needed)</br>

The first step of my script will generate a working file (testchip) from this input tab file. The working file has 22 columns, you may see $1,$2,...$22 throughout the script.

**Working file format:**
1) Name(letter only. no . , _)	
2) 5intron_sequence	
3) exon_sequence	
4) 3intron_sequence	
5) 5intron_length
6) exon_length
7) 3intron_length
8) 5intron_coord	
9) exon_coord	
10) 3intron_coord	
11) 5intron_operation (snv,1del,2del,3del,NA)	
12) exon_operation (snv,1del,2del,3del,inframe,snvre,snvfull,NA)snvfull: all possible codon; snvre: most frequent redundant codon; junction problems considered
13) 3intron_operation (snv,1del,2del,3del,NA)	
14) inframe_sequece	
15) 5additional_base_number (0,1,2)	
16) 3additional_base_number (0,1,2)
17) 5additional_base_sequence
18) 3additional_base_sequence	
19) 5constant_region for specific exon amplification (must have at least 1 letter)
20) 3constant_region for specific exon amplification (must have at least 1 letter)
21) 5adaptor_sequence	(must have at least 1 letter)
22) 3adaptor_sequence (must have at least 1 letter)


**Codon frequency table that I used:**

Comma separated: 1) Codon, 2) Amino acid, 3)Frequency in human mRNA, 4) Codon frequency rank for each amino acid, T=top, U=unique

ATA,I,0.16,RANK3</br>
ATC,I,0.48,RANKT</br>
ATT,I,0.36,RANK2</br>
ATG,M,1.00,RANKU</br>
ACA,T,0.28,RANK2</br>
ACC,T,0.36,RANKT</br>
ACG,T,0.12,RANK4</br>
ACT,T,0.24,RANK3</br>
AAC,N,0.54,RANKT</br>
AAT,N,0.46,RANK2</br>
AAA,K,0.42,RANK2</br>
AAG,K,0.58,RANKT</br>
AGC,S,0.24,RANKT</br>
AGT,S,0.15,RANK4</br>
TCA,S,0.15,RANK5</br>
TCC,S,0.22,RANK2</br>
TCG,S,0.06,RANK6</br>
TCT,S,0.18,RANK3</br>
AGA,R,0.20,RANK2</br>
AGG,R,0.20,RANK3</br>
CGA,R,0.11,RANK5</br>
CGC,R,0.19,RANK4</br>
CGG,R,0.21,RANKT</br>
CGT,R,0.08,RANK6</br>
CTA,L,0.07,RANK5</br>
CTC,L,0.20,RANK2</br>
CTG,L,0.41,RANKT</br>
CTT,L,0.13,RANK3</br>
TTA,L,0.07,RANK6</br>
TTG,L,0.13,RANK4</br>
CCA,P,0.27,RANK3</br>
CCC,P,0.33,RANKT</br>
CCG,P,0.11,RANK4</br>
CCT,P,0.28,RANK2</br>
CAC,H,0.59,RANKT</br>
CAT,H,0.41,RANK2</br>
CAA,Q,0.25,RANK2</br>
CAG,Q,0.75,RANKT</br>
GTA,V,0.11,RANK4</br>
GTC,V,0.24,RANK2</br>
GTG,V,0.47,RANKT</br>
GTT,V,0.18,RANK3</br>
GCA,A,0.23,RANK3</br>
GCC,A,0.40,RANKT</br>
GCG,A,0.11,RANK4</br>
GCT,A,0.26,RANK2</br>
GAC,D,0.54,RANKT</br>
GAT,D,0.46,RANK2</br>
GAA,E,0.42,RANK2</br>
GAG,E,0.58,RANKT</br>
GGA,G,0.25,RANK2</br>
GGC,G,0.34,RANKT</br>
GGG,G,0.25,RANK3</br>
GGT,G,0.16,RANK4</br>
TTC,F,0.55,RANKT</br>
TTT,F,0.45,RANK2</br>
TAC,Y,0.57,RANKT</br>
TAT,Y,0.43,RANK2</br>
TAA,STOP,0.28,RANK2</br>
TAG,STOP,0.20,RANK3</br>
TGA,STOP,0.52,RANKT</br>
TGC,C,0.55,RANKT</br>
TGT,C,0.45,RANK2</br>
TGG,W,1.00,RANKU</br>

