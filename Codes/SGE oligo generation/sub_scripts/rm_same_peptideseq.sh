#!/bin/bash

set -u

pname=$(basename $0)
[ $# -ne 4 ] && { printf "split saturated_mutagenesis_triplet output to synonymous (*.sym.tab) and non-synonymous (*.nonsym.tab) sequence tab file. Also, input files will be translated and output as *.protein.tab.
Dependencies: translate.sh, transeq (EMBOSS), fasta2tab.pl, RemoveInList, and RetainInList
Usage: $pname [full path for input fasta file] [full path for saturated_mutagenesis_triplet.sh output file. tab file] [output file name prefix] [FRAME]

[full path for input fasta file]	fasta only.
[full path for saturated_mutagenesis_triplet.sh output file. tab file]
[outputfile name prefix]		STRING without space; use _ for space.Tab file only.
[FRAME]					Frame 1,2,3\n
"; exit 1; }

[ ! -f $1 ] && echo "$1 does not exist" && exit 1
[ ! -f $2 ] && echo "$2 does not exist" && exit 1


translate.sh $1 1 $1.protein.tab;

cat $2| awk '{print $1"\n"$2}' - > mut.fa;
translate.sh mut.fa $4 $2.protein.tab;

cat $2.protein.tab | RemoveInList 2 $1.protein.tab 2 > notsameprotein
cat $2 | RetainInList 1 notsameprotein 1 > $3.nonsym.tab

cat $2.protein.tab | RetainInList 2 $1.protein.tab 2 > sameprotein
cat $2 | RetainInList 1 sameprotein 1 >  $3.sym.tab

rm mut.fa notsameprotein sameprotein

