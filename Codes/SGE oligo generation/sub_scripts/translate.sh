#!/bin/bash

pname=$(basename $0)
[ $# -ne 3 ] && { printf "translate fasta file and output as tab file.
Dependencies: transeq (EMBOSS), fasta2tab.pl
Usage: $pname [full path for input fasta file] [FRAME] [output file name]

[full path for input fasta file]	fasta only.
[FRAME]					Refer to transeq manual. Frame number to use
[outputfile name]			STRING without space; use _ for space.Tab file only.\n
"; exit 1; }

INPUTFILE=$1
FRAME=$2
TEMPFILE=$(echo $1.temp)
TEMPFILE2=$(echo $1.temp2)
TEMPFILE3=$(echo $1.temp3)
OUTPUTFILE=$3

fasta2tab.pl $1 | awk '{print ">"$1}'> "$TEMPFILE2"

# Translate sequence
transeq -sequence "$INPUTFILE" -outseq "$TEMPFILE" -frame "$FRAME"

fasta2tab.pl "$TEMPFILE" | cut -f 2 > "$TEMPFILE3"

paste "$TEMPFILE2" "$TEMPFILE3" > "$OUTPUTFILE"

# Remove temp file
rm "$TEMPFILE" "$TEMPFILE2" "$TEMPFILE3"
