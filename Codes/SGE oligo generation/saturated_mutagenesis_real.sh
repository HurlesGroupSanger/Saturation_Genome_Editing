#!/bin/bash

set -u

pname=$(basename $0)
[ $# -ne 2 ] && { printf "Mutate each base per input DNA/RNA sequence(s).
Usage: $pname [full path for input fasta file] [core]

[full path for input tab file]tab only.
[core]		core
Dependencies:
1) saturated_mutagenesis_subreal.sh
	Dependencies:
	1) rm_same_peptideseq.sh
		Dependencies: translate.sh -> transeq (EMBOSS), fasta2tab.pl	
	2) RetainInList, RemoveInList, InsertInList
2) gnu-parallel
3) gnu-(awk, sed, uniq, sort)\n
"; exit 1; }

[ ! -f $1 ] && echo "$1 does not exist" && exit 1

#force nonsense thread to 1
if [ $2 -lt 1 ] || [ $2 -gt 50 ]; then
thread=1
else
thread=$2
fi

#Generate working file
# 1) Name(letter only. no . , _)	
# 2) 5intron_sequence	
# 3) exon_sequence	
# 4) 3intron_sequence	
# 5) 5intron_length
# 6) exon_length
# 7) 3intron_length
# 8) 5intron_coord	
# 9) exon_coord	
# 10) 3intron_coord	
# 11) 5intron_operation (snv,1del,2del,3del,NA)	
# 12) exon_operation (snv,1del,2del,3del,inframe,snvre,snvfull,NA)snvfull: all possible codon; snvre: most frequent redundant codon; junction problems considered
# 13) 3intron_operation (snv,1del,2del,3del,NA)	
# 14) inframe_sequece	
# 15) 5additional_base_number (0,1,2)	
# 16) 3additional_base_number (0,1,2)
# 17) 5additional_base_sequence
# 18) 3additional_base_sequence	
# 19) 5constant_region for specific exon amplification (must have at least 1 letter)
# 20) 3constant_region for specific exon amplification (must have at least 1 letter)
# 21) 5adaptor_sequence	(must have at least 1 letter)
# 22) 3adaptor_sequence (must have at least 1 letter)

sed 1d $1 | awk 'BEGIN { OFS = "\t" } ; {print $1,$2,$3,$4,length($2),length($3),length($4),$5,$6,$7,$8,$9,$10,$11,$12,$13,substr($11,1,$12), substr($11,length($11)-$13+1,length($11)),$14,$15,$16,$17}' | sed -E 's/\t\t/\tNA\t/g'| sed -E 's/\t\t/\tNA\t/g' > testchip

#seq_length QC; column 2 and 3 must equal
cat testchip | awk '{print $1"\t"length($3)+$15+$16"\t"length($14)}' > input_seq_length_qc 

#Generate folder based on the column 1 "Name"
parallel -j $thread -a 'testchip' --colsep '\t' "[ ! -d \"./{1}\" ] && mkdir {1}"

#5intron; saturated_mutagenesis_subreal.sh [Name] [input sequence] [sequence length] [type] [output file name] [pos] [inframe seq] [5 addition base] [3 addition base] [5 addition base seq] [3 addition base seq]
time parallel -j $thread -a 'testchip' --colsep '\t' "cd {1}; if [ \"{11}\" != \"NA\" ]; then saturated_mutagenesis_subreal.sh {1}.5intron {2} {5} {11} {1}_5intron_{11} {8} {14} {15} {16} {17} {18}; sed -i \"s/$/$(echo {3}{4}{20}{22})/g\" {1}_5intron_{11}; sed -i \"s/\t/\t$(echo {21}{19})/g\" {1}_5intron_{11}; fi; cd .." ;

#3intron
time parallel -j $thread -a 'testchip' --colsep '\t' "cd {1}; if [ \"{13}\" != \"NA\" ]; then saturated_mutagenesis_subreal.sh {1}.3intron {4} {7} {13} {1}_3intron_{13} {10} {14} {15} {16} {17} {18}; sed -i \"s/\t/\t$(echo {21}{19}{2}{3})/g\" {1}_3intron_{13}; sed -i \"s/$/$(echo {20}{22})/g\" {1}_3intron_{13}; fi; cd .." ;

#exon
time parallel -j $thread -a 'testchip' --colsep '\t' "cd {1}; if [ \"{12}\" != \"NA\" ]; then saturated_mutagenesis_subreal.sh {1}.exon {3} {6} {12} {1}_exon_{12} {9} {14} {15} {16} {17} {18}; sed -i \"s/\t/\t$(echo {21}{19}{2})/g\" {1}_exon_{12}; sed -i \"s/$/$(echo {4}{20}{22})/g\" {1}_exon_{12}; fi; cd .." ;

#generate to_buy

parallel -j $thread -a 'testchip' --colsep '\t' "cd {1}; cat *_snv *_1del *_2del *_3del *_inframe *_snvre *_snvfull > {1}.tobuy; (sort -k 2,2 {1}.tobuy | uniq -f 1 | wc -l; wc -l {1}.tobuy; wc -l *.info; wc -l *intron*; wc -l *inframe) > {1}.qc; cd .."





