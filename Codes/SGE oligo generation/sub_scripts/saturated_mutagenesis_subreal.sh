#!/bin/bash

set -u

pname=$(basename $0)
[ $# -ne 11 ] && { printf "Mutate each base per input DNA/RNA sequence(s).
Usage: $pname [Seq Name] [input sequence] [sequence length] [type] [output file name] [pos] [inframe seq] [5 addition base] [3 addition base] [5 addition base seq] [3 addition base seq]
\n
"; exit 1; }

#[ -f $5 ] && echo "$5 exist. Change output file name" && exit 1

if [ $3 -lt 1 ]; then

printf "sequence length must be positive integer\n
Usage: $pname [Seq Name] [input sequence] [sequence length] [type] [output file name] [pos] [inframe seq] [5 addition base] [3 addition base] [5 addition base seq] [3 addition base seq]\n\n"

exit 1

fi
# $1=seqname (sequence name preffix)
# $2=input sequence (for any type of actions)
# $3=sequence length
# $4=type (action type, snv/1del/2del/3del/inframe/snvfull/snvre)
# $5=output file name (no space use underscore "_" for whitespace)
# $6=pos (start coordinate)
# $7=inframe seq (inframe CDS sequence for that exon)
# $8=5 addition base number (for adjusting frame)
# $9=3 addition base number (for adjusting frame)
# $10=5 addition base (for retaining the same bases in upstream exon)
# $11=3 addition base (for for retaining the same bases in downstream exon)

if [ "$4" = "NA" ]; then

exit 1

fi

###

if [ "$4" = "snv" ]; then
echo $2 $1> $1.snv.tmp;
pos=$6-1
for i in $(seq 1 $3);
do ((pos++));
	for n in A T C G; 
	do sed "s/./$n/$i" $1.snv.tmp | awk '{print ">"$2"_Change_position_'$pos'_to_'$n'\t" $1}'; 
	done; 
done | grep -v $2 > $5;
rm $1.snv.tmp

fi

###

if [ "$4" = "1del" ]; then
echo $2 $1> $1.1del.tmp;
pos=$6-1
for i in $(seq 1 $3);
do ((pos++));
sed "s/.//$i" $1.1del.tmp | awk '{print ">"$2"_onedelete_position_'$pos'\t" $1}'; 
done | sort -u -k 2,2 > $5;
rm $1.1del.tmp

fi

###

if [ "$4" = "2del" ]; then
echo $2 $1> $1.2del.tmp;
new=$(echo $3/2 | bc)
pos=$6
for i in $(seq 1 $new);
do sed "s/..//$i" $1.2del.tmp | awk '{print ">"$2"_twodelete_position_'$pos'\t" $1}'; pos=$(echo $pos+2| bc);
done | sort -u -k 2,2 > $5;
rm $1.2del.tmp

fi

###

if [ "$4" = "3del" ]; then
echo $2 $1> $1.3del.tmp;
new=$(echo $3/3 | bc)
pos=$6
for i in $(seq 1 $new);
do sed "s/...//$i" $1.3del.tmp | awk '{print ">"$2"_threedelete_position_'$pos'\t" $1}'; pos=$(echo $pos+3| bc);
done | sort -u -k 2,2 > $5;
rm $1.3del.tmp

fi

###

if [ "$4" = "inframe" ]; then
echo $7 $1> $1.inframe.tmp;
new=$(echo \($3 + $8 + $9\)/3 | bc)
new2=$(echo \(\($3 + $8 + $9\)/3\) -1 | bc)
pos=$(echo $6-$8+3 | bc)
pos2=$6
if [ "$9" = "1" ]; then
pos3=$(echo $6+$3-2 | bc)
elif [ "$9" = "0" ]; then
pos3=$(echo $6+$3-3 | bc)
elif [ "$9" = "2" ]; then
pos3=$(echo $6+$3-1 | bc)
fi

if [ "$8" = "0" ]; then
sed "s/...//1" $1.inframe.tmp | awk '{print ">"$2"_inframedelete_position_'$pos2'\t" $1}' | sed -E "s/.{$9}$//g" >> $5.tmp
fi

for i in $(seq 2 $new2);
do sed "s/...//$i" $1.inframe.tmp | awk '{print ">"$2"_inframedelete_position_'$pos'\t" $1}'; pos=$(echo $pos+3| bc);
done | sed -E "s/\t.{$8}/\t/g" | sed -E "s/.{$9}$//g" >> $5.tmp

if [ "$9" = "0" ]; then
sed "s/...//$new" $1.inframe.tmp | awk '{print ">"$2"_inframedelete_position_'$pos3'\t" $1}' |  sed -E "s/\t.{$8}/\t/g" >> $5.tmp
fi

sort -u -k 2,2 $5.tmp > $5;
rm $1.inframe.tmp $5.tmp

fi

###

if [ "$4" = "snvfull" ]; then
echo $7 $1> $1.snvfull.tmp;
printf "
ATA,I,0.16,RANK3
ATC,I,0.48,RANKT
ATT,I,0.36,RANK2
ATG,M,1.00,RANKU
ACA,T,0.28,RANK2
ACC,T,0.36,RANKT
ACG,T,0.12,RANK4
ACT,T,0.24,RANK3
AAC,N,0.54,RANKT
AAT,N,0.46,RANK2
AAA,K,0.42,RANK2
AAG,K,0.58,RANKT
AGC,S,0.24,RANKT
AGT,S,0.15,RANK4
TCA,S,0.15,RANK5
TCC,S,0.22,RANK2
TCG,S,0.06,RANK6
TCT,S,0.18,RANK3
AGA,R,0.20,RANK2
AGG,R,0.20,RANK3
CGA,R,0.11,RANK5
CGC,R,0.19,RANK4
CGG,R,0.21,RANKT
CGT,R,0.08,RANK6
CTA,L,0.07,RANK5
CTC,L,0.20,RANK2
CTG,L,0.41,RANKT
CTT,L,0.13,RANK3
TTA,L,0.07,RANK6
TTG,L,0.13,RANK4
CCA,P,0.27,RANK3
CCC,P,0.33,RANKT
CCG,P,0.11,RANK4
CCT,P,0.28,RANK2
CAC,H,0.59,RANKT
CAT,H,0.41,RANK2
CAA,Q,0.25,RANK2
CAG,Q,0.75,RANKT
GTA,V,0.11,RANK4
GTC,V,0.24,RANK2
GTG,V,0.47,RANKT
GTT,V,0.18,RANK3
GCA,A,0.23,RANK3
GCC,A,0.40,RANKT
GCG,A,0.11,RANK4
GCT,A,0.26,RANK2
GAC,D,0.54,RANKT
GAT,D,0.46,RANK2
GAA,E,0.42,RANK2
GAG,E,0.58,RANKT
GGA,G,0.25,RANK2
GGC,G,0.34,RANKT
GGG,G,0.25,RANK3
GGT,G,0.16,RANK4
TTC,F,0.55,RANKT
TTT,F,0.45,RANK2
TAC,Y,0.57,RANKT
TAT,Y,0.43,RANK2
TAA,STOP,0.28,RANK2
TAG,STOP,0.20,RANK3
TGA,STOP,0.52,RANKT
TGC,C,0.55,RANKT
TGT,C,0.45,RANK2
TGG,W,1.00,RANKU
" > codon_table

new=$(echo \($3 + $8 + $9\)/3 | bc)
new2=$(echo \(\($3 + $8 + $9\)/3\) -1 | bc)
pos=$(echo $6-$8+3 | bc)
pos2=$6
if [ "$9" = "1" ]; then
pos3=$(echo $6+$3-2 | bc)
elif [ "$9" = "0" ]; then
pos3=$(echo $6+$3-3 | bc)
elif [ "$9" = "2" ]; then
pos3=$(echo $6+$3-1 | bc)
fi

if [ "$8" = "0" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/1" $1.snvfull.tmp | awk '{print ">"$2"_Change_codon_position_'$pos2'_to_'$aa'\t" $1}' ; done | sed -E "s/.{$9}$//g" >> $5.tmp

elif [ "$8" = "1" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/1" $1.snvfull.tmp | grep -i ^${10} | sed 's/^.//g' | awk '{print ">"$2"_Change_codon_position_'$pos2'_to_'$aa'\t" $1}' ; done | sed -E "s/.{$9}$//g" >> $5.tmp

elif [ "$8" = "2" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/1" $1.snvfull.tmp | grep -i ^${10} | sed 's/^..//g' | awk '{print ">"$2"_Change_codon_position_'$pos2'_to_'$aa'\t" $1}' ; done | sed -E "s/.{$9}$//g" >> $5.tmp
fi


for i in $(seq 2 $new2); 
do for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$i" $1.snvfull.tmp | awk '{print ">"$2"_Change_codon_position_'$pos'_to_'$aa'\t" $1}'; 
done; pos=$(echo $pos+3| bc);
done | sed -E "s/\t.{$8}/\t/g" | sed -E "s/.{$9}$//g" >> $5.tmp; 


if [ "$9" = "0" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$new" $1.snvfull.tmp | awk '{print ">"$2"_Change_codon_position_'$pos3'_to_'$aa'\t" $1}' ; done | sed -E "s/\t.{$8}/\t/g" >> $5.tmp

elif [ "$9" = "1" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$new" $1.snvfull.tmp | awk '{print ">"$2"_Change_codon_position_'$pos3'_to_'$aa'\t" $1}' | grep -i ${11}$ | sed 's/.$//g'  ; done | sed -E "s/\t.{$8}/\t/g" >> $5.tmp

elif [ "$9" = "2" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$new" $1.snvfull.tmp | awk '{print ">"$2"_Change_codon_position_'$pos3'_to_'$aa'\t" $1}' | grep -i ${11}$ | sed 's/..$//g' ; done | sed -E "s/\t.{$8}/\t/g" >> $5.tmp

fi

grep -v $2 $5.tmp > $5

rm $1.snvfull.tmp codon_table $5.tmp

fi

###

if [ "$4" = "snvre" ]; then
echo $7 $1> $1.snvre.tmp;
printf "
ATA,I,0.16,RANK3
ATC,I,0.48,RANKT
ATT,I,0.36,RANK2
ATG,M,1.00,RANKU
ACA,T,0.28,RANK2
ACC,T,0.36,RANKT
ACG,T,0.12,RANK4
ACT,T,0.24,RANK3
AAC,N,0.54,RANKT
AAT,N,0.46,RANK2
AAA,K,0.42,RANK2
AAG,K,0.58,RANKT
AGC,S,0.24,RANKT
AGT,S,0.15,RANK4
TCA,S,0.15,RANK5
TCC,S,0.22,RANK2
TCG,S,0.06,RANK6
TCT,S,0.18,RANK3
AGA,R,0.20,RANK2
AGG,R,0.20,RANK3
CGA,R,0.11,RANK5
CGC,R,0.19,RANK4
CGG,R,0.21,RANKT
CGT,R,0.08,RANK6
CTA,L,0.07,RANK5
CTC,L,0.20,RANK2
CTG,L,0.41,RANKT
CTT,L,0.13,RANK3
TTA,L,0.07,RANK6
TTG,L,0.13,RANK4
CCA,P,0.27,RANK3
CCC,P,0.33,RANKT
CCG,P,0.11,RANK4
CCT,P,0.28,RANK2
CAC,H,0.59,RANKT
CAT,H,0.41,RANK2
CAA,Q,0.25,RANK2
CAG,Q,0.75,RANKT
GTA,V,0.11,RANK4
GTC,V,0.24,RANK2
GTG,V,0.47,RANKT
GTT,V,0.18,RANK3
GCA,A,0.23,RANK3
GCC,A,0.40,RANKT
GCG,A,0.11,RANK4
GCT,A,0.26,RANK2
GAC,D,0.54,RANKT
GAT,D,0.46,RANK2
GAA,E,0.42,RANK2
GAG,E,0.58,RANKT
GGA,G,0.25,RANK2
GGC,G,0.34,RANKT
GGG,G,0.25,RANK3
GGT,G,0.16,RANK4
TTC,F,0.55,RANKT
TTT,F,0.45,RANK2
TAC,Y,0.57,RANKT
TAT,Y,0.43,RANK2
TAA,STOP,0.28,RANK2
TAG,STOP,0.20,RANK3
TGA,STOP,0.52,RANKT
TGC,C,0.55,RANKT
TGT,C,0.45,RANK2
TGG,W,1.00,RANKU
" > codon_table

new=$(echo \($3 + $8 + $9\)/3 | bc)
new2=$(echo \(\($3 + $8 + $9\)/3\) -1 | bc)
pos=$(echo $6-$8+3 | bc)
pos2=$6
if [ "$9" = "1" ]; then
pos3=$(echo $6+$3-2 | bc)
elif [ "$9" = "0" ]; then
pos3=$(echo $6+$3-3 | bc)
elif [ "$9" = "2" ]; then
pos3=$(echo $6+$3-1 | bc)
fi

#Dealing with the 5' junction condon

if [ "$8" = "0" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/1" $1.snvre.tmp | awk '{print ">"$2"_Change_codon_position_'$pos2'_to_'$aa'\t" $1}' ; done > $5.tmptmp1 
cat $5.tmptmp1 | sed -E "s/.{$9}$//g" > $5.tmp1

elif [ "$8" = "1" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/1" $1.snvre.tmp | grep -i ^${10} | awk '{print ">"$2"_Change_codon_position_'$pos2'_to_'$aa'\t" $1}' ; done > $5.tmptmp1
cat $5.tmptmp1 | sed -E 's/\t./\t/g' |sed -E "s/.{$9}$//g" > $5.tmp1

elif [ "$8" = "2" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/1" $1.snvre.tmp | grep -i ^${10} | awk '{print ">"$2"_Change_codon_position_'$pos2'_to_'$aa'\t" $1}' ; done > $5.tmptmp1
cat $5.tmptmp1 | sed -E 's/\t../\t/g' | sed -E "s/.{$9}$//g" > $5.tmp1
fi

#Dealing with the sequence starting with codon 2 to codon n-1. ie, not the junction codon
for i in $(seq 2 $new2); 
do for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$i" $1.snvre.tmp | awk '{print ">"$2"_Change_codon_position_'$pos'_to_'$aa'\t" $1}'; 
done; pos=$(echo $pos+3| bc);
done > $5.tmptmp2
cat $5.tmptmp2 | sed -E "s/\t.{$8}/\t/g" | sed -E "s/.{$9}$//g" > $5.tmp2; 

#Dealing with the 3' junction codon
if [ "$9" = "0" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$new" $1.snvre.tmp | awk '{print ">"$2"_Change_codon_position_'$pos3'_to_'$aa'\t" $1}' ; done > $5.tmptmp3
cat $5.tmptmp3 | sed -E "s/\t.{$8}/\t/g" > $5.tmp3

elif [ "$9" = "1" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$new" $1.snvre.tmp | awk '{print ">"$2"_Change_codon_position_'$pos3'_to_'$aa'\t" $1}' | grep -i ${11}$ ; done > $5.tmptmp3
cat $5.tmptmp3 | sed 's/.$//g' | sed -E "s/\t.{$8}/\t/g" > $5.tmp3

elif [ "$9" = "2" ]; then
for n in ATA ATC ATT ATG ACA ACC ACG ACT AAC AAT AAA AAG AGC AGT TCA TCC TCG TCT AGA AGG CGA CGC CGG CGT CTA CTC CTG CTT TTA TTG CCA CCC CCG CCT CAC CAT CAA CAG GTA GTC GTG GTT GCA GCC GCG GCT GAC GAT GAA GAG GGA GGC GGG GGT TTC TTT TAC TAT TAA TAG TGA TGC TGT TGG; 
do aa=$(grep $n codon_table);
sed "s/.../$n/$new" $1.snvre.tmp | awk '{print ">"$2"_Change_codon_position_'$pos3'_to_'$aa'\t" $1}' | grep -i ${11}$ ; done > $5.tmptmp3
cat $5.tmptmp3 | sed 's/..$//g' | sed -E "s/\t.{$8}/\t/g" > $5.tmp3

fi
#tmptmp is the file contain the sequence without removing the addition base in 5' and 3' end. tmp has the additional bases removed

#The snvfull generated based on inframe CDS sequence
cat $5.tmptmp1 $5.tmptmp2 $5.tmptmp3 > $5.tmptmp
grep -v $7 $5.tmptmp > $5.tmpo

#The snvfull generated based on exon sequence ie, without the additional bases for inframe CDS
cat $5.tmp1 $5.tmp2 $5.tmp3 > $5.tmp
grep -v $2 $5.tmp > $5.tmpo.real


rm $5.tmptmp1 $5.tmptmp2 $5.tmptmp3 $5.tmp1 $5.tmp2 $5.tmp3

#Split the synonymous ; two impotant files $5.sym.tab $5.nonsym.tab

cat $1.snvre.tmp | awk '{print ">"$2"\n"$1}' > $1.snvre.fa

rm_same_peptideseq.sh $1.snvre.fa $5.tmpo $5 1

#get triplet info (inframe CDS based)

cat $5.tmpo |awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"\t"$7 }'| awk 'BEGIN { FS = "," } ; { print $0"\t"$5"\t"$7 }'| cut -f 1,2,3,6,7 | awk '{print $1"\t"$2"\t"$3"_"$4"_"$5}'> $5.tripletfull.info

#get triplet_real info (exon based)

cat $5.tmpo.real |awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"\t"$7 }'| awk 'BEGIN { FS = "," } ; { print $0"\t"$5"\t"$7 }'| cut -f 1,2,3,6,7 | awk '{print $1"\t"$2"\t"$3"_"$4"_"$5}'> $5.tripletreal.info

#get original info

grep $7 $5.tmptmp | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"\t"$7 }'| awk 'BEGIN { FS = "," } ; { print $0"\t"$5"\t"$7 }'| cut -f 1,2,3,6,7 | awk '{print $1"\t"$2"\t"$3"_"$4"_"$5}'> $5.original.info

#cat $5.tmpinfo | awk '{print $1"\t"$2"\t"$3"_"$4"\t"$5"\t"$3"_"$4"_"$5}' > $5.tmpinfo2

#generate snv sequence $5.tmplist

if [ "$8" = "0" ]; then
echo $7 $1> $1.snvre.tmp2;
pos=$6-1
to=$(echo $3+0 |bc)
for i in $(seq 1 $to);
do ((pos++));
	for n in A T C G; 
	do sed "s/./$n/$i" $1.snvre.tmp2 | awk '{print "'$pos'_to_'$n'\t" $1}'; 
	done; 
done | grep -v $2 > $5.tmplist;
rm $1.snvre.tmp2

elif [ "$8" = "1" ]; then
echo $7 $1> $1.snvre.tmp2;
pos=$6-1
to=$(echo $3+1 |bc)
for i in $(seq 2 $to);
do ((pos++));
	for n in A T C G; 
	do sed "s/./$n/$i" $1.snvre.tmp2 | awk '{print "'$pos'_to_'$n'\t" $1}'; 
	done; 
done | grep -v $2 > $5.tmplist;
rm $1.snvre.tmp2

elif [ "$8" = "2" ]; then
echo $7 $1> $1.snvre.tmp2;
pos=$6-1
to=$(echo $3+2 |bc)
for i in $(seq 3 $to);
do ((pos++));
	for n in A T C G; 
	do sed "s/./$n/$i" $1.snvre.tmp2 | awk '{print "'$pos'_to_'$n'\t" $1}'; 
	done; 
done | grep -v $2 > $5.tmplist;
rm $1.snvre.tmp2

fi

#get the SNV which is not synonymous

cat $5.nonsym.tab | RetainInList 2 $5.tmplist 2 | InsertInList 2 $5.tmplist 2 1 | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5 }' | cut -f 1-4 | awk 'BEGIN { FS = "," } ; { print $0"\t"$2"\t"$4 }' | cut -f 1-6 | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$5"_"$6"\tmissnv"}' |  sort -n -k 3,3n > $5.nosymsnv.info #| awk '{print $1"_snv\t"$2"\t"$3"\t"$4"\t"$5}' 

#get the SNV which is synonymous

cat $5.sym.tab | RetainInList 2 $5.tmplist 2 | InsertInList 2 $5.tmplist 2 1 | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5 }' | cut -f 1-4 | awk 'BEGIN { FS = "," } ; { print $0"\t"$2"\t"$4 }' | cut -f 1-6 | awk '{print $1"\t"$2"\t"$3"\t"$4"_"$5"_"$6"\tsynsnv"}' |  sort -n -k 3,3n > $5.symsnv.info

#get the synonymous which is not SNV

cat $5.sym.tab | RemoveInList 2 $5.tmplist 2 | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"\t"$7 }'| awk 'BEGIN { FS = "," } ; { print $0"\t"$5"\t"$7 }'| cut -f 1,2,3,6,7 | awk '{print $1"\t"$2"\t"$3"_"$4"_"$5}' | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"_"$6"_"$7 }' | cut -d ',' -f 1-4 | awk '{print $1"\t"$2"\t"$4"\t"$3"\tsynredundant"}' > $5.symnotsnv.info

#get the redundant from not synonymous

cat $5.nosymsnv.info | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"_"$10 }' | sort -k 6,6 | uniq -u -f 5 | grep -v RANKU > $5.nosymsnv.uniq
cat $5.nosymsnv.uniq | cut -f 4 | sed 's/RANK2/RANK1/g' | sed 's/RANK3/RANK1/g' | sed 's/RANK4/RANK1/g' | sed 's/RANK5/RANK1/g' | sed 's/RANK6/RANK1/g'|sed 's/RANKT/RANK2/g' | sed 's/RANK1/RANKT/g' > $5.nosymsnv.re
cat $5.tripletfull.info | RetainInList 3 $5.nosymsnv.re 1 | awk 'BEGIN { FS = "_" } ; { print $0"\t"$5"_"$6"_"$7 }' | cut -d ',' -f 1-4 | awk '{print $1"\t"$2"\t"$4"\t"$3"\tmisredundant"}' > $5.nosymredundant.info


#.info contains the annotation information. example:
# >DDX3XE11sg1.exon_Change_codon_position_41345180_to_AAG,K,0.58,RANKT	GTACTTGGTGTTAGATGAAGCTGATCGGATGTTGGACATGGGATTTGAGCCTCAGATTCGTAGAATAGTCGAACAAGATACTATGCCTCCAAAGGGTGTCCGCCACACTATGATGTTTAGTGCTACTTTTCCTAAGGAAATACAG	41345180_to_G	41345180_K_RANKT	synsnv
# >DDX3XE11sg1.exon_Change_codon_position_41345196_to_GCC,A,0.40,RANKT	ATACTTGGTGTTAGATGCCGCTGATCGGATGTTGGACATGGGATTTGAGCCTCAGATTCGTAGAATAGTCGAACAAGATACTATGCCTCCAAAGGGTGTCCGCCACACTATGATGTTTAGTGCTACTTTTCCTAAGGAAATACAG	41345196_to_GCC	41345196_A_RANKT	misredundant
# >DDX3XE11sg1.exon_Change_codon_position_41345322_to_CGG,R,0.21,RANKT	ATACTTGGTGTTAGATGAAGCTGATCGGATGTTGGACATGGGATTTGAGCCTCAGATTCGTAGAATAGTCGAACAAGATACTATGCCTCCAAAGGGTGTCCGCCACACTATGATGTTTAGTGCTACTTTTCCTAAGGAAATACGG	41345323_to_G	41345322_R_RANKT	missnv

# Cat and Replace the inframe sequence with real exon sequence
cat $5.symsnv.info $5.nosymsnv.info $5.symnotsnv.info $5.nosymredundant.info | InsertInList 4 $5.tripletreal.info 3 2 N/A| awk '{print $1"\t"$6"\t"$3"\t"$4"\t"$5}' > $5.final.info

sort -k 3,3 $5.final.info |  cut -f 1,2 > $5

rm $5.nosymsnv.uniq $5.nosymsnv.re

rm $1.snvre.tmp codon_table *.tab *.fa $5.tmpo $5.tmp $5.tmptmp $5.tmplist $5.tmpo.real $5.tripletreal.info

fi

#$pname [Name] [input sequence] [sequence length] [type] [output file name] [pos] [inframe seq] [5 addition base] [3 addition base]
#Alternative for more flexible replacement
#echo "AATTAACCCGT"| sed -E 's/(.{3}).{3}/\1x/'

