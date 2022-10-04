#!/bin/bash
#BSUB -P SGE
#BSUB -G ddd-grp
#BSUB -o log/plasmid_tagdust_R_%J_%I.LSFout.txt -e log/plasmid_tagdust_R_%J_%I.LSFerr.txt
#BSUB -n 15
#BSUB -R "span[hosts=1] select[mem>2200] rusage[mem=2200]" -M 2200
#BSUB -q normal
#BSUB -J "plasmid_tagdust_R[1-102]"

id=$LSB_JOBINDEX


#column1
name=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 1)
#column2
directory=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 2)
#column3
exon=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 3)
#column4
constant5=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 4)
#column5
constant3=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 5)
#column6
minlen=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 6)
#column7
ref=$(sed 1d plasmid_all_tagdust_input_tsv.txt | head -n $id |tail -n 1| cut -f 7)
#column8 comments will not be used

[ ! -d "$directory" ] && echo "$directory is not found. Exits with errors." 1>&2 && exit 1
[ ! -d "$directory/mapping" ] && mkdir -p "$directory/mapping" && echo "$directory/mapping is not found. Directory is created." 1>&2
[ ! -d "$directory/mapping/contig" ] && mkdir -p "$directory/mapping/contig" && echo "$directory/mapping/contig is not found. Directory is created." 1>&2
[ ! -d "$directory/mapping/contig/raw" ] && mkdir -p "$directory/mapping/contig/raw" && echo "$directory/mapping/contig/raw is not found. Directory is created." 1>&2
[ ! -d "$directory/mapping/contig/full" ] && mkdir -p "$directory/mapping/contig/full" && echo "$directory/mapping/contig/full is not found. Directory is created." 1>&2

#Step 1: Tagdust

[ ! -e "$directory/${name}/${name}_merge.fastq.gz" ] && echo "Step 1 Tagdust: The input file $directory/${name}/${name}_merge.fastq.gz is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/${name}_tagdust_${exon}.fq.gz" ]
then echo "Step 1: Tagdust is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

tagdust -t 15 -1 S:$constant5 -2 R:N -3 S:$constant3 -minlen $minlen -o $directory/${name}_tagdust_$exon $directory/${name}/${name}_merge.fastq.gz 1> $directory/${name}_tagdust_$exon.stdout.txt 2> $directory/${name}_tagdust_$exon.stderr.txt

parallel -j 2 "gzip {}" ::: $directory/${name}_tagdust_${exon}.fq $directory/${name}_tagdust_${exon}_un.fq

echo "Step 1: Tagdust is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 1 Tagdust: $directory/${name}_tagdust_${exon}.fq.gz exists. File is not overwritten. Skipped Step 1." 1>&2
fi

#Step 2: Mapping

[ ! -e "$directory/${name}_tagdust_${exon}.fq.gz" ] && echo "Step 2 Mapping: $directory/${name}_tagdust_${exon}.fq.gz is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/mapping/contig/full/${name}_contigfull_${exon}.txt" ]
then echo "Step 2: Mapping is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

zcat "$directory/${name}_tagdust_${exon}.fq.gz" | paste - - - - | cut -f 1,2 | RetainInList 2 "$ref" 2 | InsertInList 2 "$ref" 2 1 | gzip 1>"$directory/mapping/${name}_tagdust_${exon}_map.gz" 2> "$directory/mapping/${name}_tagdust_${exon}_mapping.stderr.txt"

zcat "$directory/${name}_tagdust_${exon}.fq.gz" | paste - - - - | cut -f 1,2 | RemoveInList 2 "$ref" 2 | gzip 1>"$directory/mapping/${name}_tagdust_${exon}_unmap.gz" 2>> "$directory/mapping/${name}_tagdust_${exon}_mapping.stderr.txt"

zcat "$directory/mapping/${name}_tagdust_${exon}_map.gz" | cut -f 2,3 | sort -k 2,2 | uniq -f 1 -c | sed 's/^\s*//g'| tr ' ' '\t' 1> "$directory/mapping/contig/raw/${name}_contig_${exon}.txt" 2>> "$directory/mapping/${name}_tagdust_${exon}_mapping.stderr.txt"

cat "$ref" | InsertInList 2 "$directory/mapping/contig/raw/${name}_contig_${exon}.txt" 2 1 0 > "$directory/mapping/contig/full/${name}_contigfull_${exon}.txt" 2>> "$directory/mapping/${name}_tagdust_${exon}_mapping.stderr.txt"

echo "Step 2: Mapping is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 2 Mapping: $directory/mapping/contig/full/${name}_contigfull_${exon}.txt exists. File is not overwritten. Skipped Step 2." 1>&2
fi

#Step 3: R Mapping

[ ! -d "$directory/mapping/contig/full/R" ] && mkdir -p "$directory/mapping/contig/full/R" && echo "$directory/mapping/contig/full/R is not found. Directory is created." 1>&2

[ ! -e "$directory/${name}_tagdust_${exon}.fq.gz" ] && echo "Step 3 R Mapping: $directory/${name}_tagdust_${exon}.fq.gz is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/mapping/contig/full/R/${name}_${exon}_mapped.txt" ]
then echo "Step 3: R Mapping is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

zcat "$directory/${name}_tagdust_${exon}.fq.gz"| paste - - - - | cut -f 2 | awk '{ cnts[$0] += 1 } END { for (v in cnts) print cnts[v]"\t"v }' 1> "$directory/mapping/contig/full/R/${name}_${exon}_all_count.txt" 2> "$directory/mapping/${name}_tagdust_${exon}_Rmapping.stderr.txt"

printf "library(tidyverse)
raw_count <- read.csv(\"$directory/mapping/contig/full/R/${name}_${exon}_all_count.txt\",sep=\"\\\t\",header =F)
oligo <- read.csv(\"$ref\",sep=\"\\\t\",header =F)
mapped <- left_join(oligo, raw_count, by=\"V2\")
unmapped <- anti_join(raw_count, oligo, by=\"V2\") %%>%% arrange(desc(V1))
colnames(mapped) <- c(\"Name\",\"Sequence\",\"Count\")
colnames(unmapped) <- c(\"Count\",\"Sequence\")
write.table(mapped ,sep=\"\\\t\", file=\"$directory/mapping/contig/full/R/${name}_${exon}_mapped.txt\", row.names=F)
write.table(unmapped ,sep=\"\\\t\", file=\"$directory/mapping/contig/full/R/${name}_${exon}_unmapped.txt\", row.names=F)
rm(list=ls())
gc()
" > "$directory/mapping/contig/full/R/${name}_${exon}.r"

source activate R

Rscript "$directory/mapping/contig/full/R/${name}_${exon}.r" 2>> "$directory/mapping/${name}_tagdust_${exon}_Rmapping.stderr.txt"

conda deactivate

sed -i 's/"//g' $directory/mapping/contig/full/R/${name}_${exon}_mapped.txt
sed -i 's/"//g' $directory/mapping/contig/full/R/${name}_${exon}_unmapped.txt

echo "Step 3: R Mapping is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 3 R Mapping: $directory/mapping/contig/full/R/${name}_${exon}_mapped.txt exists. File is not overwritten. Skipped Step 3." 1>&2
fi

