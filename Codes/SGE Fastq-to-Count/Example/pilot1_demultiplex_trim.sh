#!/bin/bash
#BSUB -P SGE
#BSUB -G ddd-grp
#BSUB -o log/pilot1_%J_%I.LSFout.txt -e log/pilot1_%J_%I.LSFerr.txt
#BSUB -n 4
#BSUB -R "span[hosts=1] select[mem>300] rusage[mem=300]" -M 300
#BSUB -q normal
#BSUB -J "pilot1[1-180]"

id=$LSB_JOBINDEX

#column1
name="$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 1)"
#column2
directory="$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 2 )"
#column3
read1path="$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 3)"
#column4
#read2path="$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 4)"
#column5
exon="$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 5)"
#column6
constant5=$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 6 | rev |cut -c 1-10| rev)
#column7
constant3=$(sed 1d demultiplex.txt | head -n $id |tail -n 1| cut -f 7| cut -c 1-10)
#column8 comments will not be used


[ ! -d "$directory/demultiplex" ] && mkdir -p "$directory/demultiplex" && echo "$directory/demultiplex is not found. Directory is created." 1>&2
[ ! -d "$directory/trim" ] && mkdir -p "$directory/trim" && echo "$directory/trim is not found. Directory is created." 1>&2

#Step 1: Demultiplex

[ ! -e "$read1path" ] && echo "Step 1 Demultiplex: $read1path is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/demultiplex/${name}_${exon}_R1.fq.gz" ]
then echo "Step 1: Demultiplex is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

#Old method: awk, perfect match only. Illumina adaptor must be trimmed before demultiplex
#zcat "$directory/${name}/${name}_L001_R1_001_trimmed.fq.gz" |paste - - - - | awk -v pat="^${constant5}.*${constant3}$" -F '\t' '$2 ~ pat' | tr '\t' '\n' | gzip > $directory/demultiplex/${newname}_${exon}_R1.fq.gz

#New method: tre-agrep. Presence of Illumina adaptor is OK.
zcat "$read1path" |paste - - - - | /nfs/users/nfs_h/hk5/miniconda2/bin/agrep -1 "${constant5}.*${constant3}" | tr '\t' '\n' | gzip > $directory/demultiplex/${name}_${exon}_R1.fq.gz

echo "Step 1: Demultiplex is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 1  Demultiplex: $directory/demultiplex/${name}_${exon}_R1.fq.gz exists. File is not overwritten. Skipped Step 1." 1>&2
fi

#Step 2: Trim-galore

[ ! -e "$directory/demultiplex/${name}_${exon}_R1.fq.gz" ] && echo "Step 2: Trim-galore: $directory/demultiplex/${name}_${exon}_R1.fq.gz  is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/trim/${name}_${exon}_R1_trimmed.fq.gz" ]
then echo "Step 2: Trim-galore is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

source activate /nfs/users/nfs_h/hk5/miniconda2/envs/trim_galore
trim_galore --fastqc --length 0 -q 20 -o $directory/$name "$directory/demultiplex/${name}_${exon}_R1.fq.gz"  1>$directory/${name}_${exon}_trimgalore.stdout.txt 2>$directory/${name}_${exon}_trimgalore.stderr.txt
conda deactivate
# -q 0 --paired --clip_R1 20 --three_prime_clip_R1 21 --clip_R2 21 --three_prime_clip_R2 20
echo "Step 2: Trim-galore is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

mv "$directory/$name/${name}_${exon}_R1_trimmed.fq.gz" "$directory/trim/${name}_${exon}_R1_trimmed.fq.gz"

else
echo "Step 2  Trim-galore: $directory/trim/${name}_${exon}_R1_trimmed.fq.gz exists. File is not overwritten. Skipped Step 2." 1>&2
fi

