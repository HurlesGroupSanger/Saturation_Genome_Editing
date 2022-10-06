#!/bin/bash
#BSUB -P SGE
#BSUB -G ddd-grp
#BSUB -o log/counting_%J_%I.LSFout.txt -e log/counting_%J_%I.LSFerr.txt
#BSUB -n 5
#BSUB -R "span[hosts=1] select[mem>2200] rusage[mem=2200]" -M 2200
#BSUB -q normal
#BSUB -J "counting[1-180]"

id=$LSB_JOBINDEX

#column1
name="$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 1)"
#column2
directory="$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 2 )"
#column3
read1path="$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 3)"
#column4
#read2path="$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 4)"
#column5
exon="$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 5)"
#column6
constant5=$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 6 )
#column7
constant3=$(sed 1d trim_primer.txt | head -n $id |tail -n 1| cut -f 7)
#column8 comments will not be used

[ ! -d "$directory/tempo" ] && mkdir -p "$directory/tempo" && echo "$directory/tempo is not found. Directory is created." 1>&2

cp "$directory/trim/${name}_${exon}_R1_trimmed.fq.gz" $directory/tempo

#Step 1: Extract with cutadapt

[ ! -d "$directory/extracted" ] && mkdir -p "$directory/extracted" && echo "$directory/extracted is not found. Directory is created." 1>&2

[ ! -e "$directory/tempo/${name}_${exon}_R1_trimmed.fq.gz" ] && echo "Step 1 Extract: $directory/tempo/${name}_${exon}_R1_trimmed.fq.gz is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/extracted/${name}_${exon}_R1_extracted.fq.gz" ]
then echo "Step 1: Extract is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

source activate /nfs/users/nfs_h/hk5/miniconda2/envs/cutadapt3

adaptor=$(echo $constant5...$constant3)

cutadapt -j 5 -m 1 -g $adaptor -o $directory/extracted/${name}_${exon}_R1_extracted.fq.gz --untrimmed-output $directory/extracted/${name}_${exon}_R1_extracted_un.fq.gz $directory/tempo/${name}_${exon}_R1_trimmed.fq.gz  1> $directory/extracted/${name}_${exon}_R1_cutadapt_extracted.stdout.txt 2> $directory/extracted/${name}_${exon}_R1_cutadapt_extracted.stderr.txt

conda deactivate

rm "$directory/tempo/${name}_${exon}_R1_trimmed.fq.gz"

echo "Step 1: Extract is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 1: Extract: $directory/extracted/${name}_${exon}_R1_extracted.fq.gz. File is not overwritten. Skipped Step 1." 1>&2
fi

#Step 2: Counting

[ ! -d "$directory/count" ] && mkdir -p "$directory/count" && echo "$directory/count is not found. Directory is created." 1>&2

[ ! -e "$directory/extracted/${name}_${exon}_R1_extracted.fq.gz" ] && echo "Step 2 Counting: $directory/extracted/${name}_${exon}_R1_extracted.fq.gz is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/count/${name}_${exon}_all_count.txt" ]
then echo "Step 2: Counting is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

zcat "$directory/extracted/${name}_${exon}_R1_extracted.fq.gz"| paste - - - - | cut -f 2 | awk '{ cnts[$0] += 1 } END { for (v in cnts) print v"\t"cnts[v] }' | sort -n -r -k 2,2 > "$directory/count/${name}_${exon}_all_count.txt"

echo "Step 2: Counting is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 2 Counting: $directory/count/${name}_${exon}_all_count.txt exists. File is not overwritten. Skipped Step 2." 1>&2
fi

