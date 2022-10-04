#!/bin/bash
#BSUB -P SGE
#BSUB -G ddd-grp
#BSUB -o log/plasmid_trim_all_%J_%I.LSFout.txt -e log/plasmid_trim_all_%J_%I.LSFerr.txt
#BSUB -n 1
#BSUB -R "span[hosts=1] select[mem>300] rusage[mem=300]" -M 300
#BSUB -q normal
#BSUB -J "plasmid_trim_all[1-102]"

id=$LSB_JOBINDEX

#column1
name="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 1)"
#column2
directory="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 2 )"
#column3
read1path="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 3)"
#column4
read2path="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 4)"
#column5
minlength="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 5)"
#column6 comments will not be used

[ ! -d "$directory/tempo" ] && mkdir -p "$directory/tempo" && echo "$directory/tempo is not found. Directory is created." 1>&2

read1name="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 3 | rev | cut -f 1 -d '/' | rev)"
read2name="$(sed 1d trim.txt | head -n $id |tail -n 1| cut -f 4 | rev | cut -f 1 -d '/' | rev)"


cp "$read1path" $directory/tempo
cp "$read2path" $directory/tempo
mv "$directory/tempo/$read1name" "$directory/tempo/${name}_L001_R1_001.fastq.gz"
mv "$directory/tempo/$read2name" "$directory/tempo/${name}_L001_R2_001.fastq.gz"

#Step 1: Trim-galore

if [ ! -e "$directory/$name/${name}_L001_R1_001_val_1.fq.gz" ]
then echo "Step 1: Trim-galore is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

source activate trim_galore
trim_galore --paired --fastqc --length $minlength -o $directory/$name "$directory/tempo/${name}_L001_R1_001.fastq.gz" "$directory/tempo/${name}_L001_R2_001.fastq.gz" 1>$directory/${name}_trimgalore.stdout.txt 2>$directory/${name}_trimgalore.stderr.txt
conda deactivate
# -q 0 --paired --clip_R1 20 --three_prime_clip_R1 21 --clip_R2 21 --three_prime_clip_R2 20
echo "Step 1: Trim-galore is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

rm "$directory/tempo/${name}_L001_R1_001.fastq.gz"
rm "$directory/tempo/${name}_L001_R2_001.fastq.gz"

else
echo "Step 1  Trim-galore: $directory/$name/${name}_L001_R1_001_val_1.fq.gz exists. File is not overwritten. Skipped Step 1." 1>&2
fi

#Step 2: SeqPrep

[ ! -e "$directory/$name/${name}_L001_R1_001_val_1.fq.gz" ] && echo "Step 2 SeqPrep: $directory/$name/${name}_L001_R1_001_val_1.fq.gz is not found. Exits with errors." 1>&2 && exit 1
[ ! -e "$directory/$name/${name}_L001_R2_001_val_2.fq.gz" ] && echo "Step 2 SeqPrep: $directory/$name/${name}_L001_R2_001_val_2.fq.gz is not found. Exits with errors." 1>&2 && exit 1

if [ ! -e "$directory/$name/${name}_merge.fastq.gz" ]
then echo "Step 2: SeqPrep is started at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

SeqPrep -f $directory/$name/${name}_L001_R1_001_val_1.fq.gz -r $directory/$name/${name}_L001_R2_001_val_2.fq.gz -1 $directory/$name/${name}_p1.fastq.gz -2 $directory/$name/${name}_p2.fastq.gz -m 0.001 -q 25 -o 15 -s $directory/$name/${name}_merge_prefilter.fastq.gz 1>$directory/${name}_seqprep.stdout.txt 2>$directory/${name}_seqprep.stderr.txt

#Remove anything that is longer than 310bp
zcat $directory/$name/${name}_merge_prefilter.fastq.gz | paste - - - - | awk  'BEGIN{FS="\t";OFS="\t";} length($2)<311' | tr '\t' '\n' | gzip > $directory/$name/${name}_merge.fastq.gz

echo "Step 2: SeqPrep is completed at $(date +%Y-%m-%d--%H:%M:%S)." 1>&2

else
echo "Step 2 SeqPrep: $directory/$name/${name}_merge.fastq.gz exists. File is not overwritten. Skipped Step 2." 1>&2
fi

