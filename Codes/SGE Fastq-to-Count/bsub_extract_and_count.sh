#!/bin/bash
#Version 25 Jan 2022, Hong Kee Tan
#Edited on 1 March 2022, edited the conda environment path. I provide .yml file now.
set -u

#PartI: Configuration

pname=$(basename $0)

G_flag='ddd-grp'
q_flag='normal'
m_flag='300'
c_flag='4'
J_flag='NULL'
I_flag='[1]'
t_flag='NULL'
s_flag='FALSE'
w_flag=''
p_flag='SE'
D_flag='1'
Q_flag='20'
sm_flag='0.001'
sq_flag='25'
so_flag='15'
sl_flag='310'
me_flag='NULL'
mc_flag='awk'

print_usage() {
  printf "LSF job array submission.
  Usage: ./$pname -G [farm user group] -m [memory] -c [core] -J [Jobname] -I [Line number] -t [PATH for input tsv] -w [pre-executed job] -s -p [SE/PE] -q [bqueues] -D [agrep errors] -Q [0-40] --sm [SeqPrep -m] --sq [SeqPrep -q] --so [SeqPrep -o] --sl [maximum merged read length] --me [NULL/cutadapt/tagdust] --mc [awk/pycroquet]
  
  Required:
  -t,--tsv                Path for the metadata file
  -J,--Jobname            Unique job name for job submission
  
  Optional:
  -G                      Farm User Group, Default: ddd-grp
  -q                      Farm queue. Check available queue with bqueues. Default: normal
  -p                      SE (single-end) or PE (pair-end). Default: SE
  -m,--memory             Memory request (in MB). Default: 300
  -c,--core               number of core. Default: 4
  -I,--Index              Line number in the inputtsv for submission (Header is Line 0). eg [1,2,3,100] or [1-20,35,40-50].
                          Default: [1]
                          For limiting job submitted, used [1-20]%%%%10. means submit 10 jobs at once. Notes: There are two \"%%\" signs.
  -s,--submission         bsub or not. If flag as -s, the job will be submitted.
  -w                      Job that has to be pre-excecuted. done(JobID|Jobname) && done(xxx). It could be either done || started || ended || exit.
                          In the terminal, the input will be like \"done\(xxx\)\". You need to escape the \"(\".

  agrep and trim:
  -Q                      Quality thresold in trim-galore for trimming. Default: 20
  -D                      agrep Levenshtein Distance. Use 10 bases of each 5 and 3 constantseq. Default: 1, ie maximum 1 error.
  
  Merged Read:
  --sm                    SeqPrep -m. Default: 0.001
  --sq                    SeqPrep -q. Default: 25
  --so                    SeqPrep -o. Default: 15
  --sl                    Maximum merged read length. Default: 310
  
  Extraction and Count:
  --me                    Method for targeton extraction. NULL or cutadapt or tagdust. Default: NULL
  --mc                    Method for counting. awk or pycroquet. Default: awk
                          
  eg. ./$pname -G ddd-grp -q normal -m 2046 --core 20 -J testing -I [1-3] -t demultiplex.txt -p SE \n
  eg. ./$pname -G ddd-grp -q normal -m 2046 --core 20 -J testing -I [1-3] -t demultiplex.txt -w done\(seq1\) -p SE -s
  
  For the constant sequence, the direction follows the R1 read. If the constant seq is xxxATTTAGATGA...GGGGAAATTCxxx in the R1 direction, then the constant5 and the constant3 which used for agrep are ATTTAGATGA and GGGGAAATTC respectively. (with regex ATTTAGATGA.*GGGGAAATTC)
  
  The current code includes a dos2unix step that converts all your input tables to ASCII text. It output a file as dos2unix.log.txt
  However, it is important for you to check that your all your tables are ASCII text before using this code start as I may miss out something.
  Run (file xxx.txt). If it showed as (ASCII text, with CRLF, LF line terminators), please run (dos2unix xxx.txt).
  You should see (ASCII text) only if you run (file xxx.txt)\n
  "
}

set -e

# read the options
TEMP=`getopt -o sm:c:t:J:I:w:G:p:q:D:Q: --long submission,memory:,core:,tsv:,Jobname:,Index:,sm:,sq:,so:,sl:,me:,mc: -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -s|--submission)
            s_flag="TRUE" ; shift ;;
        -m|--memory)
            m_flag=$2 ; shift 2;;
        -c|--core)
            c_flag=$2 ; shift 2;;
        -t|--tsv)
            t_flag=$2 ; shift 2;;
        -J|--Jobname)
            J_flag=$2 ; shift 2;;
        -I|--Index)
            I_flag=$2 ; shift 2;;
        -w)
            w_flag=$2 ; shift 2;;
        -G)
            G_flag=$2 ; shift 2;;
        -p)
            p_flag=$2 ; shift 2;;
        -q)
            q_flag=$2 ; shift 2;;
        -D)
            D_flag=$2 ; shift 2;;
        -Q)
            Q_flag=$2 ; shift 2;;
        --sm)
            sm_flag=$2 ; shift 2;;
        --sq)
            sq_flag=$2 ; shift 2;;
        --so)
            so_flag=$2 ; shift 2;;
        --sl)
            sl_flag=$2 ; shift 2;;
        --me)
            me_flag=$2 ; shift 2;;
        --mc)
            mc_flag=$2 ; shift 2;;
        --) shift ; break ;;
        *) print_usage ; exit 1 ;;
    esac
done

unset TEMP

set +e

#Rename the options and check the required flags
inputtsv=$t_flag

if [ "$inputtsv" = "NULL" ]
then
echo "./$pname -t PATH for input file. -t (--tsv) is required. Script aborted"
print_usage
exit 1
fi

[ ! -f "$inputtsv" ] && echo "input tsv does not exist. Script aborted." && exit 1 || echo "input table exist"

jobname=$J_flag

if [ "$jobname" = "NULL" ]
then
echo "./$pname -J Jobname. -J (--Jobname) is required. Script aborted"
print_usage
exit 1
fi

line=$I_flag
core=$c_flag
memory=$m_flag
submission=$s_flag
Group=$G_flag
queue=$q_flag

###Check the metadata file

#Check metadata file entry exist

dos2unix $inputtsv 2>> dos2unix.log.txt

cat $inputtsv | sed 1d > checkfile0

parallel -a 'checkfile0' -j 1 --colsep '\t' " [ ! -f {3} ] && echo {3} does not exist. Script stops."

if [ "$p_flag" = "PE" ]
then
parallel -a 'checkfile0' -j 1 --colsep '\t' " [ ! -f {4} ] && echo {4} does not exist. Script stops."
if [[ $(parallel -a 'checkfile0' -j 1 --colsep '\t' " [ ! -f {3} ] && echo {3}") ]] || [[ $(parallel -a 'checkfile0' -j 1 --colsep '\t' " [ ! -f {4} ] && echo {4}") ]]
then
exit 1
fi

else

if [[ $(parallel -a 'checkfile0' -j 1 --colsep '\t' " [ ! -f {3} ] && echo {3}") ]]
then
exit 1
fi
fi
rm checkfile0


#Check Name uniqness in metadata file
if [[ $(cat $inputtsv | sed 1d | awk '{print $1}' | sort | uniq -d) ]]
then
echo "At least two samples have same Name (column 1) in inputtsv. Script Stop."
exit 1
fi


[ ! -d "log" ] && mkdir -p "log"

#Part2: Print code

####Header

print_code_SE() {
printf "#BSUB -P SGE
#BSUB -G $Group
#BSUB -o log/${jobname}_%%J_%%I.LSFout.txt -e log/${jobname}_%%J_%%I.LSFerr.txt
#BSUB -n $core
#BSUB -R \"span[hosts=1] select[mem>$memory] rusage[mem=$memory]\" -M $memory
#BSUB -q $queue
#BSUB -J \"$jobname$line\"

id=\$LSB_JOBINDEX

#column1
name=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 1)\"
#column2
directory=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 2 )\"
#column3
read1path=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 3)\"
#column4
#read2path=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 4)\"
#column5
exon=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 5)\"
#column6
constant5=\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 6 )
#column7
constant3=\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 7)
#column8 comments will not be used

[ ! -d \"\$directory/tempo\" ] && mkdir -p \"\$directory/tempo\" && echo \"\$directory/tempo is not found. Directory is created.\" 1>&2

cp \"\$directory/trim/\${name}_\${exon}_R1_trimmed.fq.gz\" \$directory/tempo

"
}

print_code_PE() {
printf "#BSUB -P SGE
#BSUB -G $Group
#BSUB -o log/${jobname}_%%J_%%I.LSFout.txt -e log/${jobname}_%%J_%%I.LSFerr.txt
#BSUB -n $core
#BSUB -R \"span[hosts=1] select[mem>$memory] rusage[mem=$memory]\" -M $memory
#BSUB -q $queue
#BSUB -J \"$jobname$line\"

id=\$LSB_JOBINDEX

#column1
name=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 1)\"
#column2
directory=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 2 )\"
#column3
read1path=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 3)\"
#column4
#read2path=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 4)\"
#column5
exon=\"\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 5)\"
#column6
constant5=\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 6 )
#column7
constant3=\$(sed 1d $inputtsv | head -n \$id |tail -n 1| cut -f 7)
#column8 comments will not be used

[ ! -d \"\$directory/tempo\" ] && mkdir -p \"\$directory/tempo\" && echo \"\$directory/tempo is not found. Directory is created.\" 1>&2

cp \"\$directory/merge/\${name}_\${exon}_merge.fastq.gz\" \$directory/tempo
mv \"\$directory/tempo/\${name}_\${exon}_merge.fastq.gz\" \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\"

"
}


####Method for Extraction

if [ "$me_flag" = "NULL" ]
then
print_code_extract() {
printf "#Step 1: Extract

[ ! -d \"\$directory/extracted\" ] && mkdir -p \"\$directory/extracted\" && echo \"\$directory/extracted is not found. Directory is created.\" 1>&2

[ ! -e \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\" ] && echo \"Step 1 Extract: \$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz is not found. Exits with errors.\" 1>&2 && exit 1

if [ ! -e \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\" ]
then echo \"Step 1: Extract is started at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

cp \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\" \$directory/extracted
mv \"\$directory/extracted/\${name}_\${exon}_R1_trimmed.fq.gz\" \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\"
rm \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\"

echo \"Step 1: Extract is completed at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

else
echo \"Step 1: Extract: \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz. File is not overwritten. Skipped Step 1.\" 1>&2
fi

"
}

elif [ "$me_flag" = "tagdust" ]
then
print_code_extract() {
printf "#Step 1: Extract with tagdust

[ ! -d \"\$directory/extracted\" ] && mkdir -p \"\$directory/extracted\" && echo \"\$directory/extracted is not found. Directory is created.\" 1>&2

[ ! -e \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\" ] && echo \"Step 1 Extract: \$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz is not found. Exits with errors.\" 1>&2 && exit 1

if [ ! -e \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\" ]
then echo \"Step 1: Extract is started at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

/nfs/users/nfs_h/hk5/bin/tagdust -t $core -1 S:\$constant5 -2 R:N -3 S:\$constant3 -minlen 1 -o \$directory/extracted/\${name}_\${exon}_R1_extracted \$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz 1> \$directory/extracted/\${name}_\${exon}_R1_tagdust_extracted.stdout.txt 2> \$directory/extracted/\${name}_\${exon}_R1_tagdust_extracted.stderr.txt

parallel -j 2 \"gzip {}\" ::: \$directory/extracted/\${name}_\${exon}_R1_extracted.fq \$directory/extracted/\${name}_\${exon}_R1_extracted_un.fq

rm \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\"

echo \"Step 1: Extract is completed at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

else
echo \"Step 1: Extract: \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz. File is not overwritten. Skipped Step 1.\" 1>&2
fi

"
}

elif [ "$me_flag" = "cutadapt" ]
then
print_code_extract() {
printf "#Step 1: Extract with cutadapt

[ ! -d \"\$directory/extracted\" ] && mkdir -p \"\$directory/extracted\" && echo \"\$directory/extracted is not found. Directory is created.\" 1>&2

[ ! -e \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\" ] && echo \"Step 1 Extract: \$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz is not found. Exits with errors.\" 1>&2 && exit 1

if [ ! -e \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\" ]
then echo \"Step 1: Extract is started at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

source activate cutadapt3

adaptor=\$(echo \$constant5...\$constant3)

cutadapt -m 1 -g \$adaptor -o \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz --untrimmed-output \$directory/extracted/\${name}_\${exon}_R1_extracted_un.fq.gz \$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz  1> \$directory/extracted/\${name}_\${exon}_R1_cutadapt_extracted.stdout.txt 2> \$directory/extracted/\${name}_\${exon}_R1_cutadapt_extracted.stderr.txt

conda deactivate

rm \"\$directory/tempo/\${name}_\${exon}_R1_trimmed.fq.gz\"

echo \"Step 1: Extract is completed at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

else
echo \"Step 1: Extract: \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz. File is not overwritten. Skipped Step 1.\" 1>&2
fi

"
}

else
echo "--me NULL or cutadapt or tagdust only. Script stop." && exit 1
fi


####Method for Counting

if [ "$mc_flag" = "awk" ]
then
print_code_count() {
printf "#Step 2: Counting

[ ! -d \"\$directory/count\" ] && mkdir -p \"\$directory/count\" && echo \"\$directory/count is not found. Directory is created.\" 1>&2

[ ! -e \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\" ] && echo \"Step 2 Counting: \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz is not found. Exits with errors.\" 1>&2 && exit 1

if [ ! -e \"\$directory/count/\${name}_\${exon}_all_count.txt\" ]
then echo \"Step 2: Counting is started at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

zcat \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\"| paste - - - - | cut -f 2 | awk '{ cnts[\$0] += 1 } END { for (v in cnts) print v\"\\\t\"cnts[v] }' | sort -n -r -k 2,2 > \"\$directory/count/\${name}_\${exon}_all_count.txt\"

echo \"Step 2: Counting is completed at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

else
echo \"Step 2 Counting: \$directory/count/\${name}_\${exon}_all_count.txt exists. File is not overwritten. Skipped Step 2.\" 1>&2
fi

"
}

elif [ "$mc_flag" = "pycroquet" ]
then
print_code_count() {
printf "#Step 2: Counting with pycroquet

[ ! -d \"\$directory/count\" ] && mkdir -p \"\$directory/count\" && echo \"\$directory/count is not found. Directory is created.\" 1>&2

[ ! -e \"\$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz\" ] && echo \"Step 2 Counting: \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz is not found. Exits with errors.\" 1>&2 && exit 1

if [ ! -e \"\$directory/count/\${name}_\${exon}_all_count.txt\" ]
then echo \"Step 2: Counting is started at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

source activate py3.9

pycroquet long-read --unique -q \$directory/extracted/\${name}_\${exon}_R1_extracted.fq.gz -o \$directory/count/\${name}_\${exon}_all_count -s \${name}_\${exon} 1> \$directory/count/\${name}_\${exon}_pycroquet_count.stdout.txt 2> \$directory/count/\${name}_\${exon}_pycroquet_count.stderr.txt

conda deactivate

zcat \$directory/count/\${name}_\${exon}_all_count.query_counts.tsv.gz > \$directory/count/\${name}_\${exon}_all_count.txt

rm \$directory/count/\${name}_\${exon}_all_count.query_counts.tsv.gz

echo \"Step 2: Counting is completed at \$(date +%%Y-%%m-%%d--%%H:%%M:%%S).\" 1>&2

else
echo \"Step 2 Counting: \$directory/count/\${name}_\${exon}_all_count.txt exists. File is not overwritten. Skipped Step 2.\" 1>&2
fi

"
}

else
echo "--mc awk or pycroquet only. Script stop." && exit 1
fi



#Part3: Assemble the final script

#Post excecution or not
if [ "$w_flag" = "" ]
then

    if [ "$p_flag" = "PE" ]
    then
    (echo "#!/bin/bash"; print_code_PE; print_code_extract; print_code_count) > ${jobname}_extract_count.sh

    elif [ "$p_flag" = "SE" ]
    then
    (echo "#!/bin/bash"; print_code_SE; print_code_extract; print_code_count) > ${jobname}_extract_count.sh

    else
    echo "SE or PE only."
    fi
    
else

    if [ "$p_flag" = "PE" ]
    then
    (echo "#!/bin/bash"; echo "#BSUB -w $w_flag"; print_code_PE; print_code_extract; print_code_count) > ${jobname}_extract_count.sh

    elif [ "$p_flag" = "SE" ]
    then
    (echo "#!/bin/bash"; echo "#BSUB -w $w_flag"; print_code_SE; print_code_extract; print_code_count) > ${jobname}_extract_count.sh

    else
    echo "SE or PE only."
    fi
fi

#Submit or not
if [ "$submission" = "TRUE" ]
then
cat ${jobname}_extract_count.sh | bsub

elif [ "$submission" = "FALSE" ]
then
echo "Job was not submitted. Use -s (--submission) for job submission."
fi
