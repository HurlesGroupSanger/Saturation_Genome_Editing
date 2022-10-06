# SGE fastq to count
Required: <br/>
1.) [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://docs.anaconda.com/anaconda/install/index.html)<br/>
2.) [gnu-parallel](https://anaconda.org/conda-forge/parallel)<br/>
3.) Download the .yml files to your home directory and run these: <br/>
`conda env create --name trim_galore --file=trim_galore.yml`<br/>
`conda env create --name cutadapt3 --file=cutadapt.yml`<br/>
`conda env create --name py3.9 --file=pycroquet.yml`<br/>
Notes: In this version of script, all tools are directed to my workplace (with prefix `/nfs/users/nfs_h/hk5/bin/`). You do not need to install all the tools except `conda` and `gnu-parallel` if you run in Sanger Farm5. You can only run this script in Sanger Farm5 with your Farm5 user group.

- [PART1 Demultiplex and trimming](#part1-demultiplex-and-trimming)
  - [Running the code](#running-the-code) 
  - [Metadata file](#metadata-file)
  - [Examples](#examples)
  - [Outputs](#outputs)
  - [Re-run the failed runs](#re-run-the-failed-runs)
- [PART2 Targeton extraction and counting](#part2-targeton-extraction-and-counting)
  - [Running the code](#running-the-code-1)
  - [Metadata file](#metadata-file-1)
  - [Examples](#examples-1)
  - [Outputs](#outputs-1)
  - [Re-run with different parameters](#re-run-with-different-parameters)

## PART1 Demultiplex and trimming
Use [`tre-agrep`](https://anaconda.org/tsnyder/tre) for demultiplexing of a pooled sequencing results based on the constant region sequence. The reads will then be trimmed by [`trim-galore`](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) to remove the low quality base (Default q20 and lower) at the 3'end and then removed the Illumina adaptor at the 3' end. A `FastQC` will be ran after trimming. If it is a pair-end read, the trimmed read will then be merged by [`SeqPrep`](https://github.com/jstjohn/SeqPrep).

You can provide the full length of the 5'constantseq and 3'constantseq (Nucleotide sequence follow the Read1 direction) for the demultiplex. However, only the last 10 bases of the 5'constantseq and the first 10 bases of the 3'constantseq will be used. This is to avoid the bad base quality at the both ends of a read. Also, the shorter the regex, the faster the fuzzy match (`agrep`).

The 5' and 3' constantseq are not anchored. ie: it can have other sequences before the 5'constantseq and after the 3'constantseq.

For example, if the sequence is xxxATTTAGATGA...GGGGAAATTCxxx in the Read1 direction, then 5'constantseq and 3'constantseq which used for `agrep` are ATTTAGATGA and GGGGAAATTC respectively (with regex `ATTTAGATGA.*GGGGAAATTC`). xxx can be any sequence such as Illumina adaptor sequence or primer sequence.

For optimal demultiplexing, use forward primer (the primer with P5 sequence) as 5'constantseq (Column6 in metadata file) and reverse complement of the reverse primer (the primer with P7 sequence) as 3'constantseq (Column7 in metadata file) instead of the full constant sequence. This is to avoid the potential NHEJ indel that may span nearby the targeton region.

**Important:** For pair-end data, the read in read1.fastq.gz must be in the same order as the read in read2.fastq.gz. This is always true if the fastq were downloaded from Basespace or basecalled manually with `bcl2fastq`.

<br/>

#### Single-end mode

_**Raw Read1 Structure**_

[Forward Primer Sequence]-[5'Constant Sequence]-[Targeton]-[3' Constant Sequence]-[Reverse complement of Reverse primer sequence]-[Illumina adaptor]

<br/>

#### Pair-end mode

_**Raw Read1 Structure**_

[Forward Primer Sequence]-[5'Constant Sequence]-[5' end of Targeton and 3' Constant Sequence]

_**Raw Read2 Structure**_

[Reverse Primer Sequence]-[Reverse complement of 3'Constant Sequence]-[Reverse complement of 3' end of Targeton and 5' Constant Sequence]

<br/>
 

### Running the code
```
./bsub_SGE_demultiplex_and_trim_prep.sh -G [farm user group] -m [memory] -c [core] -J [Jobname] -I [Line number] -t [PATH for input tsv] -w [pre-executed job] -s -p [SE/PE] -q [bqueues] -D [agrep errors] -Q [0-40] --sm [SeqPrep -m] --sq [SeqPrep -q] --so [SeqPrep -o] --sl [maximum merged read length]
```

<br/>

**Required:**
| Option | Description |
| ------ | ------ |
| -t,--tsv | Path for the metadata file |
| -J,--Jobname | Unique job name for job submission |   

<br/>

**Optional:**
| Option          | Description |
| ------          | ------ |
| -G              | Farm User Group, Default: ddd-grp |
| -q              | Farm queue. Check available queue with `bqueues`. Default: normal |
| -p              | SE (single-end) or PE (pair-end). Default: SE |
| -m,--memory     | Memory request (in MB). Default: 300 |
| -c,--core       | number of core. Default: 4 |
| -I,--Index      | Line number in the inputtsv for submission (Header is Line 0). eg [1,2,3,100] or [1-20,35,40-50].Default: [1] <br/> For limiting job submitted, used [1-20]%%10. means submit 10 jobs at once. Notes: There are two "%" signs. |
| -s,--submission | `bsub` or not. If flag as -s, the job will be submitted. |
| -w              | Job that has to be pre-excecuted. done(JobID\|Jobname) && done(xxx). It could be either `done` or `started` or `ended` or `exit`. <br/> In the terminal, the input will be like "done\\(xxx\\)". You need to escape the "(" |

<br/>

**agrep and trim related:**
| Option          | Description |
| ------          | ------ |
| -Q              | Quality thresold in `trim-galore` (-q) for trimming. Default: 20. |
| -D              | `agrep` Levenshtein Distance. Use 10 bases of each 5 and 3 constantseq. Default: 1, ie maximum 1 error. |

<br/>

**Merged read related:**
| Option          | Description |
| ------          | ------ |
| --sm              | SeqPrep -m. Will be used only if -p PE. Default: 0.001 |
| --sq              | SeqPrep -q. Will be used only if -p PE. Default: 25 |
| --so              | SeqPrep -o. Will be used only if -p PE. Default: 15 |
| --sl              | Maximum merged read length. Will be used only if -p PE. Default: 310 |

<br/>

**Reserved Options:**
| Option          | Description |
| ------          | ------ |
| --me            | Method for targeton extraction. NULL or cutadapt or tagdust. Default: NULL |
| --mc            | Method for counting. awk or pycroquet. Default: awk |


<br/>


### Metadata file 

An 8-column tab-separated-file
| Column          | Description |
| ------          | ------ |
| Name            | Unique sample name |
| Directory       | Output directory. Must be created before running. No "/" at the end. |
| Read1           | Read1 fastq.gz full path. Must be a gzip file. |
| Read2           | Read2 fastq.gz full path. Must be a gzip file. NA if it is single end data. |
| Exon            | Exon/Targeton label |
| 5'constantseq   | 5' Constant sequence. Same direction as Read1. At least 10 bases. |
| 3'constantseq   | 3' Constant sequence. Same direction as Read1. At least 10 bases. |
| Comment         | Whatever comments. |

`./bsub_SGE_demultiplex_and_trim_prep.sh` Check the following in the metadata file: <br/>
1) Check whether the Read1 (Read1 and Read2 if `-p PE`) path exists. Return error if one of the path does not exist. <br/>
2) Check whether the entries in the Name column (Column 1) are unique. Return error if there is a duplication. 

The current code includes a `dos2unix` step that converts all your input tables to ASCII text. It output a file as dos2unix.log.txt  However, it is important for you to check that your all your tables are ASCII text before using this code start as I may miss out something.  Run (`file xxx.txt`). If it showed as (ASCII text, with CRLF, LF line terminators), please run (`dos2unix xxx.txt`).  You should see (ASCII text) only if you run (`file xxx.txt`)

If the code run successfully, it will generate a new script with name `${jobname}_demultiplex_trim.sh`. This new script can be use for `bsub` if `-s` flag was not included in the `./bsub_SGE_demultiplex_and_trim_prep.sh` run. <br/>

<br/>

### Examples 
```bash
./bsub_SGE_demultiplex_and_trim_prep.sh -J testing -I [1-3] -t demultiplex.txt
```

A `bsub` script with job array was created with job name "testing" and the single end sample in Line1 to Line3 of the "demultplex.txt" were included. Allow maximum 1 error for the `agrep`. `bsub` parameters can be changed manually in the `testing_demultiplex_trim.sh` script generated. The job can be submitted with `cat testing_demultiplex_trim.sh | bsub` 

<br/>

```bash
./bsub_SGE_demultiplex_and_trim_prep.sh -J testing -I [1-3,7-9] -t demultiplex.txt -w done\(seq1\) -p PE -D 2 -Q 0 --sl 400 -s
```

A job array with job name "testing" was submitted to "normal" queue under "ddd-grp" user group. 4 cores and 300MB were requested for each job. The job will be started after the job "seq1" completed. The job array contained the pair end sample in Line1 to Line3 and Line7 to Line9 of the "demultplex.txt". Allow maximum 2 errors for the `agrep`. Do not trim bad quality base before illumina adaptor trimming step (-Q 0). Merged read longer than 400bp will be removed. `testing_demultiplex_trim.sh` script was generated as a reference for this `bsub`.

<br/>

### Outputs
1) `${jobname}_demultiplex_trim.sh`. A `bsub` script.<br/>
2) A "demultiplex" folder will be created. The demultiplex fastq.gz files will be in this folder. <br/>
3) A "trim" folder will be created. The trim_galore output fastq.gz files will be in this folder. <br/>
4) If it is a pair end data, a "merge" folder will be created. The SeqPrep output fastq.gz files will be in this folder. <br/>
5) Folders with the "Name" (Column1 in the metadata file) will be created. The others trim_galore/SeqPrep/FastQC output files will be in these folders. <br/>
6) A log folder will be created if it was not exist before. All LSF output and error log are here.<br/>

<br/>

### Re-run the failed runs

The script runs in three steps:<br/>
1) Demultiplex with agrep. Column 3 (Read1 path) and Column 4 (Read2 path, PE only) as input. Output to "demultiplex".<br/>
2) Trim-galore. Files in "demultiplex" as input. Output to "trim". <br/>
3) PE only, SeqPrep. Files in "trim" as input. Output to "merge". <br/>

The files associate with each sample will have the prefix `${name}_${exon}`. (Column1_Column5 from metadata file.)

**Important:** <br/>
Before the script enter a new step, the script checks whether the required input file is available. If it is unavailable, the run for that sample will be ended. A message will be printed in the `LSFerr` in the log folder.

Also, the script check whether the output file for that step is available. If it is available, it will skip the current step and start the next step. A message will be printed in the `LSFerr` in the log folder.

If you want to re-run certain samples, remember to delete the output file of the sample for each step which you want to re-run.

For example, if I want to re-run SeqPrep of "SGE_sg1_Rep1_Day4_HG100", I need to delete the file with prefix "SGE_sg1_Rep1_Day4_HG100" in the "merge" folder before I re-run the script. It is OK if I accidentally (or purposely) re-run for other samples. The results of other sample will not be replaced/overwritten. The script will not be started and will end almost immediately if it detected all the output files of that samples were available. A message will be printed in the `LSFerr` in the log folder.


<br/>

## PART2 Targeton extraction and counting

Work directly to the output of the demultiplexing and trim/merge. The trimmed/merged fastq.gz will be either 1) Remove the prime/constant sequence by [`tagdust2`](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0454-y), 2) Remove the prime/constant sequence by [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/), 3) Do not remove prime/constant sequence. Then, the extracted fastq.gz will be counted by either 1) `awk` or 2) [`pycroquet`](https://github.com/cancerit/pycroquet)

**Read Extraction:**<br/>
1) NULL: Do not trim any bases (default)<br/>
2) tagdust: Use `tagdust` for trimming. Need to set `--memory 2000` when used. Use more core to speed up the process eg. `--core 10`. Slow but good for data with bad sequencing quality. The error rate of the adaptor are estimated by hidden markov model built based on the dataset provided. Minimum length of the extracted read is 1 base. The extraction might fail if you provided the adaptor sequence which are too far away from both ends of the read. The HMM model assume the 5constantseq are near to the 5'end and 3constantseq are near to the 3'end of the read. <br/>
3) cutadapt: Use `cutadapt` for trimming. Use more core to speed up the process eg. `--core 10`. Error rate used for adaptor searching is the default of the `cutadapt`. Minimum length of the extracted read is 1 base. The provided adaptor sequence are not neccessary nearby the ends of the read. see `-g ADAPTER1...ADAPTER2` in `cutadapt` manual.<br/>

**Read Counting:**<br/>
1) awk: `zcat ${input.fq.gz} | paste - - - - | cut -f 2 | awk '{ cnts[$0] += 1 } END { for (v in cnts) print v"\t"cnts[v]}' > ${output}` (default)
2) pycroquet: `pycroquet long-read --unique -q ${input.fq.gz} -o ${output} -s ${name}`

<br/>

### Running the code

```
./bsub_extract_and_count.sh -G [farm user group] -m [memory] -c [core] -J [Jobname] -I [Line number] -t [PATH for metadata file] -w [pre-executed job] -s -p [SE/PE] -q [bqueues] --me [NULL/cutadapt/tagdust] --mc [awk/pycroquet]
```

<br/>

**Required:**
| Option | Description |
| ------ | ------ |
| -t,--tsv | Path for the metadata file |
| -J,--Jobname | Unique job name for job submission |   

<br/>

**Optional:**
| Option          | Description |
| ------          | ------ |
| -G              | Farm User Group, Default: ddd-grp |
| -q              | Farm queue. Check available queue with `bqueues`. Default: normal |
| -p              | SE (single-end) or PE (pair-end). Default: SE |
| -m,--memory     | Memory request (in MB). Default: 300 |
| -c,--core       | number of core. Default: 4 |
| -I,--Index      | Line number in the inputtsv for submission (Header is Line 0). eg [1,2,3,100] or [1-20,35,40-50].Default: [1] <br/> For limiting job submitted, used [1-20]%%10. means submit 10 jobs at once. Notes: There are two "%" signs. |
| -s,--submission | `bsub` or not. If flag as -s, the job will be submitted. |
| -w              | Job that has to be pre-excecuted. done(JobID\|Jobname) && done(xxx). It could be either `done` or `started` or `ended` or `exit`. <br/> In the terminal, the input will be like "done\\(xxx\\)". You need to escape the "(" |

<br/>

**Read extraction and count related:**
| Option          | Description |
| ------          | ------ |
| --me            | Method for targeton extraction. NULL or cutadapt or tagdust. Default: NULL |
| --mc            | Method for counting. awk or pycroquet. Default: awk |

<br/>

**Reserved Options:**
| Option            | Description |
| ------            | ------ |
| -Q                | Quality thresold in `trim-galore` (-q) for trimming. Default: 20. |
| -D                | `agrep` Levenshtein Distance. Use 10 bases of each 5 and 3 constantseq. Default: 1, ie maximum 1 error. |
| --sm              | SeqPrep -m. Will be used only if -p PE. Default: 0.001 |
| --sq              | SeqPrep -q. Will be used only if -p PE. Default: 25 |
| --so              | SeqPrep -o. Will be used only if -p PE. Default: 15 |
| --sl              | Maximum merged read length. Will be used only if -p PE. Default: 310 |

<br/>

### Metadata file
This is the same as the demultiplex section. The 5' and 3' constantseq here will be the full length of primer + constant seq if you wish to trim away all the constant seq before counting.

I usually trim away all the constant sequence. So, I replace the primer sequence in 5' and 3' constantseq column to the full length of primer + constant seq. 

<br/>

### Examples 
```bash
./bsub_extract_and_count.sh -J counting -I [1-3] -t demultiplex.txt -w done\(demultiplex\) -p PE --core 10 -q long --memory 2200 --me tagdust --mc pycroquet -s
```

A job array with job name "counting" was submitted to "long" queue under "ddd-grp" user group. 10 cores and 2200MB RAM were requested for each job. The job will be started after the job "demultiplex" completed. The job array contained the pair end sample in Line1 to Line3 of the "demultplex.txt". Since this is pair end data, the merged read will be used for read extraction. The constant sequence will be removed by `tagdust2` and the unique sequence will be counted by `pycroquet`. `counting_extract_count.sh` script was generated as a reference for this `bsub`.

<br/>

### Outputs
1) `${jobname}_extract_count.sh`. A `bsub` script.<br/>
2) A "tempo" folder will be created. It stores some temporary files. It can be deleted manually after the run. <br/>
3) A "extracted" folder will be created. The extracted read (fq.gz) will be in this folder. Read without adaptor (un.fq.gz) and the extraction-related log files will be exported here as well. <br/>
4) A "count" folder will be created. The count table (.txt, tab-separated-file) will be in this folder. The `pycroquet` log files will be exported here as well. <br/>
5) A log folder will be created if it was not exist before. All LSF output and error log are here.<br/>

<br/>

### Re-run with different parameters

The script runs in three steps:<br/>
1) Files in "trim" will be moved to "tempo". If -p PE, files in "merge" will be moved to "tempo" and renamed.<br/>
2) Read Extraction. Files in "tempo" as input. Output to "extracted". <br/>
3) Count. Files in "extracted" as input. Output to "count". <br/>

The files associate with each sample will have the prefix `${name}_${exon}`. (Column1_Column5 from metadata file.)

**Important:** <br/>
Before the script enter a new step, the script checks whether the required input file is available. If it is unavailable, the run for that sample will be ended. A message will be printed in the `LSFerr` in the log folder.

Also, the script check whether the output file for that step is available. If it is available, it will skip the current step and start the next step. A message will be printed in the `LSFerr` in the log folder.

If you want to re-run certain samples, remember to delete the output file of the sample for each step which you want to re-run. OR, rename the "extracted" and "count" folder before running.

For example, if I had run the script with 'tagdust' and 'awk' and I want to re-run the script with 'cutadapt' and 'awk', I will rename the existing "extracted" folder to "tagdust_awk_extracted" and the "count" folder to "tagdust_awk_count" before rerun the script. Since you have changed the path name, the script no longer detected the outputs files, it will rerun the script with the new parameter. A new "extracted" and "count" folder will be created.
