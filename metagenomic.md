# Workflow for the metagenomic data

**Trimmomatic**
------
First trim the reads to remove the Illumina adaptors. I sumbitted this on Narval in a bash script.
I downloaded the adaptor file from the trimmomatic github page it is TruSeq3-PE.fa

```bash
#!/usr/bin/bash
for file in *R1_001.fastq.gz*
do
filetwo="${file//R1_001.fastq.gz/R2_001.fastq.gz}"
paired1="${file//R1_001.fastq.gz/paired_R1_001.fastq.gz}"
unpaired1="${file//R1_001.fastq.gz/unpaired_R1_001.fastq.gz}"
paired2="${filetwo//R2_001.fastq.gz/paired_R2_001.fastq.gz}"
unpaired2="${filetwo//R2_001.fastq.gz/unpaired_R2_001.fastq.gz}"
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 $file $
filetwo $paired1 $unpaired1 $paired2 $unpaired2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:
10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```
**MEGAHIT**
------
```bash

```
**RGI**
------
```bash

```
do this for contigs and metagenomic reads
