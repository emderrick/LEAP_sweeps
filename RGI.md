**RGI**
------

RGI on metagenomic reads. Program needs memory and cpus. Followed installation guide on rgi github for CARD database and wildcard data.

```bash
#!/usr/bin/bash
source /home/ederrick/scratch/leapmeta/goodfastq/rgi/bin/activate
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3
module load busco/5.2.2
module load scipy-stack/2021a
module load python
module load bowtie2
module load samtools/1.11
module load bamtools
module load bedtools
module load diamond
module load blast+

rgi load --card_json /home/ederrick/scratch/leapmeta/goodfastq/card.json --local
rgi load --wildcard_annotation /home/ederrick/scratch/leapmeta/goodfastq/wildcard_database_vversion_number.fasta  --wildcard_index /home/ederrick/scratch/leapmeta/goodfastq/wildcard
/index-for-model-sequences.txt --card_annotation card_database_v3.2.2.fasta --local

for file in *R1.fastq.gz*
do
filetwo="${file//R1.fastq.gz/R2.fastq.gz}"
out="${file//R1.fastq.gz/rgi}"
rgi bwt -1 $file -2 $filetwo -a bowtie2 -o $out  --local -n 16
done
deactivate
```

RGI on contigs. Program needs memory and cpus.

```bash
#!/bin/bash
source rgi/bin/activate
module load StdEnv/2020 gcc/9.3.0
module load scipy-stack/2021a
module load prodigal
module load bowtie2
module load samtools/1.10
module load diamond
module load blast+

for file in *.fa
do
out="${file//.fa/_rgioutput}"
rgi load --card_json /home/ederrick/scratch/leapmeta/individual_contigs/card.json --local
rgi main -i $file -t contig -n 40 -o $out --local --low_quality --exclude_nudge --clean
done
deactivate
```
