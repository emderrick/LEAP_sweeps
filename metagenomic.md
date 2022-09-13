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
#!/usr/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=32
module load megahit/1.2.9

megahit -1 L3_1_R1.fastq.gz,L3_2_R1.fastq.gz,L3_3_R1.fastq.gz,L3_4_R1.fastq.gz,L3_5_R1.fastq.gz -2 L3_1_R2.fastq.gz,L3_2_R2.fastq.gz,L3_3_R2.fastq.gz,L3_4_R2.fastq.gz,L3_5_R2.fastq.
gz -t 32 -o L3_coassembly
```
once I have a coassembly for each pond I will follow the anvio pipeline.

**bowtie2**
------
Once we have our coassemblies we need to map the reads at each timepoint to the coassembly.
We will use bowtie2.

```bash
#!/usr/bin/bash
module load bowtie2

for f in *_R1.fastq.gz
do
f2="${f//R1.fastq.gz/R2.fastq.gz}"
out="${f//_R1.fastq.gz/.sam}"
bowtie2 -x I4_contigs -1 $f -2 $f2 --threads 40 -S $out
done
```
**samtools**
------
After running bowtie2 we will have a .sam file for each sample. We need to convert the sam file to a bam file for anvio. Program needed lots of memory to work properly (20G).

```bash
 #!/usr/bin/bash
for f in *I4.sam*
do
out="${f//.sam/raw.bam}"
samtools view -b $f > $out
done
```

Then index the bam file. This program uses samtools and needs lots of memory (10G)
```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a

for f in *raw.bam
do
out="${f//raw.bam/.bam}"
anvi-init-bam $f -o $out
done
deactivate
```

**Anvio**
------
First reformat fasta deflines and remove short contigs (<1000bp). An example with K1

```bash
#!/usr/bin/bash
module load scipy-stack/2021a
module load python/3.7
anvi-script-reformat-fasta K1_contigs.fa -o K1_contigs_fixed.fa -l 1000 --simplify-names
```
Then rename it to get rid of fixed in title
Then make a contigs database of each. This uses prodigal for gene calling. An example with K1

```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.fa
do
out="${f//.fa/.db}"
anvi-gen-contigs-database -f $f -o $out -T 16
done
deactivate

```

Then run HMMS

```bash

#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.db
do
anvi-run-hmms -c $f -T 16
done
deactivate
```

```bash
#!/usr/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=ctb-shapiro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.db
do
anvi-run-scg-taxonomy -c $f
done
deactivate

```

Then we profile each bam file with the corresponding coassembly

```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7

for f in *_I4.bam*
do
anvi-profile -i $f -c I4_contigs.db -T 16
done
deactivate
```

Then we use cluster contigs to make a collection.

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load metabat2/2.14
anvi-cluster-contigs -p I8_merged/PROFILE.db -c I8_contigs.db -C I8_collection --driver metabat2 --just-do-it
```

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
