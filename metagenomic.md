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

```
**samtools**
------
After running bowtie2 we will have a .sam file for each sample. We need to convert the sam file to a bam file for anvio.

```bash
 samtools view -b I4_1_I4.sam > I4_1_I4_raw.bam
```
anvio needs the bam files to be indexed so we need to run this for each sample (do in loop)

```bash
anvi-init-bam K1_1_raw.bam -o K1_1.bam
```

**Anvio**
------
First reformat fasta deflines and remove short contigs (<1000bp). An example with K1

```bash
module load scipy-stack/2021a
module load python/3.7
anvi-script-reformat-fasta K1_contigs.fa -o K1_contigs_fixed.fa -l 1000 --simplify-names
```

Then make a contigs database of each. This uses prodigal for gene calling. An example with K1

```bash
anvi-gen-contigs-database -f K1_contigs-fixed.fa -o K1_contigs.db -n 'K1 contigs database'
```

Then run HMMS

```bash
anvi-run-hmms -c contigs.db
```

To take a look at number of bacteria present and other stats etc run 

```bash
anvi-display-contigs-stats K1_contigs.db
```
Then we profile each bam file with the corresponding coassembly

```bash
anvi-profile -i K1_1.bam -c contigs_K1.db --min-contig-length 2500 --output-dir K1_1_profile --sample-name K1_1
```
**RGI**
------
```bash

```
do this for contigs and metagenomic reads
