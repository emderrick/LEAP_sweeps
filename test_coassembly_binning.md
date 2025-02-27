### Complete workflow for raw reads to instrain output

#### trim adaptors and low quality seqeunces with Trimmomatic

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bashrc
conda activate trimmomatic

for f in *R1.fastq.gz
do
trimmomatic PE -threads 32 -phred33 $f ${f%*R1.fastq.gz}R2.fastq.gz \
${f%*R1.fastq.gz}P_R1.fastq.gz ${f%*R1.fastq.gz}UP_R1.fastq.gz ${f%*R1.fastq.gz}P_R2.fastq.gz ${f%*R1.fastq.gz}UP_R2.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

#### Create a co-assembly of all ponds at time point 1. Try two ways and compare quality.

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate megahit
megahit -1 I4_1_R1.fastq.gz,I8_1_R1.fastq.gz,K1_1_R1.fastq.gz,L3_1_R1.fastq.gz,L4_1_R1.fastq.gz,L2_1_R1.fastq.gz,L6_1_R1.fastq.gz,L7_1_R1.fastq.gz,L8_1_R1.fastq.gz \
-2 I4_1_R2.fastq.gz,I8_1_R2.fastq.gz,K1_1_R2.fastq.gz,L3_1_R2.fastq.gz,L4_1_R2.fastq.gz,L2_1_R2.fastq.gz,L6_1_R2.fastq.gz,L7_1_R2.fastq.gz,L8_1_R2.fastq.gz \
-o T1_coassembly --min-contig-len 1000 --threads 16
conda deactivate
```

```bash
cat *1_R1.fastq > T1_R1.fastq
cat *1_R2.fastq > T1_R2.fastq
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate spades
/mfs/ederrick/miniconda3/envs/spades/bin/metaspades.py -1 all_1_R1.fastq -2 all_1_R2.fastq -o T1_metaspades_coassembly --only-assembler -m 1000 --threads 16
conda deactivate
```

#### compare assembly qualities. They end up very similar but MEGAHIT is much faster.

```bash
seqkit seq -m 1000 metaspades_contigs.fa > metaspades_contigs_filt.fa
seqkit stats -a T1_coassembly.fa > T1_coassembly_spades_stats.txt
seqkit stats -a metaspades_contigs_filt.fa > metaspades_contigs_stats.txt
```

#### Map metagenomic reads from each pond at T1 to the co-assembly

```bash
bowtie2-build T1_coassembly.fa T1_coassembly --threads 64
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_coassembly -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/_R1.fastq.gz/.bam} --write-index -@ 8' ::: *_1_R1.fastq.gz

mkdir megahit_bam
mv *.bam megahit_bam
```

### Map metagenomic reads from each pond at T1 to the spades co-assembly

```bash
bowtie2-build metaspades_contigs_filt.fa spades_coassembly --threads 64
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x spades_coassembly -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/_R1.fastq.gz/_spades.bam} --write-index -@ 16' ::: *_1_R1.fastq.gz

mkdir spades_bam
mv *spades.bam spades_bam
```

#### bin contigs with metabat2 (bam files in their respective directories)

```bash
jgi_summarize_bam_contig_depths --outputDepth megahit_depth.txt *.bam
metabat2 -i T1_coassembly.fa -a megahit_depth.txt -o megahit_bins/bin -m 2500 -t 64

jgi_summarize_bam_contig_depths --outputDepth spades_depth.txt *.bam
metabat2 -i metaspades_contigs_filt.fa -a spades_depth.txt -o spades_bins/bins -m 2500 -t 64
```

#### check quality and dereplicate MAGs with dREP and checkM

```bash
checkm lineage_wf megahit_bins megahit_bin_quality -t 64 -x fa
checkm lineage_wf spades_bins spades_bin_quality -t 64 -x fa

dRep compare output_directory -g path/to/genomes/*.fasta
dRep dereplicate output_directory -g path/to/genomes/*.fasta
```

#### Map to the MAG database

```bash
cat *megahit_MAG* > megahit_MAGs.fa
bowtie2-build megahit_MAGs.fa megahit_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x megahit_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/megahit_MAGs.bam} --write-index -@ 8' ::: *1_R1.fastq.gz
```

```bash
cat *spades_MAG* > spades_MAGs.fa
bowtie2-build spades_MAGs.fa spades_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x spades_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/spades_MAGs.bam} --write-index -@ 8' ::: *1_R1.fastq.gz
```
