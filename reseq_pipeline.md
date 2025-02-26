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

#### Map metagenomic reads from each pond at T1 to the co-assembly

```bash
bowtie2-build T1_coassembly.fa T1_coassembly
```

```bash
#!/usr/bin/bash
for f in *_R1.fastq.gz; do bowtie2 -x T1_coassembly -1 $f -2 ${f%*R1.fastq.gz}R2.fastq.gz -S ${f%*_R1.fastq.gz}.sam --threads 32
done
```

```bash
for f in *.sam; do samtools view -b $f > ${f%*sam}bam --threads 8; done
for f in *.bam; do samtools sort $f > ${f%*.bam}_sorted.bam --threads 8; done
```

#### bin contigs with metabat2

```bash
conda create -n metabat2
conda activate metabat2
conda install metabat2
runMetaBat.sh -i T1_coassembly.fasta I4_1_sorted.bam I8_1_sorted.bam K1_1_sorted.bam L3_1_sorted.bam L4_1_sorted.bam L2_1_sorted.bam L6_1_sorted.bam L7_1_sorted.bam L8_1_sorted.bam
```

#### refine MAGs?

```bash


```

#### dereplicate MAGs with dREP

```bash

```

#### filter MAGs for quality with checkM

```bash

```

#### calls SNVs with inStrain

First create a file for inStrain that has a column with the scaffold/contig name and the second with the MAG name. 

```bash
#!/usr/bin/bash
for f in *.fa
do
out="${f//.fa/_scaffolds.txt}"
genomes="${out//scaffolds.txt/genomes.txt}"
genscaf="${genomes//genomes.txt/genomesscaffolds.txt}"
grep MAG $f | sed 's/>//' > $out
cut -d_ -f1,2,3 $out > $genomes 
paste -d"\t" $out $genomes > $genscaf
done
cat *genomesscaffolds.txt > genome_scaffold.txt
rm *scaffolds.txt
rm *genomes.txt
```

```bash
cat *.fa > reseq_MAGS.fa
```

#### Call genes with prodigal for later use with inStrain.

```bash
prodigal -i reseq_MAGS.fa -d reseq_mag_genes.fna -a reseq_mag_genes.faa -p anon
```

#### Map original metagenomic reads from each timepoint to database of non-redundant MAGs

```bash
bowtie2-build reseq_MAGS.fa reseq_MAGS
```

```bash
#!/usr/bin/bash
for f in *_R1.fastq.gz; do bowtie2 -x reseq_MAGS -1 $f -2 ${f%*R1.fastq.gz}R2.fastq.gz -S ${f%*_R1.fastq.gz}.sam --threads 32
done
```

```bash
for f in *.sam; do samtools view -b $f > ${f%*sam}bam --threads 8; done
for f in *.bam; do samtools sort $f > ${f%*.bam}_sorted.bam --threads 8; done
for f in *sorted.bam; do samtools index $f; done
```

#### run inStrain

```bash
for f in *sorted.bam
do
inStrain profile $f reseq_MAGS.fa -o ${f$*.bam}_instrain_profile -p 64 -g reseq_mag_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.95 --skip_mm_profiling --min_genome_coverage 2 --min_freq 0
done
```

#### Annotation of MAGs with bakta

```bash
parallel -j 32 'bakta {} --db /mfs/ederrick/db --out {}_bakta --threads 8' ::: *.fasta
```

#### Also annotate with EggNOG to get COG categories for each gene because bakta doesn't output it nicely

```bash
emapper.py  -m diamond --itype CDS -i reseq_mag_genes.fna -o reseq_eggnog_genes --output_dir /mfs/ederrick/chapter_1 --cpu 36
```

#### calculate number of MAGs present at TP 1 and TP 3

subsample reads to pond with least sequences. Set same seed for R1 and R2.
could maybe just subsample bam file based on proportion of total reads.

```bash
seqkit sample -p ? -s 100 I4_1_R1.fastq.gz  -o subsamp_I4_1_R1.fastq.gz
seqkit sample -p ? -s 100 I8_1_R1.fastq.gz  -o subsamp_I8_1_R1.fastq.gz
seqkit sample -p ? -s 100 K1_1_R1.fastq.gz  -o subsamp_K1_1_R1.fastq.gz
seqkit sample -p ? -s 100 L2_1_R1.fastq.gz  -o subsamp_L2_1_R1.fastq.gz
seqkit sample -p ? -s 100 L3_1_R1.fastq.gz  -o subsamp_L3_1_R1.fastq.gz
seqkit sample -p ? -s 100 L4_1_R1.fastq.gz  -o subsamp_L4_1_R1.fastq.gz
seqkit sample -p ? -s 100 L6_1_R1.fastq.gz  -o subsamp_L6_1_R1.fastq.gz
seqkit sample -p ? -s 100 L7_1_R1.fastq.gz  -o subsamp_L7_1_R1.fastq.gz
seqkit sample -p ? -s 100 L8_1_R1.fastq.gz  -o subsamp_L8_1_R1.fastq.gz

seqkit sample -p ? -s 100 I4_1_R2.fastq.gz  -o subsamp_I4_1_R2.fastq.gz
seqkit sample -p ? -s 100 I8_1_R2.fastq.gz  -o subsamp_I8_1_R2.fastq.gz
seqkit sample -p ? -s 100 K1_1_R2.fastq.gz  -o subsamp_K1_1_R2.fastq.gz
seqkit sample -p ? -s 100 L2_1_R2.fastq.gz  -o subsamp_L2_1_R2.fastq.gz
seqkit sample -p ? -s 100 L3_1_R2.fastq.gz  -o subsamp_L3_1_R2.fastq.gz
seqkit sample -p ? -s 100 L4_1_R2.fastq.gz  -o subsamp_L4_1_R2.fastq.gz
seqkit sample -p ? -s 100 L6_1_R2.fastq.gz  -o subsamp_L6_1_R2.fastq.gz
seqkit sample -p ? -s 100 L7_1_R2.fastq.gz  -o subsamp_L7_1_R2.fastq.gz
seqkit sample -p ? -s 100 L8_1_R2.fastq.gz  -o subsamp_L8_1_R2.fastq.gz
```

get new number of subsampled reads

```bash
for file in *.gz; do seqkit stats $file; done > subsampled_read_count.txt
```

map subsampled reads to MAG database

```bash
for f in *_R1.fastq.gz; do bowtie2 -x reseq_MAGS -1 $f -2 ${f%*R1.fastq.gz}R2.fastq.gz --threads 32 -S ${f%*R1.fastq.gz}subsampled.sam; done
```

```bash
for f in *.sam; do samtools view -b $f > ${f%*sam}bam --threads 8; done
```

```bash
for f in *.bam; do samtools sort $f -o ${f%*.bam}sorted.bam --threads 8; done

```

```bash
for f in *subsampled_sorted.bam
do
cargo run -- genome -d coverm_MAGs -b $f --min-read-percent-identity 95 -o ${f%*subsampled_sorted.bam}_coverm -x fa -t 32 --min-covered-fraction 0 -m relative_abundance mean trimmed_mean covered_bases count reads_per_base 
done
```
