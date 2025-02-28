#### Annotate genes with prodigal

```bash
prodigal -i T1_MAGs.fa -d T1_mag_genes.fna -a T1_mag_genes.faa -p anon
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

```bash
conda create -n instrain
conda activate instrain
conda install instrain

for f in *MAGs.bam
do
inStrain profile $f T1_MAGS.fa -o ${f%*.bam}_instrain_profile} -p 128 -g T1_mag_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.95
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
