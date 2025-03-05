### Pipeline for resequenced samples

#### trim adaptors and remove low quality seqeunces

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

#### coassemble all T1 samples

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate megahit
megahit -1 I4_1_R1.fastq.gz,I8_1_R1.fastq.gz,K1_1_R1.fastq.gz,L3_1_R1.fastq.gz,L4_1_R1.fastq.gz,L2_1_R1.fastq.gz,L6_1_R1.fastq.gz,L7_1_R1.fastq.gz,L8_1_R1.fastq.gz \
-2 I4_1_R2.fastq.gz,I8_1_R2.fastq.gz,K1_1_R2.fastq.gz,L3_1_R2.fastq.gz,L4_1_R2.fastq.gz,L2_1_R2.fastq.gz,L6_1_R2.fastq.gz,L7_1_R2.fastq.gz,L8_1_R2.fastq.gz \
-o T1_coassembly --min-contig-len 1000 -t 64 -m 1000000000000
conda deactivate
```

```bash
seqkit stats -a T1_coassembly.fa > T1_coassembly_stats.txt
```

#### Map metagenomic reads from each pond at T1 to the co-assembly

```bash
bowtie2-build T1_coassembly.fa T1_coassembly --threads 64
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_coassembly -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/_R1.fastq.gz/T1_coassembly.bam} --write-index -@ 16' ::: *1_R1.fastq.gz
```

#### bin contigs with metabat2 (bam files in their respective directories)

```bash
jgi_summarize_bam_contig_depths --outputDepth T1_coassembly_depth.txt *.bam
metabat2 -i T1_coassembly.fa -a T1_coassembly_depth.txt -o T1_bins/bin -m 2500 -t 64
```

#### check quality and then filter and dereplicate MAGs

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
checkm lineage_wf T1_bins T1_bins_quality -t 64 -x fa --tab_table -f T1_bins_checkM.txt --pplacer_threads 32
```

```bash
dRep dereplicate drep_T1_bins -g T1_bins/*.fa -l 500000 -comp 70 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### Map T1 and T3 to the MAG database and calculate quick coverage stats

```bash
cd T1_dereplicated_MAGs
cat *.fa > T1_MAGs.fa
bowtie2-build T1_MAGs.fa T1_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/T1_MAGs.bam} --write-index -@ 16' ::: *R1.fastq.gz
```

```bash
coverm genome -b *T1_MAGs.bam -d T1_dereplicated_MAGs -o T1_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 16 -x fa
```

#### Annotate genes with prodigal

```bash
prodigal -i T1_MAGs.fa -d T1_MAG_genes.fna -a T1_MAG_genes.faa -p anon
```

#### calls SNVs with inStrain

```bash
conda create -n instrain
conda activate instrain
conda install instrain

for f in *T1_MAGs.bam
do
inStrain profile $f T1_MAGS.fa -o ${f%*.bam}_instrain_profile} -p 128 -g T1_MAG_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.95
done
```

#### Annotation of MAGs with bakta

```bash
parallel -j 32 'bakta {} --db /mfs/ederrick/db --out {}_bakta --threads 8' ::: *.fasta
```

#### calculate RA of MAGs present at TP 1 and TP 3

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

seqkit sample -p ? -s 100 I4_3_R2.fastq.gz  -o subsamp_I4_3_R2.fastq.gz
seqkit sample -p ? -s 100 I8_3_R2.fastq.gz  -o subsamp_I8_3_R2.fastq.gz
seqkit sample -p ? -s 100 K1_3_R2.fastq.gz  -o subsamp_K1_3_R2.fastq.gz
seqkit sample -p ? -s 100 L2_3_R2.fastq.gz  -o subsamp_L2_3_R2.fastq.gz
seqkit sample -p ? -s 100 L3_3_R2.fastq.gz  -o subsamp_L3_3_R2.fastq.gz
seqkit sample -p ? -s 100 L4_3_R2.fastq.gz  -o subsamp_L4_3_R2.fastq.gz
seqkit sample -p ? -s 100 L6_3_R2.fastq.gz  -o subsamp_L6_3_R2.fastq.gz
seqkit sample -p ? -s 100 L7_3_R2.fastq.gz  -o subsamp_L7_3_R2.fastq.gz
seqkit sample -p ? -s 100 L8_3_R2.fastq.gz  -o subsamp_L8_3_R2.fastq.gz
```

confirm number of subsampled reads

```bash
seqkit stats *subsamp* > subsampled_read_count.txt
```

map subsampled reads to MAG database

```bash
parallel -j 9 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/subsampled.bam} --write-index -@ 16' ::: *R1.fastq.gz
```

```bash
coverm genome -b *subsampled.bam -d T1_dereplicated_MAGs --min-read-percent-identity 95 -o subsampled_reads_coverM.tsv -m mean variance covered_bases covered_fraction relative_abundance -t 16 -x fa
```
