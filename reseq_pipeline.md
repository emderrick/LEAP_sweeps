### Pipeline for resequenced samples

#### rename fastq files

```bash
for f in *.gz; do mv $f ${f:0:13}${f:28}; done
for f in *__R1_001.fastq.gz; do mv $f ${f%*__R1_001.fastq.gz}_R1.fastq.gz; done
for f in *__R2_001.fastq.gz; do mv $f ${f%*__R2_001.fastq.gz}_R2.fastq.gz; done
for f in *001.fastq.gz; do mv $f ${f%*_001.fastq.gz}.fastq.gz; done
```

#### trim adaptors and remove low quality seqeunces

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate trimmomatic

for f in *R1.fastq.gz
do
trimmomatic PE -threads 64 -phred33 $f ${f%*R1.fastq.gz}R2.fastq.gz ${f%*R1.fastq.gz}QC_R1.fastq.gz ${f%*R1.fastq.gz}UP_R1.fastq.gz ${f%*R1.fastq.gz}QC_R2.fastq.gz ${f%*R1.fastq.gz}UP
_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:12 MINLEN:36
done
```

#### check sequence quality

```bash
parallel -j 18 'fastqc {}  --threads 16' ::: *.fastq.gz
```

#### coassemble all T1 samples

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate megahit
megahit -1 LEAP_META_01_QC_R1.fastq.gz,LEAP_META_02_QC_R1.fastq.gz,LEAP_META_03_QC_R1.fastq.gz,LEAP_META_04_QC_R1.fastq.gz,LEAP_META_05_QC_R1.fastq.gz,LEAP_META_06_QC_R1.fastq.gz,LEAP_META_07_QC_R1.fastq.gz,LEAP_META_08_QC_R1.fastq.gz,LEAP_META_09_QC_R1.fastq.gz \
-2 LEAP_META_01_QC_R2.fastq.gz,LEAP_META_02_QC_R2.fastq.gz,LEAP_META_03_QC_R2.fastq.gz,LEAP_META_04_QC_R2.fastq.gz,LEAP_META_05_QC_R2.fastq.gz,LEAP_META_06_QC_R2.fastq.gz,LEAP_META_07_QC_R2.fastq.gz,LEAP_META_08_QC_R2.fastq.gz,LEAP_META_09_QC_R2.fastq.gz \
-r LEAP_META_01_UP_R1.fastq.gz,LEAP_META_02_UP_R1.fastq.gz,LEAP_META_03_UP_R1.fastq.gz,LEAP_META_04_UP_R1.fastq.gz,LEAP_META_05_UP_R1.fastq.gz,LEAP_META_06_UP_R1.fastq.gz,LEAP_META_07_UP_R1.fastq.gz,LEAP_META_08_UP_R1.fastq.gz,LEAP_META_09_UP_R1.fastq.gz,LEAP_META_01_UP_R2.fastq.gz,LEAP_META_02_UP_R2.fastq.gz,LEAP_META_03_UP_R2.fastq.gz,LEAP_META_04_UP_R2.fastq.gz,LEAP_META_05_UP_R2.fastq.gz,LEAP_META_06_UP_R2.fastq.gz,LEAP_META_07_UP_R2.fastq.gz,LEAP_META_08_UP_R2.fastq.gz,LEAP_META_09_UP_R2.fastq.gz \
-o T1_coassembly --min-contig-len 1000 -t 64 -m 1000000000000
conda deactivate
``` 

```bash
seqkit stats -a T1_coassembly.fa > T1_coassembly_stats.txt
```

```bash
seqkit seq -m 2500 T1_coassembly.fa > T1_coassembly_2500.fa
seqkit stats -a T1_coassembly_2500.fa > T1_coassembly_2500_stats.txt
```

#### Map metagenomic reads from each pond at T1 to the co-assembly

```bash
bowtie2-build T1_coassembly_2500.fa T1_coassembly_2500 --threads 64
```

Try with local and end-to-end alignment and compare
 
```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_coassembly_2500 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --local --threads 16 | samtools sort -o {/QC_R1.fastq.gz/local_T1_coass.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_coassembly_2500 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/QC_R1.fastq.gz/T1_coassembly.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

Not a huge difference. Use end-to-end mapping

#### bin contigs with metabat2 (with bam files from TP 1)

```bash
docker pull metabat/metabat
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth T1_coassembly_depth.txt *.bam
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest metabat2 -i T1_coassembly_2500.fa -a T1_coassembly_depth.txt -o T1_bins/bin -m 2500 -t 48
```

#### check quality and filter and dereplicate MAGs

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
checkm lineage_wf all_bins all_MAGs_checkM -t 64 -x fa --tab_table -f all_MAGs_checkM.tsv --pplacer_threads 64

dRep dereplicate drep_T1_bins -g T1_bins/*.fa -l 500000 -comp 70 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### get stats of bins

```bash
seqkit stats -a * > bin_stats.txt
```

#### simplify fasta headers

```bash
for f in *.fa; do cut -f1 $f > ${f%*.fa}_fix.fa; done
for f in *_fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### Map T1 and T3 to the MAG database and calculate quick coverage stats

```bash
cat *.fa > T1_MAGs.fa
bowtie2-build T1_MAGs.fa T1_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/R1.fastq.gz/T1_MAGs.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

```bash
coverm genome -b *T1_MAGs.bam -d T1_dereplicated_MAGs -o T1_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa -- 
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
inStrain profile $f T1_MAGS.fa -o ${f%*.bam}_instrain_profile} -p 128 -g T1_MAG_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.93
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


### MISC things

#### run kraken2

```bash
kraken2 --db /mfs/databases/kraken-core-nt-dec-28-2024 T1_coassembly_2500.fa --threads 8 --output T1_coassembly_2500.output.txt --report T1_coassembly_2500.report.txt
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

for f in *R1.fastq.gz
do
kraken2 --db /mfs/databases/kraken-core-nt-dec-28-2024 --threads 64 --output ${f*%_QC_R1.fastq.gz}output.txt --report ${f*%_QC_R1.fastq.gz}report.txt --paired $f ${f*%R1.fastq.gz}R2.fastq.gz
done
```

