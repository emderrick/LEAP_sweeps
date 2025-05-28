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
parallel -j 4 'fastqc {}  --threads 16' ::: *.fastq.gz
```

### using a MAG database built from timepoint 1 only

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

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_coassembly_2500 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/QC_R1.fastq.gz/T1_coassembly.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

#### bin contigs with metabat2 (with bam files from TP 1)

```bash
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth T1_coassembly_depth.txt *.bam
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest metabat2 -i T1_coassembly_2500.fa -a T1_coassembly_depth.txt -o T1_bins/bin -m 2500 -t 48
```

#### simplify fasta headers

```bash
cp T1_bins/* all_bins
for f in *.fa; do cut -f1 $f > ${f%*.fa}_fix.fa; done
for f in *_fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### check quality with checkM and filter and check that MAGs are all < 95% ANI

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_T1_bins -g all_bins/*.fa -l 500000 -comp 70 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### get stats of bins

```bash
seqkit stats -a * > MAG_stats.txt
```

#### Map T1 and T3 to the MAG database and calculate quick coverage stats

```bash
cat *.fa > T1_MAGs.fa
bowtie2-build T1_MAGs.fa T1_MAGs --threads 128
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/R1.fastq.gz/T1_MAG.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

```bash
coverm genome -b *T1_MAG.bam -d all_MAGs -o T1_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa --output-format sparse
```

#### Annotate genes with prodigal

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate checkM
prodigal -i T1_MAGs.fa -d T1_MAG_genes.fna -a T1_MAG_genes.faa -o T1_MAG_genes.gbk -p meta
```

### create .stb file for inStrain using script from dRep

```bash
conda activate drep
parse_stb.py --reverse -f all_MAGs/*.fa -o T1_MAGs.stb
```

#### calls SNVs with inStrain

inStrain needs paired reads. Go back and trim reads to keep both pairs when they overlap and remap with those.

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate trimmomatic

for f in *R1.fastq.gz
do
trimmomatic PE -threads 64 -phred33 $f ${f%*R1.fastq.gz}R2.fastq.gz ${f%*R1.fastq.gz}P_R1.fastq.gz ${f%*R1.fastq.gz}NP_R1.fastq.gz ${f%*R1.fastq.gz}P_R2.fastq.gz ${f%*R1.fastq.gz}NP
_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:12 MINLEN:36
done
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 4 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/T1_MAGs.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

```bash
coverm genome -b *T1_MAGs.bam -d /mfs/ederrick/chapter_1/04_MAGs/all_MAGs -o paired_T1_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa --output-format sparse
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 4 --plus 'inStrain profile {} T1_MAGs.fa -o {/.bam/_inStrain} -p 24 -g T1_MAG_genes.fna -s T1_MAGs.stb --min_read_ani 0.92 --min_mapq 2 --min_genome_coverage 1' ::: *P_T1_MAGs.bam
```

### using a MAG database built from both timepoints

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate megahit
megahit -1 LEAP_META_01_QC_R1.fastq.gz,LEAP_META_02_QC_R1.fastq.gz,LEAP_META_03_QC_R1.fastq.gz,LEAP_META_04_QC_R1.fastq.gz,LEAP_META_05_QC_R1.fastq.gz,LEAP_META_06_QC_R1.fastq.gz,LEAP_META_07_QC_R1.fastq.gz,LEAP_META_08_QC_R1.fastq.gz,LEAP_META_09_QC_R1.fastq.gz,LEAP_META_10_QC_R1.fastq.gz,LEAP_META_11_QC_R1.fastq.gz,LEAP_META_12_QC_R1.fastq.gz,LEAP_META_13_QC_R1.fastq.gz,LEAP_META_14_QC_R1.fastq.gz,LEAP_META_15_QC_R1.fastq.gz,LEAP_META_16_QC_R1.fastq.gz,LEAP_META_17_QC_R1.fastq.gz,LEAP_META_18_QC_R1.fastq.gz \
-2 LEAP_META_01_QC_R2.fastq.gz,LEAP_META_02_QC_R2.fastq.gz,LEAP_META_03_QC_R2.fastq.gz,LEAP_META_04_QC_R2.fastq.gz,LEAP_META_05_QC_R2.fastq.gz,LEAP_META_06_QC_R2.fastq.gz,LEAP_META_07_QC_R2.fastq.gz,LEAP_META_08_QC_R2.fastq.gz,LEAP_META_09_QC_R2.fastq.gz,LEAP_META_10_QC_R2.fastq.gz,LEAP_META_11_QC_R2.fastq.gz,LEAP_META_12_QC_R2.fastq.gz,LEAP_META_13_QC_R2.fastq.gz,LEAP_META_14_QC_R2.fastq.gz,LEAP_META_15_QC_R2.fastq.gz,LEAP_META_16_QC_R2.fastq.gz,LEAP_META_17_QC_R2.fastq.gz,LEAP_META_18_QC_R2.fastq.gz \
-r LEAP_META_01_UP_R1.fastq.gz,LEAP_META_02_UP_R1.fastq.gz,LEAP_META_03_UP_R1.fastq.gz,LEAP_META_04_UP_R1.fastq.gz,LEAP_META_05_UP_R1.fastq.gz,LEAP_META_06_UP_R1.fastq.gz,LEAP_META_07_UP_R1.fastq.gz,LEAP_META_08_UP_R1.fastq.gz,LEAP_META_09_UP_R1.fastq.gz,LEAP_META_10_UP_R1.fastq.gz,LEAP_META_11_UP_R1.fastq.gz,LEAP_META_12_UP_R1.fastq.gz,LEAP_META_13_UP_R1.fastq.gz,LEAP_META_14_UP_R1.fastq.gz,LEAP_META_15_UP_R1.fastq.gz,LEAP_META_16_UP_R1.fastq.gz,LEAP_META_17_UP_R1.fastq.gz,LEAP_META_18_UP_R1.fastq.gz,LEAP_META_01_UP_R2.fastq.gz,LEAP_META_02_UP_R2.fastq.gz,LEAP_META_03_UP_R2.fastq.gz,LEAP_META_04_UP_R2.fastq.gz,LEAP_META_05_UP_R2.fastq.gz,LEAP_META_06_UP_R2.fastq.gz,LEAP_META_07_UP_R2.fastq.gz,LEAP_META_08_UP_R2.fastq.gz,LEAP_META_09_UP_R2.fastq.gz,LEAP_META_10_UP_R2.fastq.gz,LEAP_META_11_UP_R2.fastq.gz,LEAP_META_12_UP_R2.fastq.gz,LEAP_META_13_UP_R2.fastq.gz,LEAP_META_14_UP_R2.fastq.gz,LEAP_META_15_UP_R2.fastq.gz,LEAP_META_16_UP_R2.fastq.gz,LEAP_META_17_UP_R2.fastq.gz,LEAP_META_18_UP_R2.fastq.gz \
-o full_coassembly --min-contig-len 1000 -t 128 -m 1000000000000
conda deactivate
```

```bash
seqkit stats -a full_coassembly.fa > full_coassembly_stats.txt
```

```bash
seqkit seq -m 2500 full_coassembly.fa > full_coassembly_2500.fa
seqkit stats -a full_coassembly_2500.fa > full_coassembly_2500_stats.txt
```

#### Map metagenomic reads from each pond to the co-assembly

```bash
bowtie2-build full_coassembly_2500.fa full_coassembly_2500 --threads 64 --large-index
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 4 --plus 'bowtie2 -x full_coassembly_2500 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/QC_R1.fastq.gz/full_coassembly.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

#### bin contigs with metabat2 (with bam files from TP 1 and TP 3)

```bash
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth full_coassembly_depth.txt *.bam
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest metabat2 -i full_coassembly_2500.fa -a full_coassembly_depth.txt -o full_bins/bin -m 2500
```

#### simplify fasta headers

```bash
cp -r full_bins combined_bins
for f in *.fa; do cut -f1 $f > ${f%*.fa}_fix.fa; done
for f in *_fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### check quality with checkM and filter and check that MAGs are all < 95% ANI

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_combined_bins -g combined_bins/*.fa -l 500000 -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### get stats of bins

```bash
seqkit stats -a * > combined_MAG_stats.txt
```

#### Map T1 and T3 to the MAG database and calculate quick coverage stats

```bash
cat *.fa > combined_MAGs.fa
bowtie2-build combined_MAGs.fa combined_MAGs --threads 128
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 4 --plus 'bowtie2 -x combined_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/combined_MAGs.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

```bash
coverm genome -b *combined_MAGs.bam -d combined_MAGs -o combined_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa --output-format sparse
```

```bash
prodigal -i combined_MAGs.fa -d combined_MAG_genes.fna -a combined_MAG_genes.faa -o combined_MAG_genes.gbk -p meta
```

```bash
conda activate drep
parse_stb.py --reverse -f combined_MAGs/*.fa -o combined_MAGs.stb
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 4 --plus 'inStrain profile {} combined_MAGs.fa -o {/.bam/_inStrain} -p 24 -g combined_MAG_genes.fna -s combined_MAGs.stb --min_read_ani 0.92 --min_mapq 2 --min_genome_coverage 1' ::: *P_combined_MAGs.bam
```

### Using a 50% completion cutoff for MAGs

#### check quality with checkM and filter and check that MAGs are all < 95% ANI

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_T1_50_bins -g all_bins/*.fa -l 500000 -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### get stats of bins

```bash
seqkit stats -a * > T1_MAGs_50_stats.txt
```

#### Map T1 and T3 to the MAG database and calculate quick coverage stats

```bash
cat *.fa > T1_50_MAGs.fa
bowtie2-build T1_50_MAGs.fa T1_50_MAGs --threads 64
```
```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_50_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/T1_50_MAGs.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

```bash
coverm genome -b *T1_50_MAGs.bam -d /mfs/ederrick/chapter_1/03_T1_binning/T1_50_MAG_stuff/T1_50_MAGs -o T1_50_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa --output-format sparse
```

#### Annotate genes with prodigal

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate checkM
prodigal -i T1_50_MAGs.fa -d T1_50_MAG_genes.fna -a T1_50_MAG_genes.faa -o T1_50_MAG_genes.gbk -p meta
```

### create .stb file for inStrain using script from dRep

```bash
conda activate drep
parse_stb.py --reverse -f T1_50_MAGs/*.fa -o T1_50_MAGs.stb
```

#### calls SNVs with inStrain

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain
parallel -j 4 --plus 'inStrain profile {} T1_50_MAGs.fa -o {/P_T1_50_MAGs.bam/T1_50_inStrain} -p 24 -g T1_50_MAG_genes.fna -s T1_50_MAGs.stb --min_read_ani 0.92 --min_mapq 2 --min_genome_coverage 1' ::: *T1_50_MAGs.bam
```

#### run kraken2 on paired reads

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

for f in *_P_R1.fastq.gz
do
kraken2 --db /mfs/databases/kraken-core-nt-dec-28-2024 --threads 24 --output ${f%*P_R1.fastq.gz}kraken_output.txt --report ${f%*P_R1.fastq.gz}kraken_report.txt --paired $f ${f%*R1.fastq.gz}R2.fastq.gz
done
```

#### classify MAGs with GTDB

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate gtdbtk

export GTDBTK_DATA_PATH=/mfs/ederrick/release220
gtdbtk classify_wf --genome_dir T1_50_MAGs --pplacer_cpus 64 --cpus 64 --extension fa --out_dir GTDB_T1_50_MAGs --skip_ani_screen
```

#### Annotation of MAGs with bakta

```bash
parallel -j 32 'bakta {} --db /mfs/ederrick/db --out {}_bakta --threads 8' ::: *.fa
```

