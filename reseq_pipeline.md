## Pipeline for resequenced samples

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

inStrain needs paired reads. Trim reads to keep both pairs when they overlap and use these reads with inStrain

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

#### check sequence quality

```bash
parallel -j 4 'fastqc {}  --threads 16' ::: *.fastq.gz
```

### Timepoint 1 MAG database

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
dRep dereplicate checkM_T1_bins -g all_T1_bins/*.fa -l 500000 -comp 70 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### get stats of bins

```bash
seqkit stats -a *.fa > MAG_stats.txt
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
parallel -j 4 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/T1_MAGs.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

```bash
coverm genome -b *T1_MAGs.bam -d /mfs/ederrick/chapter_1/04_MAGs/all_MAGs -o paired_T1_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa --output-format sparse
```

#### Annotate genes with prodigal

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate checkM
prodigal -i T1_MAGs.fa -d T1_MAG_genes.fna -a T1_MAG_genes.faa -o T1_MAG_genes.gbk -p meta
```

#### create .stb file for inStrain using script from dRep

```bash
conda activate drep
parse_stb.py --reverse -f all_MAGs/*.fa -o T1_MAGs.stb
```

#### calls SNVs with inStrain

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 4 --plus 'inStrain profile {} T1_MAGs.fa -o {/.bam/_inStrain} -p 24 -g T1_MAG_genes.fna -s T1_MAGs.stb --min_read_ani 0.92 --min_mapq 2 --min_genome_coverage 1' ::: *P_T1_MAGs.bam
```

### Using a 50% completion cutoff for MAGs

#### check quality with checkM and filter and check that MAGs are all < 95% ANI

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_T1_50_bins -g all_bins/*.fa -l 500000 -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### get stats of dereplicated and QC pass bins

```bash
seqkit stats -a *.fa > T1_MAGs_50_stats.txt
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

#### create .stb file for inStrain using script from dRep

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

### Build TP 1 and TP 2 MAG database

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

### TP 3 MAG database

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate megahit
megahit -1 LEAP_META_10_QC_R1.fastq.gz,LEAP_META_11_QC_R1.fastq.gz,LEAP_META_12_QC_R1.fastq.gz,LEAP_META_13_QC_R1.fastq.gz,LEAP_META_14_QC_R1.fastq.gz,LEAP_META_15_QC_R1.fastq.gz,LEAP_META_16_QC_R1.fastq.gz,LEAP_META_17_QC_R1.fastq.gz,LEAP_META_18_QC_R1.fastq.gz \
-2 LEAP_META_10_QC_R2.fastq.gz,LEAP_META_11_QC_R2.fastq.gz,LEAP_META_12_QC_R2.fastq.gz,LEAP_META_13_QC_R2.fastq.gz,LEAP_META_14_QC_R2.fastq.gz,LEAP_META_15_QC_R2.fastq.gz,LEAP_META_16_QC_R2.fastq.gz,LEAP_META_17_QC_R2.fastq.gz,LEAP_META_18_QC_R2.fastq.gz \
-r LEAP_META_10_UP_R1.fastq.gz,LEAP_META_11_UP_R1.fastq.gz,LEAP_META_12_UP_R1.fastq.gz,LEAP_META_13_UP_R1.fastq.gz,LEAP_META_14_UP_R1.fastq.gz,LEAP_META_15_UP_R1.fastq.gz,LEAP_META_16_UP_R1.fastq.gz,LEAP_META_17_UP_R1.fastq.gz,LEAP_META_18_UP_R1.fastq.gz,LEAP_META_10_UP_R2.fastq.gz,LEAP_META_11_UP_R2.fastq.gz,LEAP_META_12_UP_R2.fastq.gz,LEAP_META_13_UP_R2.fastq.gz,LEAP_META_14_UP_R2.fastq.gz,LEAP_META_15_UP_R2.fastq.gz,LEAP_META_16_UP_R2.fastq.gz,LEAP_META_17_UP_R2.fastq.gz,LEAP_META_18_UP_R2.fastq.gz \
-o T3_coassembly --min-contig-len 1000 -t 128 -m 1000000000000
conda deactivate
```

```bash
seqkit stats -a T3_coassembly.fa > T3_coassembly_stats.txt
seqkit seq -m 2500 T3_coassembly.fa > T3_coassembly_2500.fa
seqkit stats -a T3_coassembly_2500.fa > T3_coassembly_2500_stats.txt
```

```bash
bowtie2-build T3_coassembly_2500.fa T3_coassembly_2500 --threads 64
```

with T3 reads do

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T3_coassembly_2500 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/QC_R1.fastq.gz/T3_coassembly.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

#### reconstruct MAGs from T3

```bash
docker pull metabat/metabat
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth T3_coassembly_depth.txt *.bam
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest metabat2 -i T3_coassembly_2500.fa -a T3_coassembly_depth.txt -o T3_bins/bin -m 2500 -t 48
```

```bash
for f in *.fa; do cut -f1 $f > ${f%*.fa}_fix.fa; done
for f in *_fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### filter T3 bins

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_T3_bins -g T3_bins/*.fa -l 500000 -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64

seqkit stats -a * > T3_MAG_stats.txt
```

#### map back to MAGs

```bash
cat *.fa > T3_MAGs.fa
bowtie2-build T3_MAGs.fa T3_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T3_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/R1.fastq.gz/T3_MAGs.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

#### estimate quick coverage of MAGs
```bash
coverm genome -b *T3_MAGs.bam -d T3_MAGs -o T3_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 64 -x fa --output-format sparse
```

### dereplicate T1 MAG database and T3 MAG database to get list of duplicate species. Then manually create non-overlapping database while selecting T1 MAGs over T3 MAGs.

in T1 dereplicated MAGs

```bash
for f in *.fa; do cp $f T1_$f; done
```

in T3 dereplicated MAGs

```bash
for f in *.fa; do cp $f T3_$f; done
```

```bash
dRep dereplicate drep_T1_T3_MAGs -g T1_T3_MAGs/*.fa -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

### community composition stuff

```bash
seqkit stats -j 18 *.gz -a -T > paired_read_stats.tsv
```

#### downsample fastq reads to lowest sample depth

```bash
parallel -j 18 --plus 'seqtk sample -s100 {} 87942340 | gzip > {/R1.fastq.gz/sub_R1.fastq.gz}' ::: *R1.fastq.gz
parallel -j 18 --plus 'seqtk sample -s100 {} 87942340 | gzip > {/R2.fastq.gz/sub_R2.fastq.gz}' ::: *R2.fastq.gz
```

```bash
seqkit stats -j 18 *sub* -a -T > subsamp_read_stats.tsv
```

#### try metaphlan

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate metaphlan

parallel -j 9 --plus 'metaphlan {},{/R1_fastq.gz/R2.fastq.gz} --bowtie2out {/_R1.fastq.gz/.bowtie2.bz2} --nproc 24 --input_type fastq --unclassified_estimation -o {/R1.fastq.gz/metaphlan.txt}' ::: *sub_R1.fastq.gz
```

```bash
merge_metaphlan_tables.py *metaphlan.txt > merged_metaphlan_abundance_table.tsv
```

#### try kraken2

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

for f in *sub_R1.fastq.gz
do
kraken2 --db /mfs/databases/kraken-core-nt-dec-28-2024 --threads 128 --output ${f%*R1.fastq.gz}kraken_output.txt --report ${f%*R1.fastq.gz}kraken_report.txt --paired $f ${f%*R1.fastq.gz}R2.fastq.gz
done
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bracken

for f in *report.txt
do
bracken -d /mfs/databases/kraken-core-nt-dec-28-2024 -i $f -o ${f%*_kraken_report.txt}.bracken -r 100 -l S -t 10
done
```
