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

#### Coassemble all T1 samples

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

#### bin contigs with metabat2 (with bam files from TP 1)

Map metagenomic reads from each pond at T1 to the co-assembly

```bash
bowtie2-build T1_coassembly_2500.fa T1_coassembly_2500 --threads 64
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_coassembly_2500 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} -U {/QC_R1.fastq.gz/UP_R1.fastq.gz},{/QC_R2.fastq.gz/UP_R2.fastq.gz} --threads 16 | samtools sort -o {/QC_R1.fastq.gz/T1_coassembly.bam} --write-index -@ 16' ::: *QC_R1.fastq.gz
```

```bash
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth T1_coassembly_depth.txt *.bam
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest metabat2 -i T1_coassembly_2500.fa -a T1_coassembly_depth.txt -o T1_bins/bin -m 2500 -t 48
```

simplify fasta headers

```bash
cp T1_bins/* all_T1_bins
for f in *.fa; do cut -f1 $f > ${f%*.fa}_fix.fa; done
for f in *_fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### filter T1 bins

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_T1_50_bins -g all_T1_bins/*.fa -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

```bash
seqkit stats -a *.fa > T1_MAGs_50_stats.txt
```

### Timepoint 3 MAG database

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

#### bin contigs

```bash
docker pull metabat/metabat
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest jgi_summarize_bam_contig_depths --outputDepth T3_coassembly_depth.txt *.bam
docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat/metabat:latest metabat2 -i T3_coassembly_2500.fa -a T3_coassembly_depth.txt -o T3_bins/bin -m 2500 -t 48
```

```bash
cp T3_bins/* all_T3_bins
for f in *.fa; do cut -f1 $f > ${f%*.fa}_fix.fa; done
for f in *_fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### filter T3 bins

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
dRep dereplicate checkM_T3_bins -g all_T3_bins/*.fa -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64

seqkit stats -a * > T3_MAG_stats.txt
```

### dereplicate T1 MAG database and T3 MAG database

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

Then I went through drep_T1_T3_MAGs/figures/secondary_clustering_dendrograms.pdf and made a list of all the T3 MAGs that clustered with a T1 MAG

Then remove those MAGs to make a non-redundant MAG database

```bash
cp -r T1_T3_MAGs nonred_T1_T3_MAGs
cat to_remove.txt | xargs -I {} mv {} ..
```

```bash
seqkit stats -a *.fa > nonred_T1_T3_stats.txt
```

Then add T1 or T3 to the start of contig names to make sure they are all unique

```bash
for f in *T1*; do sed 's/>/>T1_/' $f > ${f%*.fa}_fix.fa; done
for f in *T3*; do sed 's/>/>T3_/' $f > ${f%*.fa}_fix.fa; done
for f in *fix.fa; do mv $f ${f%*_fix.fa}.fa; done
```

#### profile MAGs with inStrain

map to nonredundant MAGs

```bash
cat *.fa > nonred_T1_T3.fa
bowtie2-build nonred_T1_T3.fa nonred_T1_T3 --threads 128

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 6 --plus 'bowtie2 -x nonred_T1_T3 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/nonred_T1_T3.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

Annotate genes with prodigal

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate checkM
prodigal -i nonred_T1_T3.fa -d nonred_T1_T3_genes.fna -a nonred_T1_T3_genes.faa -o nonred_T1_T3_genes.gbk -p meta
```

create .stb file for inStrain using script from dRep

```bash
conda activate drep
parse_stb.py --reverse -f nonred_T1_T3_MAGs/*.fa -o nonred_T1_T3_MAGs.stb
```

calls SNVs with inStrain

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 6 --plus 'inStrain profile {} nonred_T1_T3.fa -o {/.bam/_inStrain} -p 24 -g nonred_T1_T3_genes.fna -s nonred_T1_T3_MAGs.stb --min_read_ani 0.92 --min_mapq 2 --min_genome_coverage 
1' ::: *nonred_T1_T3.bam
```

get depth of each position

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 18 'samtools depth -a -Q 1 LEAP_META_{}_P_nonred_T1_T3.bam -o LEAP_META_{}_nonred_depth.txt' ::: {01..18}
```

```bash
Rscript get_depth.R
```

#### annotate genes with eggnog to get COG categories

```bash
emapper.py -m diamond --itype CDS -i T1_50_MAG_genes.fna -o eggnog_genes --output_dir /mfs/ederrick/chapter_1/06_inStrain/T1_50_inStrain/ --cpu 72
```

#### also annotate with bakta?

```bash
for f in *.fa; do bakta --db /mfs/ederrick/db $f ${f%*.fa}_bakta --threads 8; done
```

### community composition of ponds

#### use reads with stricter filtering. Subsample first.

find lowest read depth

```bash
seqkit stats -j 18 *.gz -a -T > QC_read_stats.tsv
```

subsample

```bash
parallel -j 18 --plus 'seqtk sample -s100 {} 68107586 | gzip > {/R1.fastq.gz/sub_R1.fastq.gz}' ::: *QC_R1.fastq.gz
parallel -j 18 --plus 'seqtk sample -s100 {} 68107586 | gzip > {/R2.fastq.gz/sub_R2.fastq.gz}' ::: *QC_R2.fastq.gz
parallel -j 18 --plus 'seqtk sample -s100 {} 20194431 | gzip > {/R1.fastq.gz/sub_R1.fastq.gz}' ::: *UP_R1.fastq.gz
```

#### run metaphlan with paired reads and unpaired

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate metaphlan

parallel -j 18 --plus 'metaphlan {},{/R1_fastq.gz/R2.fastq.gz},{/QC_sub_R1.fastq.gz/UP_sub_R1.fastq.gz} --bowtie2out {/_R1.fastq.gz/.bowtie2.bz2} --nproc 12 --input_type fastq --unclassified_estimation -o {/R1.fastq.gz/metaphlan.txt}' ::: *QC_sub_R1.fastq.gz
```

```bash
merge_metaphlan_tables.py *metaphlan.txt > QC_merged_metaphlan_abundance_table.tsv
```

see without unclassified 

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate metaphlan

parallel -j 6 --plus 'metaphlan {},{/R1_fastq.gz/R2.fastq.gz},{/QC_sub_R1.fastq.gz/UP_sub_R1.fastq.gz} --bowtie2out {/R1.fastq.gz/nounclass.bowtie2.bz2} --nproc 12 --input_type fastq -o {/R1.fastq.gz/no_unclass_metaphlan.txt}' ::: *QC_sub_R1.fastq.gz
```

```bash
merge_metaphlan_tables.py *metaphlan.txt > QC_merged_metaphlan_abundance_table.tsv
```

#### one last try kraken with the QC reads. Does not support the SE library with the paired.

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

for f in *QC_sub_R1.fastq.gz
do
kraken2 --db /mfs/databases/kraken-core-nt-dec-28-2024 --confidence 0.1 --threads 72 --output ${f%*R1.fastq.gz}kraken_output.txt --report ${f%*R1.fastq.gz}kraken_report.txt --paired $f ${f%*R1.fastq.gz}R2.fastq.gz
done
```

extract only bacterial reads

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

parallel -j 6 'python /mfs/ederrick/miniconda3/envs/kraken2/bin/extract_kraken_reads.py -k LEAP_META_{}_QC_sub_kraken_output.txt -r LEAP_META_{}_QC_sub_kraken_report.txt -s LEAP_META_{}_QC_sub_R1.fastq.gz -s2 LEAP_META_{}_QC_sub_R2.fastq.gz -o LEAP_META_{}_QC_bac_R1.fastq.gz -o2 LEAP_META_{}_QC_bac_R2.fastq.gz -t 2 --include-children' ::: {01..18}
```

run bracken

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bracken

for f in *report.txt
do
bracken -d /mfs/databases/kraken-core-nt-dec-28-2024 -i $f -o ${f%*kraken_report.txt}phylum.bracken -r 100 -l P -t 10
done
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bracken

for f in *report.txt
do
bracken -d /mfs/databases/kraken-core-nt-dec-28-2024 -i $f -o ${f%*kraken_report.txt}class.bracken -r 100 -l C -t 10
done
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bracken

for f in *report.txt
do
bracken -d /mfs/databases/kraken-core-nt-dec-28-2024 -i $f -o ${f%*kraken_report.txt}genus.bracken -r 100 -l G -t 10
done
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bracken

for f in *report.txt
do
bracken -d /mfs/databases/kraken-core-nt-dec-28-2024 -i $f -o ${f%*kraken_report.txt}species.bracken -r 100 -l S -t 10
done
```

```bash
for f in *.bracken; do  python /mfs/ederrick/miniconda3/envs/kraken2/bin/alpha_diversity.py -f $f -a Si; done
```


#### try with no confidence threshold

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

parallel -j 4 'kraken2 --db /mfs/databases/kraken-core-nt-dec-28-2024 --threads 24 --output LEAP_META_{}_QC_sub_kraken_output.txt --report LEAP_META_{}_QC_sub_kraken_report.txt --paired LEAP_META_{}_QC_sub_R1.fastq.gz LEAP_META_{}_QC_sub_R2.fastq.gz' ::: {01..18}
```

extract only bacterial reads

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate kraken2

parallel -j 6 'python /mfs/ederrick/miniconda3/envs/kraken2/bin/extract_kraken_reads.py -k LEAP_META_{}_QC_sub_kraken_output.txt -r LEAP_META_{}_QC_sub_kraken_report.txt -s LEAP_META_{}_QC_sub_R1.fastq.gz -s2 LEAP_META_{}_QC_sub_R2.fastq.gz -o LEAP_META_{}_QC_bac_R1.fastq.gz -o2 LEAP_META_{}_QC_bac_R2.fastq.gz -t 2 --include-children' ::: {01..18}
```


