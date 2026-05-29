### Try with adding in Timepoint 3 MAGs

```bash
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

with T3 reads do

```bash
conda activate bowtie2
bowtie2-build T3_coassembly_2500.fa T3_coassembly_2500 --threads 64

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
conda activate bowtie2
cat *.fa > nonred_T1_T3.fa
bowtie2-build nonred_T1_T3.fa nonred_T1_T3 --threads 128

parallel -j 6 --plus 'bowtie2 -x nonred_T1_T3 -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/nonred_T1_T3.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

Annotate genes with prodigal

```bash
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
conda activate instrain

parallel -j 6 --plus 'inStrain profile {} nonred_T1_T3.fa -o {/.bam/_inStrain} -p 24 -g nonred_T1_T3_genes.fna -s nonred_T1_T3_MAGs.stb --min_read_ani 0.92 --min_mapq 2 --min_genome_coverage 
1' ::: *nonred_T1_T3.bam
```
