### Testing two coassembly methods for binning without anvi'o manual refinement 

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
-o T1_coassembly --min-contig-len 1000 -t 16 -m 1000000000000
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
/mfs/ederrick/miniconda3/envs/spades/bin/metaspades.py -1 T1_R1.fastq.gz -2 T1_R2.fastq.gz -o T1_metaspades_coassembly --only-assembler -m 1000 --threads 16
conda deactivate
```

#### compare assembly qualities. They end up very similar but MEGAHIT is much faster.

```bash
seqkit seq -m 1000 metaspades_contigs.fa > metaspades_contigs_filt.fa
seqkit stats -a T1_coassembly.fa > T1_coassembly_stats.txt
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

#### check quality and dereplicate MAGs with dREP and checkM. aftering filtering checkM results in R I have a list of good bins from both assemblies

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
checkm lineage_wf megahit_bins megahit_bin_quality -t 64 -x fa --tab_table -f megahit_checkM.txt --pplacer_threads 16
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate drep
checkm lineage_wf spades_bins spades_bin_quality -t 64 -x fa --tab_table -f spades_checkM.txt --pplacer_threads 16
```

```bash
cp $(<megahit_good_bins.txt) megahit_good_bins
for f in *.fa; do mv $f megahit_${f%*}; done
cp megahit_good_bins/* good_bins

cp $(<spades_good_bins.txt) spades_good_bins
for f in *.fa; do mv $f spades_${f%*}; done
cp spades_good_bins/* good_bins

checkm lineage_wf good_bins good_bin_quality -t 64 -x fa --tab_table -f good_checkM.txt --pplacer_threads 16
```

#### then do dereplication and filtering of all bins and see if bins from one method are chosen more. Spades bins are chosen more and have less contamination on average but assembly much slower)

```bash
for f in *.fa; do cp $f /mfs/ederrick/chapter_1/03_binning/megahit_${f%*}; done
for f in *.fa; do cp $f /mfs/ederrick/chapter_1/03_binning/spades_${f%*}; done
dRep dereplicate drep_all_bins -g all_bins/*.fa -l 500000 -comp 70 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### Map T1 and T3 to the MAG database to see how many MAGs are present in both at 1X (will hopefully be high enough to profile in reseq data)

```bash
cat *.fa > T1_MAGs.fa
bowtie2-build T1_MAGs.fa T1_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/T1_MAGs.bam} --write-index -@ 16' ::: *1_R1.fastq.gz
parallel -j 9 --plus 'bowtie2 -x T1_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/T1_MAGs.bam} --write-index -@ 16' ::: *3_R1.fastq.gz
```

```bash
conda create -n coverM
conda activate coverM
conda install coverM

coverm genome -b *.bam -d T1_dereplicated_MAGs -o T1_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 16 -x fa
```

#### try opposite way (assemble T3). Using MEGAHIT for speed.

#### Create a co-assembly of all ponds at time point 3

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate megahit
megahit -1 I4_3_R1.fastq.gz,I8_3_R1.fastq.gz,K1_3_R1.fastq.gz,L3_3_R1.fastq.gz,L4_3_R1.fastq.gz,L2_3_R1.fastq.gz,L6_3_R1.fastq.gz,L7_3_R1.fastq.gz,L8_3_R1.fastq.gz \
-2 I4_3_R2.fastq.gz,I8_3_R2.fastq.gz,K1_3_R2.fastq.gz,L3_3_R2.fastq.gz,L4_3_R2.fastq.gz,L2_3_R2.fastq.gz,L6_3_R2.fastq.gz,L7_3_R2.fastq.gz,L8_3_R2.fastq.gz \
-o T3_coassembly --min-contig-len 1000 --threads 64
conda deactivate
```

#### check assembly quality

```bash
seqkit stats -a T3_coassembly.fa > T3_coassembly_stats.txt
```

#### Map metagenomic reads from each pond at T3 to the co-assembly

```bash
bowtie2-build T3_coassembly.fa T3_coassembly --threads 64
```

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T3_coassembly -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/_R1.fastq.gz/.bam} --write-index -@ 8' ::: *_3_R1.fastq.gz

mkdir T3_bam
mv *.bam T3_bam
```

#### bin contigs with metabat2 (bam files in their respective directories)

```bash
jgi_summarize_bam_contig_depths --outputDepth T3_depth.txt *.bam
metabat2 -i T3_coassembly.fa -a T3_depth.txt -o T3_bins/bin -m 2500 -t 64
```

#### check quality and dereplicate MAGs with dRep and checkM

```bash
dRep dereplicate drep_T3_bins -g T3_bins/*.fa -l 500000 -comp 70 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

#### Map T1 and T3 to the MAG database to see how many MAGs are present in both at 1X (will hopefully be high enough to profile in reseq data)

```bash
cat *.fa > T3_MAGs.fa
bowtie2-build T3_MAGs.fa T3_MAGs --threads 64

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 9 --plus 'bowtie2 -x T3_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/T3_MAGs.bam} --write-index -@ 16' ::: *_1_R1.fastq.gz
parallel -j 9 --plus 'bowtie2 -x T3_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --local --threads 16 | samtools sort -o {/R1.fastq.gz/T3_MAGs.bam} --write-index -@ 16' ::: *_3_R1.fastq.gz
```

```bash
coverm genome -b *.bam -d T3_dereplicated_MAGs -o T3_MAGs_coverM.tsv -m mean variance covered_fraction relative_abundance -t 16 -x fa
```
