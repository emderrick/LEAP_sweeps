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

map to T1 MAGs

```bash
cat *.fa > T1_50_MAGs.fa
bowtie2-build T1_50_MAGs.fa T1_50_MAGs --threads 128

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 6 --plus 'bowtie2 -x T1_50_MAGs -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/T1_50_MAGs.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

Annotate genes with prodigal

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate checkM
prodigal -i T1_50_MAGs.fa -d T1_50_MAG_genes.fna -a T1_50_MAG_genes.faa -o T1_50_MAG_genes.gbk -p meta
```

create .stb file for inStrain using script from dRep

```bash
conda activate drep
parse_stb.py --reverse -f T1_50_MAGs/*.fa -o T1_50_MAGs.stb
```

calls SNVs with inStrain

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 18 --plus 'inStrain profile {} T1_50_MAGs.fa -o {/P_T1_50_MAGs.bam/T1_95_inStrain} -p 12 -g T1_50_MAG_genes.fna -s T1_50_MAGs.stb  --min_read_ani 0.92 --min_mapq 1' ::: *T1_50_MAGs.bam
```

#### annotate genes with eggnog to get COG categories

```bash
emapper.py -m diamond --itype CDS -i T1_50_MAG_genes.fna -o eggnog_genes --output_dir /mfs/ederrick/chapter_1/06_inStrain/T1_50_inStrain/ --cpu 72
```

#### also annotate with bakta?

```bash
for f in *.fa; do bakta --db /mfs/ederrick/db $f ${f%*.fa}_bakta --threads 8; done
```

### Try with adding in Timepoint 3 MAGs

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

### MAGs kind of suck. run anvi'o pipeline to manually check and refine MAGs

```bash
cut -d" " -f1 T1_coassembly_2500.fa > T1_coassembly_anvio.fa
anvi-gen-contigs-database -f T1_coassembly_anvio.fa -o T1_coassembly_anvio.db -T 74
```

```bash
bowtie2-build T1_coassembly_anvio.fa T1_coassembly_anvio --threads 72
parallel -j 9 'bowtie2 -x T1_coassembly_anvio -1 LEAP_META_{}_QC_R1.fastq.gz -2 LEAP_META_{}_QC_R2.fastq.gz -U LEAP_META_{}_UP_R1.fastq.gz,LEAP_META_{}_UP_R2.fastq.gz --threads 12 | samtools sort -o LEAP_META_{}_T1_anvio.bam --write-index -@ 12' ::: {01..09}
```

```bash
for f in *T1_anvio.bam; do anvi-profile -i $f -c T1_coassembly_anvio.db -o ${f%*.bam}_profile --min-coverage-for-variability 5 -T 72; done
```

```bash
anvi-merge -c T1_coassembly_anvio.db LEAP_META_01_T1_anvio_profile/PROFILE.db LEAP_META_02_T1_anvio_profile/PROFILE.db LEAP_META_03_T1_anvio_profile/PROFILE.db LEAP_META_04_T1_anvio_profile/PROFILE.db LEAP_META_05_T1_anvio_profile/PROFILE.db LEAP_META_06_T1_anvio_profile/PROFILE.db LEAP_META_07_T1_anvio_profile/PROFILE.db LEAP_META_08_T1_anvio_profile/PROFILE.db LEAP_META_09_T1_anvio_profile/PROFILE.db --enforce-hierarchical-clustering -o T1_merged_profile -T 64
```

```bash
anvi-import-collection metabat2_bins.txt -p T1_merged_nclu_profile/PROFILE.db -c T1_coassembly_anvio.db -C metabat2_bin_collection --contigs-mode
anvi-import-collection vamb_scaffold_bins.txt -p T1_merged_nclu_profile/PROFILE.db -c T1_coassembly_anvio.db -C vamb_bins --contigs-mode
```

```bash
anvi-rename-bins -c T1_coassembly_anvio.db -p T1_merged_nclu_profile/PROFILE.db --prefix mbat --collection-to-read metabat2_bin_collection --collection-to-write metabat2_50_bins --report-file rename_mbat50.txt --call-MAGs --min-completion-for-MAG 50 --max-redundancy-for-MAG 100 --exclude-bins

ssh -L localhost:8080:localhost:8080 ederrick@10.140.2.26
anvi-refine -p T1_merged_nclu_profile/PROFILE.db -c T1_coassembly_anvio.db -C metabat2_50_bins -b bin_XX
anvi-summarize -c T1_coassembly_anvio.db -p T1_merged_nclu_profile/PROFILE.db -o metabat_refined_bins -C metabat2_50_bins
```

realized this was missing some MAGs that checkM thought were fine

```bash
anvi-rename-bins -c T1_coassembly_anvio.db -p T1_merged_nclu_profile/PROFILE.db --prefix ex_mbat --collection-to-read metabat2_bin_collection --collection-to-write extra_metabat2_bins --report-file extra_mbat.txt --call-MAGs --min-completion-for-MAG 50 --max-redundancy-for-MAG 100
anvi-summarize -c T1_coassembly_anvio.db -p T1_merged_nclu_profile/PROFILE.db -o extra_metabat_refined_bins -C extra_metabat2_bins
```

in first set of refined MAGs
```bash
for f in *.fa; do mv $f ${f#mbat_}; done
for f in *.fa; do mv $f ${f%*-contigs.fa}.fa; done
```

in second set of refined MAGs. then move to one directory refined_MAGs
```bash
for f in *.fa; do mv $f MAG${f#ex_mbat_Bin}; done
for f in *.fa; do mv $f ${f%*-contigs.fa}.fa; done
```

check quality and dereplicate

```bash
dRep dereplicate refined_T1_MAGs -g refined_MAGs/*.fa -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 128
checkm lineage_wf refined_good_MAGs/ refined_mags_checkM/ -t 64 -x fa
```

map to refined MAGs

```bash
cd refined_good_MAGs
cat *.fa > T1_refined.fa
bowtie2-build T1_refined.fa T1_refined --threads 128

#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate bowtie2
parallel -j 6 --plus 'bowtie2 -x T1_refined -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/T1_refined.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

Annotate genes with prodigal

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate checkM
prodigal -i T1_refined.fa -d T1_refined_genes.fna -a T1_refined_genes.faa -o T1_refined_genes.gbk -p meta
```

create .stb file for inStrain using script from dRep

```bash
conda activate drep
parse_stb.py --reverse -f T1_refined_MAGs/*.fa -o T1_refined.stb
```

calls SNVs with inStrain

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 18 --plus 'inStrain profile {} T1_refined.fa -o {/P_T1_refined.bam/T1_refined_inStrain} -p 12 -g T1_refined_genes.fna -s T1_refined.stb  --min_read_ani 0.92 --min_mapq 1 --min_genome_coverage 1' ::: *T1_refined.bam
```

#### subsample to 5x (R script subsampling.R to get list of samtools commands to run (each mag x pond combination) and rerun instrain for some analyses. 

merge mags back together by timepoint

```bash
samtools merge LEAP_META_01_sub.bam *LEAP_META_01*
samtools merge LEAP_META_02_sub.bam *LEAP_META_02*
samtools merge LEAP_META_03_sub.bam *LEAP_META_03*
samtools merge LEAP_META_04_sub.bam *LEAP_META_04*
samtools merge LEAP_META_05_sub.bam *LEAP_META_05*
samtools merge LEAP_META_06_sub.bam *LEAP_META_06*
samtools merge LEAP_META_07_sub.bam *LEAP_META_07*
samtools merge LEAP_META_08_sub.bam *LEAP_META_08*
samtools merge LEAP_META_09_sub.bam *LEAP_META_09*
samtools merge LEAP_META_10_sub.bam *LEAP_META_10*
samtools merge LEAP_META_11_sub.bam *LEAP_META_11*
samtools merge LEAP_META_12_sub.bam *LEAP_META_12*
samtools merge LEAP_META_13_sub.bam *LEAP_META_13*
samtools merge LEAP_META_14_sub.bam *LEAP_META_14*
samtools merge LEAP_META_15_sub.bam *LEAP_META_15*
samtools merge LEAP_META_16_sub.bam *LEAP_META_16*
samtools merge LEAP_META_17_sub.bam *LEAP_META_17*
samtools merge LEAP_META_18_sub.bam *LEAP_META_18*
```

rerun instrain

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 18 --plus 'inStrain profile {} T1_refined.fa -o {/sub.bam/T1_subsamp_inStrain} -p 12 -g T1_refined_genes.fna -s T1_refined.stb  --min_read_ani 0.92 --min_mapq 1' ::: *sub.bam
```

#### annotate MAGs 

run eggnog

```bash
emapper.py -m diamond --itype CDS -i T1_refined_genes.fna -o eggnog_genes --output_dir /mfs/ederrick/chapter_1/09_anvio_binning/ --cpu 72
```

also annotate with bakta?

```bash
parallel -j 24 --plus 'bakta --db /mfs/ederrick/db {} -o {/.fa/_bakta} --threads 8' ::: *.fa
```

get EPSPS gene from annotation

```bash
#!/usr/bin/bash
for f in *_bakta
do
cd $f
sed -n '/3-phosphoshikimate 1-carboxyvinyltransferase/,/[^>]/p' ${f%*_bakta}.faa > ${f%*bakta}EPSPS.faa
cd ..
done
```

make fasta of all EPSPS then classify on webserver

```bash
for f in *.faa; do sed "s/>/>${f%*_EPSPS.faa}_/" $f > ${f%*.faa}_fix.faa; done
for f in *_fix.faa; do mv $f ${f%*_fix.faa}.faa; done
cat *.faa > MAG_EPSPS_genes.faa
```

annotate ARGs

```bash
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json /mfs/ederrick/card.json --local
parallel -j 12 --plus 'rgi main -i {} -n 8 --low_quality --clean -o {/.fa/_RGI}' ::: *.fa
parallel -j 12 --plus 'rgi main -i {} -n 8 --low_quality --clean --include_loose -o {/.fa/_loose_RGI}' ::: *.fa
```

#### calculate relative abundance of MAGs

with coverM

```bash
conda activate coverM
coverm genome --bam-files *.bam --genome-fasta-directory refined_good_MAGs --min-read-percent-identity 0.95 -m relative_abundance mean covered_bases count reads_per_base rpkm -o T1_refined_coverM.tsv --output-format sparse --min-covered-fraction 0 -tmx fa -t 64
```

with inStrain but profile all MAGs regardless of coverage

```bash
#!/usr/bin/bash
source /mfs/ederrick/.bash_profile
conda activate instrain

parallel -j 18 --plus 'inStrain profile {} T1_refined.fa -o {/P_T1_refined.bam/T1_refined_all_inStrain} -p 12 -g T1_refined_genes.fna -s T1_refined.stb  --min_read_ani 0.95 --min_mapq 1' ::: *T1_refined.bam
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


