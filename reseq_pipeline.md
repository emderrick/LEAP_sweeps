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
conda activate trimmomatic

for f in *R1.fastq.gz
do
trimmomatic PE -threads 64 -phred33 $f ${f%*R1.fastq.gz}R2.fastq.gz ${f%*R1.fastq.gz}QC_R1.fastq.gz ${f%*R1.fastq.gz}UP_R1.fastq.gz ${f%*R1.fastq.gz}QC_R2.fastq.gz ${f%*R1.fastq.gz}UP
_R2.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:12 MINLEN:36
done
```

inStrain needs paired reads. Trim reads to keep both pairs when they overlap and use these reads with inStrain

```bash
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
conda activate drep
dRep dereplicate checkM_T1_50_bins -g all_T1_bins/*.fa -comp 50 -con 10 --checkM_method lineage_wf --warn_aln 0.50 -p 64
```

```bash
seqkit stats -a *.fa > T1_MAGs_50_stats.txt
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
conda activate bowtie2

cd refined_good_MAGs
cat *.fa > T1_refined.fa
bowtie2-build T1_refined.fa T1_refined --threads 128

parallel -j 6 --plus 'bowtie2 -x T1_refined -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/T1_refined.bam} --write-index -@ 24' ::: *_P_R1.fastq.gz
```

Annotate genes with prodigal

```bash
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
conda activate instrain

parallel -j 18 --plus 'inStrain profile {} T1_refined.fa -o {/P_T1_refined.bam/T1_refined_inStrain} -p 12 -g T1_refined_genes.fna -s T1_refined.stb  --min_read_ani 0.92 --min_mapq 1 --min_genome_coverage 1' ::: *T1_refined.bam
```

#### subsample to 5x (R script subsampling.R to get list of samtools commands to run (each mag x pond combination) and rerun instrain. Do subsampling 5 times to see if we get the same results. sample 03 and 12 have been removed.

for each subsampling directory merge mags back together by timepoint

```bash
samtools merge LEAP_META_01_sub.bam *LEAP_META_01*
samtools merge LEAP_META_02_sub.bam *LEAP_META_02*
samtools merge LEAP_META_04_sub.bam *LEAP_META_04*
samtools merge LEAP_META_05_sub.bam *LEAP_META_05*
samtools merge LEAP_META_06_sub.bam *LEAP_META_06*
samtools merge LEAP_META_07_sub.bam *LEAP_META_07*
samtools merge LEAP_META_08_sub.bam *LEAP_META_08*
samtools merge LEAP_META_09_sub.bam *LEAP_META_09*
samtools merge LEAP_META_10_sub.bam *LEAP_META_10*
samtools merge LEAP_META_11_sub.bam *LEAP_META_11*
samtools merge LEAP_META_13_sub.bam *LEAP_META_13*
samtools merge LEAP_META_14_sub.bam *LEAP_META_14*
samtools merge LEAP_META_15_sub.bam *LEAP_META_15*
samtools merge LEAP_META_16_sub.bam *LEAP_META_16*
samtools merge LEAP_META_17_sub.bam *LEAP_META_17*
samtools merge LEAP_META_18_sub.bam *LEAP_META_18*
```

```bash
samtools merge LEAP_META_01_sub_02.bam *LEAP_META_01*
samtools merge LEAP_META_02_sub_02.bam *LEAP_META_02*
samtools merge LEAP_META_04_sub_02.bam *LEAP_META_04*
samtools merge LEAP_META_05_sub_02.bam *LEAP_META_05*
samtools merge LEAP_META_06_sub_02.bam *LEAP_META_06*
samtools merge LEAP_META_07_sub_02.bam *LEAP_META_07*
samtools merge LEAP_META_08_sub_02.bam *LEAP_META_08*
samtools merge LEAP_META_09_sub_02.bam *LEAP_META_09*
samtools merge LEAP_META_10_sub_02.bam *LEAP_META_10*
samtools merge LEAP_META_11_sub_02.bam *LEAP_META_11*
samtools merge LEAP_META_13_sub_02.bam *LEAP_META_13*
samtools merge LEAP_META_14_sub_02.bam *LEAP_META_14*
samtools merge LEAP_META_15_sub_02.bam *LEAP_META_15*
samtools merge LEAP_META_16_sub_02.bam *LEAP_META_16*
samtools merge LEAP_META_17_sub_02.bam *LEAP_META_17*
samtools merge LEAP_META_18_sub_02.bam *LEAP_META_18*
```

```bash
samtools merge LEAP_META_01_sub_03.bam *LEAP_META_01*
samtools merge LEAP_META_02_sub_03.bam *LEAP_META_02*
samtools merge LEAP_META_04_sub_03.bam *LEAP_META_04*
samtools merge LEAP_META_05_sub_03.bam *LEAP_META_05*
samtools merge LEAP_META_06_sub_03.bam *LEAP_META_06*
samtools merge LEAP_META_07_sub_03.bam *LEAP_META_07*
samtools merge LEAP_META_08_sub_03.bam *LEAP_META_08*
samtools merge LEAP_META_09_sub_03.bam *LEAP_META_09*
samtools merge LEAP_META_10_sub_03.bam *LEAP_META_10*
samtools merge LEAP_META_11_sub_03.bam *LEAP_META_11*
samtools merge LEAP_META_13_sub_03.bam *LEAP_META_13*
samtools merge LEAP_META_14_sub_03.bam *LEAP_META_14*
samtools merge LEAP_META_15_sub_03.bam *LEAP_META_15*
samtools merge LEAP_META_16_sub_03.bam *LEAP_META_16*
samtools merge LEAP_META_17_sub_03.bam *LEAP_META_17*
samtools merge LEAP_META_18_sub_03.bam *LEAP_META_18*
```

```bash
samtools merge LEAP_META_01_sub_04.bam *LEAP_META_01*
samtools merge LEAP_META_02_sub_04.bam *LEAP_META_02*
samtools merge LEAP_META_04_sub_04.bam *LEAP_META_04*
samtools merge LEAP_META_05_sub_04.bam *LEAP_META_05*
samtools merge LEAP_META_06_sub_04.bam *LEAP_META_06*
samtools merge LEAP_META_07_sub_04.bam *LEAP_META_07*
samtools merge LEAP_META_08_sub_04.bam *LEAP_META_08*
samtools merge LEAP_META_09_sub_04.bam *LEAP_META_09*
samtools merge LEAP_META_10_sub_04.bam *LEAP_META_10*
samtools merge LEAP_META_11_sub_04.bam *LEAP_META_11*
samtools merge LEAP_META_13_sub_04.bam *LEAP_META_13*
samtools merge LEAP_META_14_sub_04.bam *LEAP_META_14*
samtools merge LEAP_META_15_sub_04.bam *LEAP_META_15*
samtools merge LEAP_META_16_sub_04.bam *LEAP_META_16*
samtools merge LEAP_META_17_sub_04.bam *LEAP_META_17*
samtools merge LEAP_META_18_sub_04.bam *LEAP_META_18*
```

```bash
samtools merge LEAP_META_01_sub_05.bam *LEAP_META_01*
samtools merge LEAP_META_02_sub_05.bam *LEAP_META_02*
samtools merge LEAP_META_04_sub_05.bam *LEAP_META_04*
samtools merge LEAP_META_05_sub_05.bam *LEAP_META_05*
samtools merge LEAP_META_06_sub_05.bam *LEAP_META_06*
samtools merge LEAP_META_07_sub_05.bam *LEAP_META_07*
samtools merge LEAP_META_08_sub_05.bam *LEAP_META_08*
samtools merge LEAP_META_09_sub_05.bam *LEAP_META_09*
samtools merge LEAP_META_10_sub_05.bam *LEAP_META_10*
samtools merge LEAP_META_11_sub_05.bam *LEAP_META_11*
samtools merge LEAP_META_13_sub_05.bam *LEAP_META_13*
samtools merge LEAP_META_14_sub_05.bam *LEAP_META_14*
samtools merge LEAP_META_15_sub_05.bam *LEAP_META_15*
samtools merge LEAP_META_16_sub_05.bam *LEAP_META_16*
samtools merge LEAP_META_17_sub_05.bam *LEAP_META_17*
samtools merge LEAP_META_18_sub_05.bam *LEAP_META_18*
```

for each subsampling directory rerun instrain

```bash
conda activate instrain
parallel -j 18 --plus 'inStrain profile {} T1_refined.fa -o {/sub.bam/T1_subsamp_inStrain} -p 12 -g T1_refined_genes.fna -s T1_refined.stb  --min_read_ani 0.92 --min_mapq 1' ::: *sub.bam
```

```bash
parallel -j 4 --plus 'inStrain profile {} /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.fa -o {/.bam/_inStrain} -p 12 -g /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined_genes.fna -s /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.stb  --min_read_ani 0.92 --min_mapq 1 --min_genome_coverage 1' ::: *sub_02.bam
```

```bash
parallel -j 4 --plus 'inStrain profile {} /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.fa -o {/.bam/_inStrain} -p 12 -g /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined_genes.fna -s /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.stb  --min_read_ani 0.92 --min_mapq 1 --min_genome_coverage 1' ::: *sub_03.bam
```

```bash
parallel -j 4 --plus 'inStrain profile {} /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.fa -o {/.bam/_inStrain} -p 12 -g /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined_genes.fna -s /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.stb  --min_read_ani 0.92 --min_mapq 1 --min_genome_coverage 1' ::: *sub_04.bam
```

```bash
parallel -j 4 --plus 'inStrain profile {} /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.fa -o {/.bam/_inStrain} -p 12 -g /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined_genes.fna -s /mfs/ederrick/chapter_1/09_anvio_binning/T1_refined.stb  --min_read_ani 0.92 --min_mapq 1 --min_genome_coverage 1' ::: *sub_05.bam
```

#### annotate MAGs 

run eggnog

```bash
emapper.py -m diamond --itype CDS -i T1_refined_genes.fna -o eggnog_genes --output_dir /mfs/ederrick/chapter_1/09_anvio_binning/ --cpu 72
```

also annotate with bakta

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

### community composition of ponds

#### Subsample first.

find lowest read depth

```bash
seqkit stats -j 18 *.gz -a -T > paired_read_stats.tsv
```

subsample

```bash
parallel -j 18 --plus 'seqtk sample -s100 {} 87942340 | gzip > {/R1.fastq.gz/sub_R1.fastq.gz}' ::: *P_R1.fastq.gz
parallel -j 18 --plus 'seqtk sample -s100 {} 87942340 | gzip > {/R2.fastq.gz/sub_R2.fastq.gz}' ::: *P_R2.fastq.gz
```

#### run kraken

```bash
conda activate kraken2

for f in *P_sub_R1.fastq.gz
do
kraken2 --db /mfs/ederrick/k2_standard_db --threads 64 --output ${f%*R1.fastq.gz}kraken_output.txt --report ${f%*R1.fastq.gz}kraken_report.txt --paired $f ${f%*R1.fastq.gz}R2.fastq.gz
done
```

extract only bacterial reads

```bash
parallel -j 18 'python /mfs/ederrick/miniconda3/envs/kraken2/bin/extract_kraken_reads.py -k LEAP_META_{}_P_sub_kraken_output.txt -r LEAP_META_{}_P_sub_kraken_report.txt -s LEAP_META_{}_P_sub_R1.fastq.gz -s2 LEAP_META_{}_P_sub_R2.fastq.gz -o LEAP_META_{}_bac_R1.fastq.gz -o2 LEAP_META_{}_bac_R2.fastq.gz -t 2 --include-children' ::: {01..18}
```

```bash
for f in *bac_R1.fastq.gz
do
kraken2 --db /mfs/ederrick/k2_standard_db --threads 64 --output ${f%*R1.fastq.gz}kraken_output.txt --report ${f%*R1.fastq.gz}kraken_report.txt --paired $f ${f%*R1.fastq.gz}R2.fastq.gz
done
```

run bracken

```bash
conda activate bracken

for f in *bac_kraken_report.txt
do
bracken -d /mfs/ederrick/k2_standard_db -i $f -o ${f%*kraken_report.txt}species.bracken -w ${f%*kraken_report.txt}species_bracken_report.txt -r 100 -l S -t 10
bracken -d /mfs/ederrick/k2_standard_db -i $f -o ${f%*kraken_report.txt}genus.bracken -w ${f%*kraken_report.txt}genus_bracken_report.txt -r 100 -l G -t 10
bracken -d /mfs/ederrick/k2_standard_db -i $f -o ${f%*kraken_report.txt}class.bracken -w ${f%*kraken_report.txt}class_bracken_report.txt -r 100 -l C -t 10
bracken -d /mfs/ederrick/k2_standard_db -i $f -o ${f%*kraken_report.txt}phylum.bracken -w ${f%*kraken_report.txt}phylum_bracken_report.txt -r 100 -l P -t 10
done
```

calculate beta div - first separate files by day

```bash
conda activate kraken2

python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *species.bracken --type bracken > beta_div_species_day_0.txt
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *genus.bracken --type bracken > beta_div_genus_day_0.txt
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *class.bracken --type bracken > beta_div_class_day_0.txt
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *phylum.bracken --type bracken > beta_div_phylum_day_0.txt
```

```bash
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *species.bracken --type bracken > beta_div_species_day_28.txt
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *genus.bracken --type bracken > beta_div_genus_day_28.txt
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *class.bracken --type bracken > beta_div_class_day_28.txt
python /mfs/ederrick/miniconda3/envs/kraken2/bin/beta_diversity.py -i *phylum.bracken --type bracken > beta_div_phylum_day_28.txt
```

#### calculate relative abundance of MAGs

remap subsampled reads

```
conda activate bowtie2
parallel -j 6 --plus 'bowtie2 -x T1_refined -1 {} -2 {/R1.fastq.gz/R2.fastq.gz} --threads 24 | samtools sort -o {/R1.fastq.gz/sub_coverM.bam} --write-index -@ 24' ::: *sub_R1.fastq.gz
```

with coverM

```bash
conda activate coverM
coverm genome --bam-files *sub_coverM.bam --genome-fasta-directory refined_good_MAGs --min-read-percent-identity 0.95 -m relative_abundance mean covered_bases count reads_per_base rpkm -o T1_refined_sub_coverM.tsv --output-format sparse --min-covered-fraction 0 -x fa -t 64
```

