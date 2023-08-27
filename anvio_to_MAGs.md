# Workflow for the metagenomic data

**Trimmomatic**

First trim the reads to remove the Illumina adaptors. I sumbitted this on Narval in a bash script.
I downloaded the adaptor file from the trimmomatic github page it is TruSeq3-PE.fa

```bash
#!/usr/bin/bash
for file in *R1_001.fastq.gz*
do
filetwo="${file//R1_001.fastq.gz/R2_001.fastq.gz}"
paired1="${file//R1_001.fastq.gz/paired_R1_001.fastq.gz}"
unpaired1="${file//R1_001.fastq.gz/unpaired_R1_001.fastq.gz}"
paired2="${filetwo//R2_001.fastq.gz/paired_R2_001.fastq.gz}"
unpaired2="${filetwo//R2_001.fastq.gz/unpaired_R2_001.fastq.gz}"
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 $file $
filetwo $paired1 $unpaired1 $paired2 $unpaired2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:
10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```
**MEGAHIT**

```bash
#!/usr/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=32
module load megahit/1.2.9

megahit -1 I4_1_R1.fastq.gz,I4_2_R1.fastq.gz,I4_3_R1.fastq.gz,I4_4_R1.fastq.gz,I4_5_R1.fastq.gz -2 I4_1_R2.fastq.gz,I4_2_R2.fastq.gz,I4_3_R2.fastq.gz,I4_4_R2.fastq.gz,I4_5_R2.fastq.
gz -t 32 -o L3_coassembly
```
once I have a coassembly for each pond I will follow the anvio pipeline.

**Anvio**

First reformat fasta deflines and remove short contigs (<1000bp)

```bash
#!/usr/bin/bash
module load scipy-stack/2021a
module load python/3.7
anvi-script-reformat-fasta I4_contigs.fa -o I4_contigs_fixed.fa -l 1000 --simplify-names
```
Then rename it to get rid of fixed in title
Then make a contigs database of each. This uses prodigal for gene calling.

**bowtie2**

Once we have our coassemblies we need to map the reads at each timepoint to the coassembly.
We will use bowtie2.

```bash
#!/usr/bin/bash
module load bowtie2

for f in *_R1.fastq.gz
do
f2="${f//R1.fastq.gz/R2.fastq.gz}"
out="${f//_R1.fastq.gz/.sam}"
bowtie2 -x I4_contigs -1 $f -2 $f2 --threads 40 -S $out
done
```
**samtools**

After running bowtie2 we will have a .sam file for each sample. We need to convert the sam file to a bam file for anvio. Program needed lots of memory to work properly (20G).

```bash
 #!/usr/bin/bash
for f in *I4.sam*
do
out="${f//.sam/raw.bam}"
samtools view -b $f > $out
done
```

Then index the bam file. This program uses samtools and needs lots of memory (10G)
```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a

for f in *raw.bam
do
out="${f//raw.bam/.bam}"
anvi-init-bam $f -o $out
done
deactivate
```


```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.fa
do
out="${f//.fa/.db}"
anvi-gen-contigs-database -f $f -o $out -T 16
done
deactivate

```

Then run HMMS

```bash

#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.db
do
anvi-run-hmms -c $f -T 16
done
deactivate
```
The run scg-taxonomy

```bash
#!/usr/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=1

source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.db
do
anvi-run-scg-taxonomy -c $f
done
deactivate

```

Then we profile each bam file with the corresponding coassembly

```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7

for f in *_I4.bam*
do
anvi-profile -i $f -c I4_contigs.db -T 16
done
deactivate
```

Then we use cluster contigs to make a collection.

```bash

source anvio/bin/activate
module load scipy-stack/2021a
anvi-cluster-contigs -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection --driver concoct --just-do-it

```
Visualizing bins using anvi-interactive

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-interactive -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection --server-only -P 8081

```
Summarize bins

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-summarize -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection -o I4_bins_summary

```

Then refine each bin of each pond with anvi-refine

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-refine -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection --server-only -P 8081 -b Bin_01

```

Then rename all bins to MAGs that are <10 red and >70 comp and include pond name in the MAG

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-rename-bins -c I4_contigs.db -p I4_profiles_merged/PROFILE.db --prefix I4 --collection-to-read I4_collection --collection-to-write I4_MAGS_FINAL --report-file I4_rename.txt --call-MAGs --min-completion-for-MAG 70 --max-redundancy-for-MAG 10 

```

Then summarize final collection of refined bins

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-summarize -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_MAGS_FINAL -o I4_BINS_FINAL_SUMMARY

```
Then extract a fasta file of each MAG from each pond and put in new directory called redundant_MAGS. Use salloc and 10G of memory.

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

#!/usr/bin/bash
for f in I4_BINS_FINAL_SUMMARY/bin_by_bin/*MAG*
 do
	name="${f//I4_BINS_FINAL_SUMMARY}"
	name1="${name///bin_by_bin}"
	name2="${name1///}"
	anvi-script-reformat-fasta $f/$name2-contigs.fa --simplify-names --prefix $name2 -o redundant_MAGS/$name2.fa
done

```
Then I created a tab deliminated text file "all_redundant_mags.txt" that contained the name of each MAG and its file path.

Then I dereplicated the MAGs with an anvio program using fastANI. I used 16 threads. I used a 98% ANI.

```bash

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer
module load fastani/1.32

anvi-dereplicate-genomes -f all_redundant_mags.txt -o dereplicate_MAGS --similarity-threshold 0.98 --program fastANI -T 16
```
