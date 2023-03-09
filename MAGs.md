**Prepping files for InStrain**

I now have my dereplicated MAGs. I also created a file with scaffolds_to_genomes.sh that has two columns, the first is the scaffold/contig name and the second column is the MAG name. 

Create a fasta file of all MAGs. (In directory containing all non-redundant MAGs)

```bash
cat *.fa > ALL_MAGS.fa
```

Build a bowtie index of all MAGs. I did this with salloc.

```bash
module load bowtie2
bowtie2-build ALL_MAGS.fa
```

I will map the metagenomic reads from each timepoint to a database containing all the non-redudant MAGs.

```bash
#!/usr/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32

module load bowtie2

for f in *_R1.fastq.gz
do
f2="${f//R1.fastq.gz/R2.fastq.gz}"
out="${f//_R1.fastq.gz/default.sam}"
bowtie2 -x ALL_MAGS -1 $f -2 $f2 --threads 32 -S $out
done
```

I will call genes with prodigal for later use with instrain. I did this using salloc. Use -p meta to call genes in metagenomic mode. Otherwise the program will think all the genomes are the same genome.

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=1
module load prodigal

prodigal -i ALL_MAGS.fa -d mag_genes.fna -a mag_genes.faa -p meta
```

Need to covert sam to bam

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

module load samtools
for f in *default.sam
do
out="${f//default.sam/default.bam}"
samtools view -b $f > $out
done
```

Then need to sort bam files

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

module load samtools
for f in *default.bam
do
out="${f//default.bam/_sorted.bam}"
samtools sort $f > $out
done
```

Then need to index bam files (inStrain will also do this for you if you don't)

```bash
#!/usr/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G

module load samtools
for f in *sorted.bam
do
samtools index $f
done
```
**InStrain**
used inStrain version 1.6.4

"By default, if a read maps equally well to multiple genomes, Bowtie2 will pick one of the positions randomly and give the read a MAPQ score of 1. Thus, if you’d like to remove multi-mapped reads, you can set the minimum mapQ score to 2. 

There is also an option to set the minimum coverage for a genome to be profiled. I will set this to 5X which will save a lot of time.

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G

module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0
module load prodigal samtools/1.16

for f in *sorted.bam
do
out="${f//sorted.bam/cov5_instrain_profile}"
inStrain profile $f ALL_MAGS.fa -o $out -p 32 -g mag_genes.fna -s genome_scaffold.stb --min_mapq 1 --skip_mm_profiling --min_genome_coverage 5
done
```
This didn't leave me with many MAGs with over 5X coverage present in multiple ponds. 

I am going to merge the timepoints after the first pulse and the timepoints after the second pulse. 

There is no merging for I4 pulse 2 because there is no timepoint 5. There is no pulse 2 for L7.

In interactive node with 8 cpu and 4G memory.

```bash
module load samtools
samtools merge -o I4_pulse1.bam I4_2_sorted.bam I4_3_sorted.bam --threads 8
samtools merge -o I8_pulse1.bam I8_2_sorted.bam I8_3_sorted.bam --threads 8
samtools merge -o I8_pulse2.bam I8_4_sorted.bam I8_5_sorted.bam --threads 8
samtools merge -o K1_pulse1.bam K1_2_sorted.bam K1_3_sorted.bam --threads 8
samtools merge -o K1_pulse2.bam K1_4_sorted.bam K1_5_sorted.bam --threads 8
samtools merge -o L2_pulse1.bam L2_2_sorted.bam L2_3_sorted.bam --threads 8
samtools merge -o L2_pulse2.bam L2_4_sorted.bam L2_5_sorted.bam --threads 8
samtools merge -o L3_pulse1.bam L3_2_sorted.bam L3_3_sorted.bam --threads 8
samtools merge -o L3_pulse2.bam L3_4_sorted.bam L3_5_sorted.bam --threads 8
samtools merge -o L4_pulse1.bam L4_2_sorted.bam L4_3_sorted.bam --threads 8
samtools merge -o L4_pulse2.bam L4_4_sorted.bam L4_5_sorted.bam --threads 8
samtools merge -o L6_pulse1.bam L6_2_sorted.bam L6_3_sorted.bam --threads 8
samtools merge -o L6_pulse2.bam L6_4_sorted.bam L6_5_sorted.bam --threads 8
samtools merge -o L7_pulse1.bam L7_2_sorted.bam L7_3_sorted.bam --threads 8 
samtools merge -o L8_pulse1.bam L8_2_sorted.bam L8_3_sorted.bam --threads 8
samtools merge -o L8_pulse2.bam L8_4_sorted.bam L8_5_sorted.bam --threads 8
```

InStrain profile on merged files. Changing some parameters from earlier. The --database_mode flag does --min_read_ani 0.92 --skip_mm_profiling --min_genome_coverage 1 but I want --min_genome_coverage 5 so I will use those three parameters individually instead of using --database_mode

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G

module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0
module load prodigal samtools/1.16

for f in *pulse*.bam
do
out="${f//sorted.bam/cov5_instrain_profile}"
inStrain profile $f ALL_MAGS.fa -o $out -p 32 -g mag_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.92 --skip_mm_profiling --min_genome_coverage 5
done
```

InStrain released verson 1.7.1 so going to try that. Also changing min_read_ani to 0.99.

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G

source /home/ederrick/virtual_envs/instrain/bin/activate
module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0
module load prodigal samtools/1.16

for f in *pulse*.bam
do
out="${f//.bam/_ANI_99_instrain_profile}"
inStrain profile $f ALL_MAGS.fa -o $out -p 32 -g mag_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.99 --skip_mm_profiling --min_genome_coverage 5
done
```
