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

**InStrain**

"By default, if a read maps equally well to multiple genomes, Bowtie2 will pick one of the positions randomly and give the read a MAPQ score of 1. Thus, if you’d like to remove multi-mapped reads, you can set the minimum mapQ score to 2. Must be >, not >=" - I will set MAPQ to 1 which means must be greater than 1 so 2.

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
