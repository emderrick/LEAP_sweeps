I now have my dereplicated MAGs. I also created a file with scaffolds_to_genomes.sh that has two columns, the first is the scaffold/contig name and the second column is the MAG name. I will map the metagenomic reads from each timepoint to a database containing all the non-redudant MAGs.

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

I will call genes with prodigal for later use with instrain. I did this using salloc.

```bash
module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0
module load prodigal

prodigal -i ALL_MAGS.fa -d mag_genes.fna -a mag_genes.faa -o mag_genes.out
```
Need to covert sam to bam

need to sort bam files
