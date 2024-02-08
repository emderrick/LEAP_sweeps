# Complete workflow for raw reads to instrain output

**Trimmomatic**

First trim the reads to remove the Illumina adaptors.
I downloaded the adaptor file from the trimmomatic github page it is TruSeq3-PE.fa

```bash
#!/usr/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=16
module load trimmomatic

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
Create a co-assembly for each pond. Example below if for pond I4. Repeated for each of the nine ponds.

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

First reformat fasta deflines and remove short contigs (<1000bp). Done in CC interactive node.

```bash
#!/usr/bin/bash
module load scipy-stack/2021a
module load python/3.7
anvi-script-reformat-fasta I4_contigs.fa -o I4_contigs_fixed.fa -l 1000 --simplify-names

```
Then rename it to get rid of fixed in title

```bash
for f in *fixed.txt
do
out="${f//_fixed.txt/.txt}"
mv $f $out
done

```

**bowtie2**

Once we have our coassemblies we need to map the reads at each timepoint to the coassembly. Need to do this for each pond. 5 timepoints from pond I4 map to I4 etc.
We will use bowtie2.

```bash
#!/usr/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=40

module load bowtie2/2.4.4

for f in *_R1.fastq.gz
do
f2="${f//R1.fastq.gz/R2.fastq.gz}"
out="${f//R1.fastq.gz/I4.sam}"
bowtie2 -x I4_contigs -1 $f -2 $f2 --threads 40 -S $out
done

```

**samtools**

After running bowtie2 we will have a .sam file for each sample. We need to convert the sam file to a bam file for anvio.

```bash
#!/usr/bin/bash
#SBATCH --account=
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G

module load samtools
for f in *.sam
do
out="${f//.sam/raw.bam}"
samtools view -b $f > $out
done

```

Then index the bam file. did in CC interactive node.

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

Then generate contig database for each pond.

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=36

source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7
module load prodigal
module load diamond
module load hmmer

for f in *.fa
do
out="${f//.fa/.db}"
anvi-gen-contigs-database -f $f -o $out -T 36
done
deactivate

```

Then run HMMS

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G

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

Then run scg-taxonomy

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
#SBATCH --time=01:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G

source anvio/bin/activate
module load scipy-stack/2021a
module load python/3.7

for f in *I4.bam
do
out="${f//_.bam/_profile}"
anvi-profile -i $f -c I4_contigs.db -T 16 --output-dir $out --min-contig-length 2500
done
deactivate

```

Then we merge profiles from each pond to make a merged profile for each pond.

```bash
#!/usr/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G

source anvio/bin/activate
module load scipy-stack/2021a
anvi-merge I4_1_I4_profile/PROFILE.db I4_2_I4_profile/PROFILE.db I4_3_I4_profile/PROFILE.db I4_4_I4_profile/PROFILE.db -o I4_profiles_merged -c I4_contigs.db --sample-name I4
deactivate
```

Then we use cluster contigs to make a collection for each pond. did in CC interactive node.

```bash
source anvio/bin/activate
module load scipy-stack/2021a
anvi-cluster-contigs -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection --driver concoct --just-do-it
deactivate

```

Then we can visualize bins for each pond using anvi-interactive.

```bash
source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-interactive -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection --server-only -P 8081
deactivate

```

Summarize bins to get a list of bins and their completion and redundancy so I can choose which to manually refine. Do for each pond. did in CC interactive node.

```bash
source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-summarize -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection -o I4_bins_summary
deactivate

```

Then refine each bin of each pond with anvi-refine. This is an example command to refine one bin. I did this for a thousand bins.

```bash
source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-refine -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_collection --server-only -P 8081 -b Bin_01
deactivate

```

Then rename all bins and name it a MAG if bin is <10 redundant and >70 complete and include pond name in the MAG. Do for each pond. did in CC interactive node.

```bash
source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-rename-bins -c I4_contigs.db -p I4_profiles_merged/PROFILE.db --prefix I4 --collection-to-read I4_collection --collection-to-write I4_MAGS_FINAL --report-file I4_rename.txt --call-MAGs --min-completion-for-MAG 70 --max-redundancy-for-MAG 10 
deactivate

```

Then summarize final collection of refined bins. Do for each pond. did in CC interactive node.

```bash
source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

anvi-summarize -p I4_profiles_merged/PROFILE.db -c I4_contigs.db -C I4_MAGS_FINAL -o I4_BINS_FINAL_SUMMARY
deactivate

```

Then extract a fasta file of each MAG from each pond and put in new directory called redundant_MAGS. Do for each pond. did in CC interactive node.

```bash
#!/usr/bin/bash
source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer

for f in I4_BINS_FINAL_SUMMARY/bin_by_bin/*MAG*
do
name="${f//I4_BINS_FINAL_SUMMARY}"
name1="${name///bin_by_bin}"
anvi-script-reformat-fasta $f/$name1-contigs.fa --simplify-names --prefix $name1 -o redundant_MAGS/$name1.fa
done
deactivate

```
Then I created a tab deliminated text file "all_redundant_mags.txt" that contained the name of each MAG and its file path.

Then I dereplicated the MAGs through anvio program using fastANI at 98% ANI.

```bash
#!/usr/bin/bash
#SBATCH --time=00:20:00
#SBATCH --account=
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

source anvio/bin/activate
module load scipy-stack/2021a
module load diamond
module load prodigal
module load hmmer
module load fastani/1.32

anvi-dereplicate-genomes -f all_redundant_mags.txt -o dereplicate_MAGS --similarity-threshold 0.98 --program fastANI -T 8
deactivate

```
The non-redundant genomes output with strange names so need to fix them.

```bash
#!/usr/bin/bash
for f in *.fa
do
out="${f:0:15}"
mv "$f" $out
done

```

**Prepping files for InStrain**

Now create a file for inStrain that has a column with the scaffold/contig name and the second with the MAG name. 

```bash
#!/usr/bin/bash
for f in *.fa
do
out="${f//.fa/_scaffolds.txt}"
genomes="${out//scaffolds.txt/genomes.txt}"
genscaf="${genomes//genomes.txt/genomesscaffolds.txt}"
grep MAG $f | sed 's/>//' > $out
cut -d_ -f1,2,3 $out > $genomes 
paste -d"\t" $out $genomes > $genscaf
done
cat *genomesscaffolds.txt > genome_scaffold.txt
rm *scaffolds.txt
rm *genomes.txt

```

Then create a fasta file of all non-redundant MAGs

```bash
cat *.fa > ALL_MAGS.fa

```

Then I called genes with prodigal for later use with inStrain. I used -p meta to call genes in metagenomic mode (in future versions this is -anon but CC version is old). Otherwise the program will think all the genomes are the same genome.

```bash
#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=1
module load prodigal

prodigal -i ALL_MAGS.fa -d mag_genes.fna -a mag_genes.faa -p meta

```

Then build a bowtie index of all MAGs. 

```bash
#!/usr/bin/bash
#SBATCH --time=00:30:00
#SBATCH --account=
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G

module load bowtie2/2.4.4
bowtie2-build ALL_MAGS.fa ALL_MAGS

```

Then I mapped the metagenomic reads from each timepoint to the bowtie index of the database containing all the non-redundant MAGs.

```bash
#!/usr/bin/bash
#SBATCH --time=8:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=32

module load bowtie2

for f in *_R1.fastq.gz
do
f2="${f//R1.fastq.gz/R2.fastq.gz}"
out="${f//_R1.fastq.gz/.sam}"
bowtie2 -x ALL_MAGS -1 $f -2 $f2 --threads 32 -S $out
done

```

Then convert sam to bam

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

module load samtools
for f in *.sam
do
out="${f//.sam/.bam}"
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
for f in *.bam
do
out="${f//.bam/_sorted.bam}"
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

Coverage for each MAG is low so I merged the timepoints after the first pulse and the timepoints after the second pulse together to increase the coverage. 
There is no merging for I4 pulse 2 because there is no timepoint 5 (never got sequenced?). There is no pulse 2 for L7 (pond leaked). Did this in CC interactive node.

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

**InStrain**

"By default, if a read maps equally well to multiple genomes, Bowtie2 will pick one of the positions randomly and give the read a MAPQ score of 1. Thus, if you’d like to remove multi-mapped reads, you can set the minimum mapQ score to 2. 
There is also an option to set the minimum coverage for a genome to be profiled. I will set this to 5X which will save a lot of time to run inStrain.

Then I ran inStrain profile on merged files. The --database_mode flag does --min_read_ani 0.92 --skip_mm_profiling --min_genome_coverage 1 but I want --min_genome_coverage 5 so I will use those three parameters individually instead of using --database_mode

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G

source /home/ederrick/virtual_envs/instrain/bin/activate
module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0
module load prodigal samtools

for f in *pulse*.bam
do
out="${f//sorted.bam/cov5_instrain_profile}"
inStrain profile $f ALL_MAGS.fa -o $out -p 32 -g mag_genes.fna -s genome_scaffold.stb--min_mapq 2 --min_read_ani 0.95 --skip_mm_profiling --min_genome_coverage 5 --min_freq 0
done
```

**Annotation of MAG genes**

Annotation of my candidate MAGs with bakta. I installed bakta in a virtual environment. This is bakta version 1.8.1.

```bash
module load python/3.10
virtualenv --no-download bakta
source bakta/bin/activate
pip install --no-index --upgrade pip
pip install bakta

wget www.drive5.com/pilercr/pilercr1.06.tar.gz # install pilercr
tar -zxvf pilercr1.06.tar.gz
cd pilercr1.06
make
cd ..

git clone https://github.com/ncbi/amr.git
cd amr
git checkout master
make
make clean
make DEFAULT_DB_DIR=/home/ederrick/scratch/dbs/bakta/amrfinderplus-db/
make install INSTALL_DIR=/home/ederrick/venvs/

```

Then I downloaded the database with

```bash
wget https://zenodo.org/record/7669534
tar -xzf db.tar.gz
rm db.tar.gz

```

Then I updated the database to include the amrfinder db with bakta's internal command.

```bash
export PATH=/home/ederrick/venvs/amr:/home/ederrick/venvs/pilercr1.06:$PATH
amrfinder_update --force_update --database db/amrfinderplus-db

```

Then I ran bakta on each MAG.

```bash
#!/usr/bin/bash
#SBATCH --account=
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=16

source /home/ederrick/venvs/bakta/bin/activate

module load StdEnv/2020  gcc/9.3.0 trnascan-se/2.0.12
module load aragorn/1.2.41
module load infernal/1.1.4
module load hmmer/3.3.2
module load diamond/2.0.15
module load blast+/2.13.0
module load circos/0.69-9

export PATH=/home/ederrick/venvs/amr:/home/ederrick/venvs/pilercr1.06:$PATH

for file in candidate_MAGs/*.fa
do
out="${file//.fa/_full_bakta_output}"
bakta $file --db /home/ederrick/scratch/dbs/bakta  --output $out --threads 16
done
deactivate

```

Also annotating with EggNOG to get COG categories for each gene because bakta doesn't output it nicely. Don't end up using bakta annotation.

```bash
#!/usr/bin/bash
#SBATCH --account=
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=36

source /home/ederrick/virtual_envs/eggnog/bin/activate
export EGGNOG_DATA_DIR=/home/ederrick/virtual_envs/eggnog/

emapper.py  -m diamond --itype CDS -i oct_mag_genes.fna -o eggnog_genes_oct_nt --output_dir /home/ederrick/scratch/eggnog_output_oct --cpu 36

deactivate

```

**Downsampling MAGs to even out coverage**

First I made a file of contig names for each MAG I need to extract from bam files

```bash
cat all_mags_copy.fa | grep "I4_MAG_00006" > I4_MAG_00006_names.txt
cat all_mags_copy.fa | grep "I4_MAG_00065" > I4_MAG_00065_names.txt
cat all_mags_copy.fa | grep "L2_MAG_00052" > L2_MAG_00052_names.txt
cat all_mags_copy.fa | grep "L3_MAG_00058" > L3_MAG_00058_names.txt
cat all_mags_copy.fa | grep "L4_MAG_00099" > L4_MAG_00099_names.txt
cat all_mags_copy.fa | grep "L7_MAG_00020" > L7_MAG_00020_names.txt
cat all_mags_copy.fa | grep "L7_MAG_00028" > L7_MAG_00028_names.txt
cat all_mags_copy.fa | grep "L7_MAG_00043" > L7_MAG_00043_names.txt
cat all_mags_copy.fa | grep "L8_MAG_00011" > L8_MAG_00011_names.txt
cat all_mags_copy.fa | grep "L8_MAG_00019" > L8_MAG_00019_names.txt
cat all_mags_copy.fa | grep "L8_MAG_00042" > L8_MAG_00042_names.txt
for file in *names.txt; do cat $file | tr -d '>' > "${file%_names.txt}"_contigs.txt; done

```

Then extract the reads mapping to each MAG in each pond

```bash
#!/usr/bin/bash
#SBATCH --time=00:15:00
#SBATCH --account=
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

module load StdEnv/2023
module load samtools/1.18

for magfile in *contigs.txt
do
	for file in *pulse1.bam
	do
	mag="${magfile//_contigs.txt/}"
  out="${file//.bam/_$mag.bam}"
	cat $magfile | tr "\n" " " | xargs samtools view -bh $file  > $out
	done
done

```

Then subsammple each mag/pond combination to the lowest coverage mag/pond combination. did in CC interactive node.

```bash
samtools view -bh --subsample 0.228 I4_pulse1_I4_MAG_00006.bam > I4_pulse1_I4_MAG_00006_subsamp.bam
samtools view -bh --subsample 1.000 I8_pulse1_I4_MAG_00006.bam > I8_pulse1_I4_MAG_00006_subsamp.bam
samtools view -bh --subsample 0.676 L2_pulse1_I4_MAG_00006.bam > L2_pulse1_I4_MAG_00006_subsamp.bam
samtools view -bh --subsample 0.896 L8_pulse1_I4_MAG_00006.bam > L8_pulse1_I4_MAG_00006_subsamp.bam

samtools view -bh --subsample 1.000 K1_pulse1_I4_MAG_00065.bam > K1_pulse1_I4_MAG_00065_subsamp.bam 
samtools view -bh --subsample 0.491 I4_pulse1_I4_MAG_00065.bam > I4_pulse1_I4_MAG_00065_subsamp.bam
samtools view -bh --subsample 0.206 L3_pulse1_I4_MAG_00065.bam > L3_pulse1_I4_MAG_00065_subsamp.bam  
samtools view -bh --subsample 0.464 L4_pulse1_I4_MAG_00065.bam > L4_pulse1_I4_MAG_00065_subsamp.bam
samtools view -bh --subsample 0.497 I8_pulse1_I4_MAG_00065.bam > I8_pulse1_I4_MAG_00065_subsamp.bam  
samtools view -bh --subsample 0.833 L6_pulse1_I4_MAG_00065.bam > L6_pulse1_I4_MAG_00065_subsamp.bam 

samtools view -bh --subsample 1.000 K1_pulse1_L2_MAG_00052.bam > K1_pulse1_L2_MAG_00052_subsamp.bam
samtools view -bh --subsample 0.716 I4_pulse1_L2_MAG_00052.bam > I4_pulse1_L2_MAG_00052_subsamp.bam
samtools view -bh --subsample 0.472 L4_pulse1_L2_MAG_00052.bam > L4_pulse1_L2_MAG_00052_subsamp.bam
samtools view -bh --subsample 0.208 I8_pulse1_L2_MAG_00052.bam > I8_pulse1_L2_MAG_00052_subsamp.bam  
samtools view -bh --subsample 0.197 L2_pulse1_L2_MAG_00052.bam > L2_pulse1_L2_MAG_00052_subsamp.bam

samtools view -bh --subsample 1.000 L3_pulse1_L3_MAG_00058.bam > L3_pulse1_L3_MAG_00058_subsamp.bam 
samtools view -bh --subsample 0.250 L4_pulse1_L3_MAG_00058.bam > L4_pulse1_L3_MAG_00058_subsamp.bam 
samtools view -bh --subsample 0.571 L7_pulse1_L3_MAG_00058.bam > L7_pulse1_L3_MAG_00058_subsamp.bam 
samtools view -bh --subsample 0.362 L8_pulse1_L3_MAG_00058.bam > L8_pulse1_L3_MAG_00058_subsamp.bam

samtools view -bh --subsample 0.771 L4_pulse1_L4_MAG_00099.bam > L4_pulse1_L4_MAG_00099_subsamp.bam
samtools view -bh --subsample 1.000 L2_pulse1_L4_MAG_00099.bam > L2_pulse1_L4_MAG_00099_subsamp.bam 
samtools view -bh --subsample 0.227 L6_pulse1_L4_MAG_00099.bam > L6_pulse1_L4_MAG_00099_subsamp.bam 
samtools view -bh --subsample 0.152 L7_pulse1_L4_MAG_00099.bam > L7_pulse1_L4_MAG_00099_subsamp.bam  
samtools view -bh --subsample 0.172 L8_pulse1_L4_MAG_00099.bam > L8_pulse1_L4_MAG_00099_subsamp.bam

samtools view -bh --subsample 1.000 K1_pulse1_L7_MAG_00020.bam > K1_pulse1_L7_MAG_00020_subsamp.bam
samtools view -bh --subsample 0.169 L3_pulse1_L7_MAG_00020.bam > L3_pulse1_L7_MAG_00020_subsamp.bam 
samtools view -bh --subsample 0.044 L4_pulse1_L7_MAG_00020.bam > L4_pulse1_L7_MAG_00020_subsamp.bam 
samtools view -bh --subsample 0.612 L2_pulse1_L7_MAG_00020.bam > L2_pulse1_L7_MAG_00020_subsamp.bam

samtools view -bh --subsample 0.407 I8_pulse1_L7_MAG_00028.bam > I8_pulse1_L7_MAG_00028_subsamp.bam
samtools view -bh --subsample 0.934 L2_pulse1_L7_MAG_00028.bam > L2_pulse1_L7_MAG_00028_subsamp.bam  
samtools view -bh --subsample 1.000 L6_pulse1_L7_MAG_00028.bam > L6_pulse1_L7_MAG_00028_subsamp.bam
samtools view -bh --subsample 0.409 L7_pulse1_L7_MAG_00028.bam > L7_pulse1_L7_MAG_00028_subsamp.bam 

samtools view -bh --subsample 0.182 L4_pulse1_L7_MAG_00043.bam > L4_pulse1_L7_MAG_00043_subsamp.bam  
samtools view -bh --subsample 0.331 L6_pulse1_L7_MAG_00043.bam > L6_pulse1_L7_MAG_00043_subsamp.bam
samtools view -bh --subsample 1.000 L7_pulse1_L7_MAG_00043.bam > L7_pulse1_L7_MAG_00043_subsamp.bam

samtools view -bh --subsample 0.031 I8_pulse1_L8_MAG_00011.bam > I8_pulse1_L8_MAG_00011_subsamp.bam
samtools view -bh --subsample 1.000 L2_pulse1_L8_MAG_00011.bam > L2_pulse1_L8_MAG_00011_subsamp.bam  
samtools view -bh --subsample 0.100 L8_pulse1_L8_MAG_00011.bam > L8_pulse1_L8_MAG_00011_subsamp.bam

samtools view -bh --subsample 0.448 I8_pulse1_L8_MAG_00019.bam > I8_pulse1_L8_MAG_00019_subsamp.bam 
samtools view -bh --subsample 0.408 L2_pulse1_L8_MAG_00019.bam > L2_pulse1_L8_MAG_00019_subsamp.bam  
samtools view -bh --subsample 0.934 L6_pulse1_L8_MAG_00019.bam > L6_pulse1_L8_MAG_00019_subsamp.bam  
samtools view -bh --subsample 1.000 L8_pulse1_L8_MAG_00019.bam > L8_pulse1_L8_MAG_00019_subsamp.bam

samtools view -bh --subsample 0.760 K1_pulse1_L8_MAG_00042.bam > K1_pulse1_L8_MAG_00042_subsamp.bam  
samtools view -bh --subsample 0.346 L3_pulse1_L8_MAG_00042.bam > L3_pulse1_L8_MAG_00042_subsamp.bam
samtools view -bh --subsample 0.159 L4_pulse1_L8_MAG_00042.bam > L4_pulse1_L8_MAG_00042_subsamp.bam  
samtools view -bh --subsample 0.686 I8_pulse1_L8_MAG_00042.bam > I8_pulse1_L8_MAG_00042_subsamp.bam
samtools view -bh --subsample 1.000 L8_pulse1_L8_MAG_00042.bam > L8_pulse1_L8_MAG_00042_subsamp.bam
```
Then merge the bam files for each pond back together. Did in CC interactive node.

```bash
samtools merge subsamp_K1_pulse1.bam *K1_pulse1*
samtools merge subsamp_I4_pulse1.bam *I4_pulse1*
samtools merge subsamp_L3_pulse1.bam *L3_pulse1*
samtools merge subsamp_L4_pulse1.bam *L4_pulse1*
samtools merge subsamp_I8_pulse1.bam *I8_pulse1*
samtools merge subsamp_L2_pulse1.bam *L2_pulse1*
samtools merge subsamp_L6_pulse1.bam *L6_pulse1*
samtools merge subsamp_L7_pulse1.bam *L7_pulse1*
samtools merge subsamp_L8_pulse1.bam *L8_pulse1*

```

Then rerun instrain with downsampled bam files.

```bash
#!/usr/bin/bash
#SBATCH --time=4:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G

source /home/ederrick/virtual_envs/instrain/bin/activate
module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0
module load prodigal samtools

for f in *pulse*.bam
do
out="${f//sorted.bam/cov5_instrain_profile}"
inStrain profile $f ALL_MAGS.fa -o $out -p 32 -g oct_mag_genes.fna -s genome_scaffold.stb --min_mapq 2 --min_read_ani 0.95
done

```
