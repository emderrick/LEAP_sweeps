### SNV and Gene Analysis of my 12 MAGs

#### Target Genes

After some processing of SNV and samtools files in R I have a list of positions and genes that may be targets of selection from R in a file called gene_locations.tsv 
I made a new file with just the first columnn

```bash
cut -f 1 gene_locations.tsv | sed 's/$/ /' > unique_genes.txt
```

Then I extracted the genes listed in that file from from the mag protein sequences previously called by prodigal. 

```bash
#!/usr/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=1
cat unique_genes.txt | xargs -i sed -n '/{}/,/--/p' mag_proteins.fa > gene_prot_seq.fa
sed -i 's/*//g' gene_prot_seq.fa 
sed -i '/--/d' gene_prot_seq.fa 
```
Then I clustered the protein sequecnes with cd-hit to see if any are related

```bash
#!/usr/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=
#SBATCH --cpus-per-task=1

module load cd-hit/4.8.1
cd-hit -i gene_prot_seq.fa -o clustered_mag_prot_0.7 -c 0.7 -sc 1 -n 5 -d 30
```
```bash
#!/usr/bin/bash
#SBATCH --time=00:05:00
#SBATCH --account=
#SBATCH --cpus-per-task=1

module load cd-hit/4.8.1
cd-hit -i gene_prot_seq.fa -o clustered_mag_prot_0.5 -c 0.5 -sc 1 -n 2 -d 30
```

I'm going to reannotate my candidate MAGs with bakta since I had previously done it with prokka. First I installed bakta with apptainer. This is bakta version 1.8.1.

```bash
module load apptainer/1.1.8
apptainer build bakta.sif docker://oschwengers/bakta:latest
```

Then I downloaded the database with

```bash
wget https://zenodo.org/record/7669534
tar -xzf db.tar.gz
rm db.tar.gz
```

Then I updated the database to include the amrfinder db with bakta's internal command. I couldn't do this to download the main db it didn't work.
```bash
apptainer shell -B /lustre07/scratch/ederrick bakta.sif
amrfinder_update --force_update --database db/amrfinderplus-db
```

Then I ran bakta on each MAG.

```bash
#!/usr/bin/bash
#SBATCH --time=00:30:00
#SBATCH --account=
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

module load apptainer/1.1.8
for file in candididate_mags/*.fa
do
out="${f//.fa/_bakta_output}"
apptainer run -B /lustre07/scratch/ederrick bakta.sif $file --db /home/ederrick/scratch/db --output $out --threads 4 
done
```
