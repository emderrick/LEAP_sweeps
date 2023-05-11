### SNV and Gene Analysis of my 12 MAGs

#### Target Genes

After some processing of SNV and samtools files in R I have a list of positions and genes that may be targets of selection. I imported a list of just those genes into narval.

I extracted only the genes I'm interested in. 4th command took awhile so I submitted it in a compute node.
```bash
cat only_genes.csv | tr -d '"' > only_genes.txt
uniq -i only_genes.txt > unique_genes.txt
sed  '/>/i--' mag_genes.fna > mag_genes.fasta
cat unique_genes.txt | xargs -i sed -n '/{}/,/--/p' mag_genes.fasta > genes.fa
sed -i '/--/d' genes.fa 
```
I'm going to use eggNOG to annotate my genes because prokka can't use just genes it needs to be contigs
I installed eggNOG in a virtual environment
```bash
python3 -m venv eggnog
source eggnog/bin/activate
python3 -m pip install eggnog-mapper
```
then I'll install the required databases in the eggnog folder
```bash
download_eggnog_data.py --data_dir /home/ederrick/eggnog/
```
to run eggNOG on my gene that are nucleotide sequences
```bash
salloc --time=01:00:00 --account= --cpus-per-task=4 --mem-per-cpu=6G
source /home/ederrick/eggnog/bin/activate
export EGGNOG_DATA_DIR=/home/ederrick/eggnog
mkdir eggNOG_genes_output
emapper.py -m diamond --itype CDS -i genes.fa --output_dir eggNOG_genes_output --cpu 4 -o eggNOG_genes
deactivate
```
