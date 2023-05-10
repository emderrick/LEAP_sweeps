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
