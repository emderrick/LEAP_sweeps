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
