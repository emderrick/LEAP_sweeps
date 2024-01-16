**Downsampling MAGs to even out coverage**

First made a file of contig names for each MAG I need to extract from bam files
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

Extract the reads mapping to each MAG in each pond

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
Subsammple each mag/pond combination to the lowest coverage mag/pond combination

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
Then merge the bam files for each pond back together

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

Count mapped reads at each TP

```bash
#!/usr/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

module load StdEnv/2023
module load samtools/1.18

for magfile in *contigs.txt
do
  	for file in *.bam
        do
	mag="${magfile//_contigs.txt/}"
        total=$(samtools view -c $file)
        mapped=$(cat $magfile | tr "\n" " " | xargs samtools view -c -f 1 $file)
        echo $mag, $file, $mapped, $total
        done
done
```
