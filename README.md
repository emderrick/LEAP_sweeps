# Workflow for the reference sequences

SPAdes
------
The first thing I did was with the whole-genome sequences of the original spike-in bacteria.
I assembled the raw fastq reads using spades 
I submitted this job to the cluster because it is RAM heavy.

	#!/usr/bin/bash
	for file in *R1.fastq*
	do
	filetwo="${file//R1.fastq/R2.fastq}"
	out="${file//R1.fastq/spadesoutput}"
	spades.py -o $out -1 $file -2 $filetwo --isolate
	done

After assembling with SPAdes you get output files of contigs.

PROKKA
------
I annotated these contigs using the program PROKKA.
I downloaded PROKKA through singluarity and ran the program with salloc asking for 10 CPUS
I used the --centre X --compliant flags to generate clean contig names for contigs_23.fa because the names were too long.

	#!/usr/bin/bash
	module load bioperl/1.7.7
	module load java/13.0.2
	module load hmmer/3.2.1
	module load singularity

	for file in *fa 
	do
	out="${file//.fa/prokkaoutput}"
	singularity exec prokka.sif prokka --outdir $out --cpus 10 $file
	done

GTDBtk
------
To classify the genomes, I used GTDBtk. This is also a RAM heavy program so I submitted it through the cluster.


	#!/usr/bin/bash
	
	module load python/3.7.7
	module load prodigal/2.6.3
	module load hmmer/3.2.1
	module load pplacer/1.1.alpha19
	module load fastani/1.32

	export GTDBTK_DATA_PATH=/home/ederrick/scratch/leap/release202

	gtdbtk classify_wf --genome_dir /home/ederrick/scratch/leap/spades_contigs --pplacer_cpus 1 --cpus 40 --extension fa --out_dir /home/ederrick/scratch/leap/gtdboutput

Resistance Gene Identifier
--------------------------

Next I downloaded the resistance gene identifier through a virtual environment to identify antibiotic resistance genes in the original WGS. This program uses the
CARD database from McMaster.
I did not have enough memory to run this on my computer so I asked for 10G through salloc.


	source /home/ederrick/scratch/virtual_envs/rgi/bin/activate
	
	module load nixpkgs/16.09
	module load gcc/7.3.0
	module load python/3.7.4
	module load scipy-stack
	module load bowtie2/2.3.5.1 
	module load samtools/1.10
	module load bamtools/2.4.1 
	module load diamond/0.9.32
	module load blast+

	rgi load --card_json /home/ederrick/scratch/virtual_environments/card.json

	rgi main -i contig_file -t contig -n 1 -o rgi_output -d wgs --include_loose
