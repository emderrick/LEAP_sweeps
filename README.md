## LEAP 2021
This is a place for all analyses of the 2021 LEAP experiement. In this experiment we collected samples for metagenomic and whole-genome sequencing.

### Table of Contents
* RGI.md
  + running RGI with CARD database on metagenomic data
* Anvio_to_MAGs.md
  + starts with processing raw metagenomic sequencing reads
  + ends with dereplicating refined MAGs
  + output of one step is in renaming_bins_output.txt
* inStrain_MAGs.md
  + starts with building bowtie index of all MAGs for competetive mapping with Bowtie2
  + ends with running samtools depth to get coverage at all positions
  + uses scaffold_to_genome.sh
* mag_presence.R
  + makes big plot of all MAGs in each pond at each timepint
* graphing_presence.R
  + plots summary of what ponds my MAGs of interest are in
* combine_mag_files.R
  + combines SNV file, scaffold file, and genome file from instrain output together for downstream analysis
* filter_mag_SNVs_95.R
  + filteres SNVs based on position and coverage
  + filteres out timepoints I won't use downstream
* graphing_coverage.R
  + plots scaffold coverage for each mag of interest
* combine_depth.R
  + combines the samtools depth output with the inStrain output to get depth at all positions
* graphing_heatmaps.R
  + plots heatmaps of the reference frequency at each position there is an SNV
* find_snvs.R
  + finds positions where ref is above 0.5 in all ctl and above 0.5 in all gly
  + saves list of genes these positions are in
* SNVs_genes.md
  + starts with extracting genes of interest in list from all MAG genes
  + ends with ..?
* getEPSPS.sh
  + extracts EPSPS from prokka annotation 
  + need to add this part to another section
* wgs.md
  + processing of whole genome sequences
  + probably won't use this for anything





