## LEAP 2021
This is a place for all analyses of the 2021 LEAP experiement. In this experiment we collected samples for metagenomic and whole-genome sequencing.

### what each file does for analysis for metagenomic data
Anvio_to_MAGs.md
      + starts with processing raw metagenomic sequencing reads
      + ends with dereplicating refined MAGs
      + output of one step is in renaming_bins_output.txt
inStrain_MAGs.md
      + starts with building bowtie index of all MAGs for competetive mapping with Bowtie2
      + ends with running samtools depth to get coverage at all positions
      + uses scaffold_to_genome.sh

combine_filter_mag_files.R
      + combines SNV file, scaffold file, and genome file from instrain output together for downstream analysis
      + filters SNVs based on position from ends and coverage
      + filters out timepoints I won't use downstream
graphing_coverage.R
      + plots scaffold coverage for each mag of interest
combine_depth.R
      + combines the samtools depth output with the inStrain output to get depth at all positions
SNVs.R
      + finds positions where ref has avergae abs difference of at least 0.5 ctl vs gly
      + makes file for graphing heatmaps
sum_snvs.R
      + makes simple tables of all SNV and SNS sums for each MAG
grouped_heatmaps.R
      + makes the figures with the heatmaps, snvs, sns
get_genes.md
      + merges significant genes with eggnog output for COG enrichment
COG_enrichment_test.R
      +
COG_summary.R
      +
fisher_test.R
      +
gene_cov.R
      +
MAG_overview_plot.R
      +
NS_ratio.R
      +
snvs_to_matrix.R
      +

