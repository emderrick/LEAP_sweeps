## LEAP 2021
This is a place for all analyses of the 2021 LEAP experiement. In this experiment we collected samples for metagenomic and whole-genome sequencing.

### what each file does for analysis for metagenomic data not in proper order
1. Anvio_to_MAGs.md
      + starts with processing raw metagenomic sequencing reads
      + ends with dereplicating refined MAGs
      + output of one step is in renaming_bins_output.txt
2. inStrain_MAGs.md
      + starts with building bowtie index of all MAGs for competetive mapping with Bowtie2
      + ends with running samtools depth to get coverage at all positions
      + uses scaffold_to_genome.sh
3. combine_filter_mag_files.R
      + combines SNV file, scaffold file, and genome file from instrain output together for downstream analysis
      + filters SNVs based on position from ends and coverage
      + filters out timepoints I won't use downstream
4. graphing_coverage.R
      + plots scaffold coverage for each mag of interest
5. combine_depth.R
      + combines the samtools depth output with the inStrain output to get depth at all positions
6. SNVs.R
      + finds positions where ref has avergae abs difference of at least 0.5 ctl vs gly
      + makes file for graphing heatmaps
7. sum_snvs.R
      + makes simple tables of all SNV and SNS sums for each MAG
8. grouped_heatmaps.R
      + makes the figures with the heatmaps, snvs, sns
9. get_genes.md
      + merges significant genes with eggnog output for COG enrichment
10. COG_enrichment_test.R
      +
11. COG_summary.R
      +
12. fisher_test.R
      +
13. gene_cov.R
      +
14. MAG_overview_plot.R
      +
15. NS_ratio.R
      +
16. snvs_to_matrix.R
      +

