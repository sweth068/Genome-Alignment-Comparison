# Genome-Alignment-Comparison
Comparative Analysis of Alignment Tools on Diverse Genomic Data: Bowtie2, BWA, Minimap

This README file provides an overview of the scripts used for generating plots and tables from alignments performed using three different tools (Bowtie2, BWA, Minimap) on two reference genomes and a negative control.

Bash Scripts
Plot Generation Scripts
1. Bowtie2 Plot Script: bowtie_plot_script.sh
Generates plots from Bowtie2 alignment results.
2. BWA Plot Script: bwa_plot_script.sh
Generates plots from BWA alignment results.
3. Minimap Plot Script: mini_plot_script.sh
Generates plots from Minimap alignment results.
Note: The plots generated by these scripts include Mapping Quality vs. Coverage, Coverage, and Percent Identity.

Table Generation Scripts
1. Bowtie2 Table Script: bowtie_table_script.sh
Generates tables from Bowtie2 alignment results.
2. BWA Table Script: bwa_table_script.sh
Generates tables from BWA alignment results.
3. Minimap Table Script: mini_table_script.sh
Generates tables from Minimap alignment results.
Note: The tables generated by these scripts include Alignment Summary table, Coverage Statistics and Mapping Distribution Table.


Python Scripts
Plot Generation Scripts
1. Mapping Quality vs. Coverage Plot: cov_map.py
Generates a plot illustrating the relationship between mapping quality and coverage.
2. Coverage Plot: coverage_plot.py
Generates a coverage plot.
3. GC Content Plot: gc_content_plot.py
Generates a GC content plot. (This script needs to be run separately and is not included in the bash script.)
4. Percent Identity Plot: ident_perc.py
Generates a plot showing the percent identity.
Table Generation Scripts
1. Base Quality Distribution Table: base_quality.py
Generates a table illustrating the distribution of base quality. (This script needs to be run separately and is not included in the bash script.)
2. Mismatch Distribution Table: mismatch_distribution.py
Generates a table showing the distribution of mismatches.  (This script needs to be run separately and is not included in the bash script.)
3. Mismatch Count Table: mismatch_count.py  (This script needs to be run separately and is not included in the bash script.)
Generates a table indicating the count of mismatches.
Note: Ensure that Python and the required libraries are installed to run the Python scripts successfully.

These scripts provide a comprehensive analysis of alignment results, generating plots and tables for further evaluation of the data. Follow the instructions provided above to run the scripts and generate the desired visualizations and summaries.