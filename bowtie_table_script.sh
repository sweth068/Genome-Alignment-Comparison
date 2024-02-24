#!/bin/bash

#Directory containing sorted BAM files
bam_dir="/home/swetha/final/bowtie_plot_output"

#Directory containing Python files
python_files_dir="/home/swetha/final"

#Output directory for tables
output_table_dir="/home/swetha/final/bowtie_table_output"

#Create the output directory if it doesn't exist
mkdir -p "$output_table_dir"

#Iterate over each sorted BAM file in the directory
for bam_file in "$bam_dir"/sorted_*.bam; do
    #Extract the base name without extension
    base_name=$(basename "$bam_file" .bam)

    #Alignment Summary Table
    echo "Alignment Summary Table for $base_name"
    samtools flagstat "$bam_file" > "$output_table_dir/${base_name}_alignment_summary.txt"
    echo "Alignment summary written to ${base_name}_alignment_summary.txt"
    echo ""

    #Coverage Statistics
    echo "Coverage Statistics for $base_name"
    samtools depth -a "$bam_file" > "$output_table_dir/${base_name}_coverage_stats.txt"
    echo "Coverage statistics written to ${base_name}_coverage_stats.txt"
    echo ""

    #Mapping Distribution Table
    echo "Mapping Distribution Table for $base_name"
    mapped_reads=$(samtools view -c -F 4 -F 256 "$bam_file")
    unmapped_reads=$(samtools view -c -f 4 "$bam_file")
    total_reads=$(samtools view -c "$bam_file")
    
    echo "Mapped Reads: $mapped_reads" > "$output_table_dir/${base_name}_mapping_distribution.txt"
    echo "Unmapped Reads: $unmapped_reads" >> "$output_table_dir/${base_name}_mapping_distribution.txt"
    echo "Total Reads: $total_reads" >> "$output_table_dir/${base_name}_mapping_distribution.txt"
    echo ""



done

