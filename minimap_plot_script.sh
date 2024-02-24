#!/bin/bash

#Directories
input_dir="/home/swetha/final/refgen"
fastq_dir="/home/swetha/final/meta_data"
output_dir="/home/swetha/final/minimap_plot_output"
result_dir="/home/swetha/final/test_output"
python_files_dir="/home/swetha/final"
output_image_dir="/home/swetha/final/minimap_plot_output"

#Minimap2 reference genome build
for fasta_file in "$input_dir"/*.fasta; do
    #Extract the base name without extension
    base_name=$(basename "$fasta_file" .fasta)

    #Build Minimap2 reference genome for each FASTA file
    minimap2 -d "$output_dir/${base_name}_reference.mmi" "$fasta_file"

    #Perform sequence alignment using Minimap2
    minimap2 -ax sr "$output_dir/${base_name}_reference.mmi" "$fastq_dir/SRR3961906_1.fastq" -o "$output_dir/${base_name}_align.sam"

    #Convert SAM to BAM
    samtools view -b -o "$output_dir/${base_name}_align.bam" "$output_dir/${base_name}_align.sam"

    #Sort the BAM file
    samtools sort -o "$output_dir/sorted_${base_name}_align.bam" "$output_dir/${base_name}_align.bam"

    #Index the sorted BAM file
    samtools index "$output_dir/sorted_${base_name}_align.bam"

    #Execute multiple Python scripts 
    python3 "$python_files_dir/cov_map.py" "$output_dir/sorted_${base_name}_align.bam" "$output_image_dir/${base_name}_cov_map.png"
    python3 "$python_files_dir/coverage_plot.py" "$output_dir/sorted_${base_name}_align.bam" "$output_image_dir/${base_name}_coverage_plot.png"
    python3 "$python_files_dir/ident_perc.py" "$output_dir/sorted_${base_name}_align.bam" "$output_image_dir/${base_name}_ident_perc.png"
done
