import sys
import pysam
import matplotlib.pyplot as plt
import numpy as np

def mapping_quality_vs_coverage(bam_file, output_image_path):
    bam = pysam.AlignmentFile(bam_file, "rb")

    #Initialize lists to store coverage and mapping quality
    coverage = []
    mapping_quality = []

    #Iterate through each reference sequence in the BAM file
    for seq in bam.header['SQ']:
        seq_name = seq['SN']
        seq_length = seq['LN']

        #Initialize arrays for coverage and mapping quality for each position
        coverage_array = np.zeros(seq_length, dtype=int)
        mapping_quality_array = np.zeros(seq_length, dtype=float)

        #Iterate through reads and update coverage and mapping quality arrays
        for read in bam.fetch(seq_name):
            start = read.reference_start
            end = read.reference_end
            mapq = read.mapping_quality

            coverage_array[start:end] += 1
            mapping_quality_array[start:end] += mapq

        #Calculate average mapping quality for each position
        average_mapping_quality = np.divide(mapping_quality_array, coverage_array, 
                                            out=np.zeros_like(mapping_quality_array), 
                                            where=coverage_array!=0)

        #Append values to the overall lists
        coverage.extend(coverage_array)
        mapping_quality.extend(average_mapping_quality)

    #Create a scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(coverage, mapping_quality, alpha=0.5)
    plt.title('Mapping Quality vs. Coverage')
    plt.xlabel('Coverage')
    plt.ylabel('Average Mapping Quality')
    plt.savefig(output_image_path)
    plt.show()

if __name__ == "__main__":
    #Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python mapping_quality_vs_coverage.py <input_bam_file> <output_image_path>")
        sys.exit(1)

    #Extract command-line arguments
    input_bam_file = sys.argv[1]
    output_image_path = sys.argv[2]

    #Call the function with provided arguments
    mapping_quality_vs_coverage(input_bam_file, output_image_path)
