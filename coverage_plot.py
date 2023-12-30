#coverage_plot.py

import sys
import pysam
import matplotlib.pyplot as plt

def get_coverage(bam_file):
    #Open the BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")

    #Initialize variables to store coverage information
    chromosome_sizes = {ref['SN']: ref['LN'] for ref in samfile.header['SQ']}
    coverage = {chrom: [0] * size for chrom, size in chromosome_sizes.items()}

    #Iterate over each read in the BAM file and update coverage
    for read in samfile.fetch():
        for i in range(read.reference_start, read.reference_end):
            coverage[read.reference_name][i] += 1

    #Close the BAM file
    samfile.close()

    return coverage

def plot_coverage(coverage, output_plot):
    #Plot the coverage for each chromosome
    for chrom, counts in coverage.items():
        plt.plot(range(len(counts)), counts, label=chrom)

    #Customize plot
    plt.title('Genomic Coverage Plot')
    plt.xlabel('Genomic Position')
    plt.ylabel('Coverage')
    plt.legend()
    plt.savefig(output_plot)
    plt.show()

if __name__ == "__main__":
    #Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python coverage_plot.py <input_bam_file> <output_plot_path>")
        sys.exit(1)

    #Extract command-line arguments
    bam_file_path = sys.argv[1]
    output_plot_path = sys.argv[2]

    #Get coverage information
    coverage_data = get_coverage(bam_file_path)

    #Plot coverage
    plot_coverage(coverage_data, output_plot_path)
