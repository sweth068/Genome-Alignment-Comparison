#percent_identity_plot.py

import sys
import pysam
import matplotlib.pyplot as plt

def generate_percent_identity_plot(bam_file, output_plot):
    #Open the BAM file
    bamfile = pysam.AlignmentFile(bam_file, 'rb')

    #Initialize lists to store position and percent identity
    positions = []
    percent_identities = []

    #Iterate through the reads in the BAM file
    for read in bamfile:
        if not read.is_unmapped:  #Ignore unmapped reads
            positions.append(read.reference_start)  #Start position of the read
            percent_identities.append(read.get_tag('NM'))  #Percent identity from NM tag

    #Close the BAM file
    bamfile.close()

    #Create the percent identity plot
    plt.figure(figsize=(10, 6))
    plt.plot(positions, percent_identities, label='Percent Identity', color='blue')
    plt.title('Percent Identity Plot')
    plt.xlabel('Genomic Position')
    plt.ylabel('Percent Identity')
    plt.legend()
    plt.grid(True)

    #Set y-axis limits to ensure it extends up to 100
    plt.ylim(0, 100)

    #Save the plot as an image file
    plt.savefig(output_plot)

 

if __name__ == "__main__":
    #Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python percent_identity_plot.py <input_bam_file> <output_plot_path>")
        sys.exit(1)

    #Extract command-line arguments
    bam_file_path = sys.argv[1]
    output_plot_path = sys.argv[2]

    #Generate percent identity plot and save the image
    generate_percent_identity_plot(bam_file_path, output_plot_path)
