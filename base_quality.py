import pysam
import numpy as np

def calculate_base_quality_distribution(bam_file):
    #Open the BAM file
    bamfile = pysam.AlignmentFile(bam_file, 'rb')

    #Initialize an array to store base qualities
    base_qualities = []

    #Iterate through reads in the BAM file
    for read in bamfile:
        #Ignore unmapped reads
        if not read.is_unmapped:
            #Append base qualities to the array
            base_qualities.extend(read.query_qualities)

    #Close the BAM file
    bamfile.close()

    #Calculate the frequency of each base quality
    quality_counts = np.bincount(base_qualities, minlength=256)

    #Print the table
    print("Base Quality\tFrequency")
    for quality, count in enumerate(quality_counts):
        if count > 0:
            print(f"{quality}\t{count}")

#Calculate base quality distribution
calculate_base_quality_distribution('/home/swetha/final/minimap_plot_output/sorted_myco_12603_align.bam')
