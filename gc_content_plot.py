import pysam
import matplotlib.pyplot as plt

def calculate_gc_content(read):
    gc_count = 0
    total_count = 0

    for base in read.seq:
        if base.upper() in {'G', 'C'}:
            gc_count += 1
        total_count += 1

    return gc_count / total_count if total_count > 0 else 0

def plot_gc_content_bias(bam_file, window_size=1000, output_plot='mini_gc_negative.png'):
    samfile = pysam.AlignmentFile(bam_file, "rb")

    gc_content = []
    positions = []

    for pileupcolumn in samfile.pileup():
        position = pileupcolumn.pos
        gc_sum = 0
        count = 0

        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                gc_sum += calculate_gc_content(pileupread.alignment)
                count += 1

        gc_content.append(gc_sum / count if count > 0 else 0)
        positions.append(position)

    samfile.close()

    plt.plot(positions, gc_content, label="GC Content")
    plt.title('GC Content Bias Plot')
    plt.xlabel('Genomic Position')
    plt.ylabel('GC Content')
    plt.legend(loc='upper right')

    #Save the plot as an image file
    plt.savefig(output_plot)

    #Display the plot
    plt.show()

if __name__ == "__main__":
    #Path to BAM file
    bam_file_path = '/home/swetha/final/minimap_plot_output/sorted_Negative_Thermus_align.bam'

    #Set the window size for plotting (default: 1000)
    window_size = 1000

    #Set the output plot file name
    output_plot_path = 'mini_gc_negative.png'

    #Plot GC content bias and save the image
    plot_gc_content_bias(bam_file_path, window_size, output_plot_path)
