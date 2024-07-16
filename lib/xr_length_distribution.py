import pysam
import matplotlib.pyplot as plt

import os
import pysam

def plot_read_length_distribution(bam_file, max_length=500, minimum_length=100, include_long_reads=True, label=None):
    """
    Extract read lengths from a BAM file.
    
    Parameters:
    - bam_file: Path to the BAM file.
    - max_length: Maximum read length to consider (default is 500).
    - minimum_length: Minimum read length to be considered for bad data count.
    - include_long_reads: If True, include reads longer than max_length in the last bin.
    - label: Label for the dataset (for printing purposes).
    
    Returns:
    - read_lengths: List of read lengths.
    - total_reads: Total number of reads.
    - bad_reads: Total number of bad reads (below minimum_length and above max_length).
    """
    
    # Check if the BAM file is indexed
    if not os.path.exists(bam_file + '.bai'):
        try:
            # Generate the BAM index using samtools
            exit_code = os.system(f'samtools index {bam_file}')
            if exit_code != 0:
                raise Exception("samtools index command failed.")
            print(f"Index for {bam_file} created successfully.")
        except Exception as e:
            print(f"Error indexing BAM file: {e}")
            return None, None, None

    try:
        bam = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    except ValueError as e:
        print(f"Error opening BAM file: {e}")
        return None, None, None

    # Initialize a list to store read lengths
    read_lengths = []
    bad_reads = 0
    total_reads = 0
    
    # Iterate through each read in the BAM file
    for read in bam:
        length = read.query_length
        if length <= max_length:
            read_lengths.append(length)
        elif include_long_reads:
            read_lengths.append(max_length)
        total_reads += 1
        if length > max_length or length < minimum_length:
            bad_reads += 1
    
    # Close the BAM file
    bam.close()

    # Print the total number of reads
    print(f"Total number of reads in {label}: {total_reads}")
    print(f"Number of reads over {max_length} nt or below {minimum_length} nt long in {label}: {bad_reads}/{total_reads}")
    
    return read_lengths, total_reads, bad_reads


def plot_combined_read_length_distributions(bam_file_1, bam_file_2, output_file, max_length=500, minimum_length=100, include_long_reads=True, highlight_range=(200, 300)):
    """
    Plot the combined distribution of read lengths from two BAM files and save to a file.
    
    Parameters:
    - bam_file_1: Path to the first BAM file.
    - bam_file_2: Path to the second BAM file.
    - output_file: Path to the output PNG file.
    - max_length: Maximum read length to consider (default is 500).
    - minimum_length: Minimum read length to be considered for bad data count.
    - include_long_reads: If True, include reads longer than max_length in the last bin.
    - highlight_range: Tuple indicating the range to be highlighted in a different color.
    """
    read_lengths_1, total_reads_1, bad_reads_1 = plot_read_length_distribution(bam_file_1, max_length, minimum_length, include_long_reads, label="Batch-1_N44-mer")
    read_lengths_2, total_reads_2, bad_reads_2 = plot_read_length_distribution(bam_file_2, max_length, minimum_length, include_long_reads, label="Batch-50_N44-mer")

    if read_lengths_1 is None or read_lengths_2 is None:
        return

    # Calculate percentages of reads within the highlight range
    percent_within_range_1 = len([length for length in read_lengths_1 if highlight_range[0] <= length <= highlight_range[1]]) / total_reads_1 * 100
    percent_within_range_2 = len([length for length in read_lengths_2 if highlight_range[0] <= length <= highlight_range[1]]) / total_reads_2 * 100

    # Normalize weights by total number of reads
    weights_1 = [1 / total_reads_1] * len(read_lengths_1)
    weights_2 = [1 / total_reads_2] * len(read_lengths_2)

    # Separate data into within and outside the highlight range
    within_range_1 = [length for length in read_lengths_1 if highlight_range[0] <= length <= highlight_range[1]]
    outside_range_1 = [length for length in read_lengths_1 if length < highlight_range[0] or length > highlight_range[1]]

    within_range_2 = [length for length in read_lengths_2 if highlight_range[0] <= length <= highlight_range[1]]
    outside_range_2 = [length for length in read_lengths_2 if length < highlight_range[0] or length > highlight_range[1]]

    # Weights for within and outside range
    weights_within_1 = [1 / total_reads_1] * len(within_range_1)
    weights_outside_1 = [1 / total_reads_1] * len(outside_range_1)

    weights_within_2 = [1 / total_reads_2] * len(within_range_2)
    weights_outside_2 = [1 / total_reads_2] * len(outside_range_2)

    # Plot the distribution of read lengths
    plt.figure(figsize=(14, 6))

    plt.subplot(1, 2, 1)
    plt.hist(outside_range_1, bins=max_length + 1, range=(0, max_length), weights=weights_outside_1, alpha=0.7, color='red', label='Outside 200-300')
    plt.hist(within_range_1, bins=max_length + 1, range=(0, max_length), weights=weights_within_1, alpha=0.7, color='blue', label='Within 200-300')
    plt.title(f'Normalized Distribution of Read Lengths - N27-mer\nTotal Reads: {total_reads_1}\n{percent_within_range_1:.2f}% of reads within 200-300')
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.hist(outside_range_2, bins=max_length + 1, range=(0, max_length), weights=weights_outside_2, alpha=0.7, color='red', label='Outside 200-300')
    plt.hist(within_range_2, bins=max_length + 1, range=(0, max_length), weights=weights_within_2, alpha=0.7, color='blue', label='Within 200-300')
    plt.title(f'Normalized Distribution of Read Lengths - N44-mer\nTotal Reads: {total_reads_2}\n{percent_within_range_2:.2f}% of reads within 200-300')
    plt.xlabel('Read Length')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

# Example usage
bam_file_1 = "/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240617_N44_p2_chunks/N_5000_pod5/canonical/bam/3/final.bam"
bam_file_2 = "/home/marchandlab/github/jay/xemora_randomer_troubleshooting/double_alignment_development/240617_N44_p2_chunks/N_5000_pod5/canonical/bam/48/final.bam"
output_file = "/home/marchandlab/github/jay/xemora_randomer_troubleshooting/visuals_for_randomer_analysis/N44_P2_3_48.png"

plot_combined_read_length_distributions(bam_file_1, bam_file_2, output_file, max_length=500, minimum_length=100, include_long_reads=True, highlight_range=(200, 300))

