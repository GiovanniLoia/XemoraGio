"""
xr_consensus_modifications.py

Inputs: per-read_modifications.tsv path

Output: file with probability at each analyzed position of modification or no modification 
"""

########################################################################
# Imports 
import os 
import numpy as np
import pandas as pd
from collections import Counter
import math
import sys

########################################################################

# Inputs 
working_dir = os.path.expanduser(sys.argv[1])
per_read_tsv_filepath = os.path.expanduser(sys.argv[2])
bed_file_path = os.path.expanduser(sys.argv[3])

########################################################################

# Making output tsv file path 
output_tsv_filepath = os.path.join(working_dir, 'output_modification_results.tsv')


def filter_read_tsv(read_tsv_path, bed_file_path):
    """
    filter_read_tsv: extract bed file indices and remove rows from the 
    per-read_modifications.tsv where the focus base is not in this range
    
    Parameters:
    read_tsv_path: the per-read_modification.tsv file path
    bed_file_path: the bed file path
    
    Returns: 
    filt_df: filtered dataframe without rows outside of indices
    """
    
    print('Xemora [STATUS] - loading per-read_modifications.tsv and bed file')
    # Read the per-read modifications file
    filt_df = pd.read_csv(read_tsv_path, sep='\t')
    
    # Read the BED file
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])
    
    # Adjust the start positions (subtract 1)
    bed_df['chromStart'] = bed_df['chromStart'] - 1
    bed_df['chromEnd'] = bed_df['chromEnd'] - 2
    
    # Initialize an empty list to store filtered DataFrames
    filtered_dfs = []
    
    print('Xemora [STATUS] - filtering reads in per-read_modifcations.tsv not in bed file range')
    # Iterate through each row in the BED file
    for _, row in bed_df.iterrows():
        start, end = row['chromStart'], row['chromEnd']
        filtered_df = filt_df[(filt_df['read_focus_base'] >= start) & 
                              (filt_df['read_focus_base'] <= end)]
        filtered_dfs.append(filtered_df)
    
    # Concatenate all filtered DataFrames
    filt_df = pd.concat(filtered_dfs, ignore_index=True)
    
    print(filt_df)
    return filt_df

def calculate_metrics(group):
    """
    Calculate various metrics for a group of data.

    Parameters:
    group (DataFrame): A pandas DataFrame group for which to calculate metrics.

    Returns:
    dict: A dictionary containing calculated metrics.
    """
    # Calculate fractions
    total = len(group)
    counts = Counter(group['class_pred'])
    fraction_called_0 = counts[0] / total
    fraction_called_1 = counts[1] / total

    # Calculate final predicted class (mode)
    final_predicted_class = group['class_pred'].mode()[0]

    # Calculate average class probability for each class
    class_probs = group['class_probs'].apply(eval).tolist()
    avg_prob_0 = sum(prob[0] for prob in class_probs) / total
    avg_prob_1 = sum(prob[1] for prob in class_probs) / total

    # Calculate confidence score
    confidence_score = abs(avg_prob_0 - avg_prob_1)

    # Calculate entropy
    entropy = -sum(prob[0] * math.log2(prob[0]) + prob[1] * math.log2(prob[1]) for prob in class_probs if prob[0] > 0 and prob[1] > 0) / total

    # Calculate consensus probability
    consensus_prob = avg_prob_0 if final_predicted_class == 0 else avg_prob_1

    # Calculate prediction variability
    prob_std_dev_0 = pd.Series([prob[0] for prob in class_probs]).std()
    prob_std_dev_1 = pd.Series([prob[1] for prob in class_probs]).std()

    return {
        'fraction_called_dna': fraction_called_0,
        'fraction_called_xna': fraction_called_1,
        'final_predicted_class': final_predicted_class,
        'avg_prob_dna': avg_prob_0,
        'avg_prob_xna': avg_prob_1,
        'confidence_score': confidence_score,
        'entropy': entropy,
        'consensus_prob': consensus_prob,
        'prob_std_dev_dna': prob_std_dev_0,
        'prob_std_dev_xna': prob_std_dev_1
    }

def process_dataframe(data):
    """
    Process the data to calculate metrics for each base of interest.

    Parameters:
    data (DataFrame): A pandas DataFrame containing the input data.

    Returns:
    dict: A dictionary containing the calculated metrics for each base of interest.
    """
    
    results = {}
    base_groups = data.groupby('read_focus_base')

    for base, group in base_groups:
        results[base] = calculate_metrics(group)

    return results

def save_results(results, output_file):
    """
    Save the calculated results to a TSV file.

    Parameters:
    results (dict): The dictionary containing the calculated metrics.
    output_file (str): The path to the output TSV file.
    """
    
    # Convert the results dictionary to a DataFrame
    results_df = pd.DataFrame.from_dict(results, orient='index')
    results_df.index.name = 'read_focus_base'
    results_df.reset_index(inplace=True)

    # Save the DataFrame to a TSV file
    results_df.to_csv(output_file, sep='\t', index=False)

def main(working_dir, per_read_tsv_filepath, bed_file_path, output_tsv_filepath):
    """
    Main function to run the script.

    Parameters:
    working_dir (str): The working directory path.
    per_read_tsv_filepath (str): The path to the input TSV file.
    output_tsv_filepath (str): The path to the output TSV file.
    """
    filt_df = filter_read_tsv(per_read_tsv_filepath, bed_file_path)

    # Print unique read_focus_base values
    unique_bases = filt_df['read_focus_base'].unique()
    print("Unique read_focus_base values:", unique_bases)
    
    # Print sum of occurrences for each unique read_focus_base
    base_counts = filt_df['read_focus_base'].value_counts().reset_index()
    base_counts.columns = ['read_focus_base', 'count']
    print("Sum of occurrences for each read_focus_base:")
    print(base_counts)

    results = process_dataframe(filt_df)
    save_results(results, output_tsv_filepath)

    print(f"Xemora [STATUS] - Consensus modification results saved to {output_tsv_filepath}")

# Run the script
if __name__ == "__main__":
    main(working_dir, per_read_tsv_filepath, bed_file_path, output_tsv_filepath)

