""" Helper functions """

# Standard library
import time
import os
import logging

# Third party libraries
import pandas as pd
from scipy.stats import ttest_ind


def array_file_names(path):
    """
    Given a path to a folder (array library), returns dictionary
    containing each directory (experiment folder) as a key and a
    list of array data files (analysis of array containing Log2Ratio
    data) as its value.
    """
    array_dict = {}
    # Get all directory names in the path (array database)
    for d in os.listdir(path):
        if os.path.isdir(os.path.join(path, d)):
            # Get all array file names (ignore dirs, hidden files, and file names with and extension)
            ffs = []
            for f in os.listdir(os.path.join(path, d)):
                if os.path.isfile(os.path.join(path, d, f)) and not f[0] == '.' and not os.path.splitext(f)[1]:
                    ffs.append(f)
                    array_dict[d] = ffs
    
    return array_dict


def list_array_paths(path, array_dict):
    """
    Given a dictionary containing each directory (experiment folder) as a
    key and a list of array data files (analysis of array containing Log2Ratio
    data) as its value (i.e., the output of 'find_arrays'), returns a list of
    full paths to array files.
    """
    array_path_list = []
    
    for key_folder in array_dict:
        for array_value in array_dict[key_folder]:
            array_path_list.append(path + key_folder + "/" + array_value)
    
    return array_path_list


def avrg_array_signal(array_df, signal_column, exp_folder, array_file_name):
    """
    Given a pandas dataframe containing an array analysis,
    calculates and returns a pandas data frame containing
    total length covered by probes, total signal (typically
    ratio after converted back fromLog2Ratio), and
    total signal/total length for each chromosome.
    """
    # Loop over all chromosomes in turn
    d = []
    for chr_num in array_df['chr'].unique():
        chr_data = array_df.loc[(array_df.chr == chr_num), :]
        # Calculate total signal
        total_signal = array_df.loc[(array_df.chr == chr_num), signal_column].sum()
        mean_signal = array_df.loc[(array_df.chr == chr_num), signal_column].mean()
        # Loop over all rows for each chromosome and calculate
        # sum of probe lengths
        total_length = 0
        for _, row in array_df.loc[(array_df.chr == chr_num), :].iterrows():
            # This length is usualy constant, so I could also just
            # calculate it first and then count probes and multiply by length
            total_length += row['end'] - (row['start'] - 1)
        
        d.append({'exp_folder': exp_folder, 'array': array_file_name,
                  'chr': chr_num, 'total_length': total_length,
                  'total_signal': total_signal,
                  'mean_signal': mean_signal})
    
    return pd.DataFrame(d)


def small_vs_large(chr_df):
    """
    Given a pandas dataframe containing mean signal per chromosome ('chr_df),
    calculates and returns the ratio between mean signal for small versus
    large chromosomes (TO-DO: add t-test / Mann-Whitney test P value calc.).
    """
    # Get means and std
    small = [1, 3, 6]
    small_data = chr_df.loc[chr_df.chr.isin(small), 'mean_signal']
    large_data = chr_df.loc[~chr_df.chr.isin(small), 'mean_signal']
    
    small_avrg = small_data.mean()
    small_std = small_data.std()
    large_avrg = large_data.mean()
    large_std = large_data.std()
    ratio = small_avrg / large_avrg
    
    # Student's t-test (do small chromosomes show higher ssDNA signal?)
    two_tailed_p_val = ttest_ind(small_data, large_data)
    
    #return small_avrg, small_std, large_avrg, large_std, ratio
    return {'exp_folder': chr_df['exp_folder'][0], # names repeated in all 16 lines
            'array': chr_df['array'][0],           # names repeated in all 16 lines
            'small_chr_avrg': small_avrg, 'small_chr_sd': small_std,
            'large_chr_avrg': large_avrg, 'large_chr_sd': large_std,
            'ratio_small_vs_large': ratio, 'ttest_2-tail_p_val': two_tailed_p_val[1]}



def print_elapsed_time(start_time):
    """
    Computes and prints time elapsed since a provided start time.

    Keyword arguments:
    :param start_time: Start time (start_time = time.time()) to compute
    elapsed time from (no default)
    :return: Prints elapsed time since 'start_time'
    """
    elapsed_time = time.time() - start_time
    
    print("\n---")
    if elapsed_time < 60:
        print("Completed in {:2.1f} sec.".format(elapsed_time))
    elif 60 < elapsed_time < 3600:
        print("Completed in {:2.1f} min.".format(elapsed_time / 60))
    else:
        print("Completed in {:2.1f} hr.".format(elapsed_time / 3600))
