"""
Screen Hochwagen lab's ssDNA ChIP-chip library for chromosome size bias effects.

August 2016

The lab has an extensive library of ssDNA ChIP-chip data, used for mapping DSB hotspots
in meiosis. Not all these data have been collected (or even later examined) to detect
changes in chromosome size bias effects in meiotic recombination. This way, it is worth
screening the entire library for such effects.

Strategy: In order to screen the library, we will use a simple statistical test to compare
meiotic recombination levels between small and large chromosomes and then compare the output
between different samples. The approach is the following:
    - Convert Log2Ratio signal back to ratio
    - Calculate average ChIP-chip signal for each chromosome
    - Consider individual chrs as observations in two groups:
            . small chrs (chr I, III, and VI)
            . large chrs (remaining 13 chrs)
    - Use a statistical test to compare the two groups



For command line usage run:
python chr_size_bias_in_ssDNA_chip-chip.py --help
"""

import time
import pandas as pd
import os
# Import dedicated helper functions from local module
import helpers


# Function to screen microarray experiment data library
def screen_chr_size_bias(path, outfile=False):
    """
    Calculates total signal
    
    Keyword arguments:
    :param path: path to the microarray experiment folders (our microarray data base; no default)
    :param write_to_file: boolean indicating whther to write data to file in cwd (defaults to False)
    :return: Pandas data frame containing all read files
    """
    
    start_time = time.time()
    
    # Get list of all files in the microarray experiment library
    print("Collecting array file names...")
    
    array_name_dict = helpers.array_file_names(path)
    array_path_list = helpers.list_array_paths(path, array_name_dict)
    
    print("Reading array files (with Pandas) from a total of {}:".format(len(array_path_list)))
    
    num_files = 0
    all_array_avrgs_list = []
    
    # Loop through array experiments in folder, Load all array data and collect mean and sd for each chr
    for array_path in array_path_list:
        try:
            print(".", end="")
            # Get data and convert from Log2Ratio back to ratio
            array_df = pd.read_table(array_path)
            array_df['ratio'] = 2**array_df['Log2Ratio']
            # Get pandas df of mean signal for each chromosome
            array_file_name = os.path.split(array_path)[-1]
            exp_folder = os.path.split(os.path.split(array_path)[-2])[-1]
            chr_df = helpers.avrg_array_signal(array_df=array_df, signal_column='ratio',
                                               exp_folder=exp_folder, array_file_name=array_file_name)
            
            all_array_avrgs_list.append(chr_df)
            
            num_files += 1
            if num_files == len(array_path_list):
                print(" Read all files: {}".format(num_files))
            elif num_files % 50 == 0:
                print("({})".format(num_files))
            
        except ValueError:
            print("ERROR: Could not read entry at\n\"{}\"".format(array_path),
                  "\nIs this a valid array analysis file?\nContinuing...\n")
    
    # Calculate signal ratio between small and large chromosomes
    print("Calculating signal ratio for small versus large chromosomes...")
    final_data = []
    for array in all_array_avrgs_list:
        array_data = helpers.small_vs_large(array)
        final_data.append(array_data)
    
    if outfile:
        #file_name = os.getcwd() + '/microarray_chromosome_size_bias.csv'
        file_name = os.getenv("HOME") + '/Desktop/microarray_chromosome_size_bias.csv'
        print("Writing data to '{}'".format(file_name))
        pd.DataFrame(final_data).to_csv(file_name, index=False)
        helpers.print_elapsed_time(start_time)
    else:
        helpers.print_elapsed_time(start_time)
        return(pd.DataFrame(final_data))



def main():
    print()
    print("--------------------------------------------")
    print(" Screen ssDNA ChIP-chip data for changes in")
    print("          chromosome size bias")
    print("--------------------------------------------")
    print()
    screen_chr_size_bias(path=args.path, outfile=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Screen Hochwagen lab's ssDNA ChIP-chip library" +
                                                 " for chromosome size bias effects.")
    parser.add_argument("-p", "--path", default='/Volumes/LabShare/HTGenomics/Microarray_database/arrays/',
                        type=str, help="path to microarray experiment folders (default: '%(default)s)'")
    args = parser.parse_args()
    
    main()
