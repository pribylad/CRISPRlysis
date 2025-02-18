# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
This script processes CRISPRESSO2 output files to create a master table for CRISPR analysis. 
It includes functions to:
- Write text to a report file.
- Copy files with their parent folders to another location.
- Generate the reverse complement of a DNA sequence.
- Find file paths matching specific criteria.
- Check if a path exists and create it if not.
- Filter out specified values from a DataFrame.
- Save statistical information to an output file.
- Convert Windows-style paths to WSL paths.
- Classify sequences based on counts of deletions, insertions, and mutations.
- Identify ranges of consecutive numbers.
- Count the number of gaps in a sequence.
- Find the closest gap to a cut site.
- Determine the range of a deletion around a cut site.
- Identify microhomology around a deletion.
- Extract insertion events.
- Identify differences between aligned and reference sequences.
- Create a dictionary mapping gRNA sequences to reference sequences.
- Extract sequences around a cut site.
- Calculate the size of a deletion around a cut.
- Generate a summary of gRNA samples and biological replicates.
- Classify deletion events based on their position relative to a cut site.
- Process CRISPRESSO2 output tables and add necessary columns for analysis.
- Create a master table from individual sample tables.
- Add additional columns to the master table.
- Create subtables based on specified parameters.

The script is designed to facilitate the analysis of CRISPR editing events by organizing and processing data from CRISPRESSO2 outputs.

Made by Adela Pribylova, Charles University, Prague, Czech Republic. 2025
"""
##########################################################################################
# To make mastertable from CRISPRESSO2 output files
##########################################################################################
import csv
from datetime import datetime
from pathlib import Path
import os
from os.path import dirname, split, basename
import more_itertools as mit
import pandas as pd

from general_functions import write_into_report, find_files_paths_with_criteria

def filter_out_values_from_df(df, dict_to_filter, project_output_dir='./report.txt'):
    """
    Filters out specified values from a given DataFrame based on the key in the given dictionary.

    Parameters:
        df (pd.DataFrame): DataFrame to be filtered.
        dict_to_filter (dic): Dictionary, where key is the header of the DataFrame and value list of items to be filtered out.
        project_output_dir (str): Path to the report file, default current folder.

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    try:
        for key in dict_to_filter:
            list_to_filter = dict_to_filter[key]
            # Check if the values are as lists, if not, make a list from it.
            if not isinstance(list_to_filter, list):
                write_into_report(f'NOTE: When filtering the data, the ... {list_to_filter} ... from filter_out_values was not given as a list, the code made the conversion. Please, check the output, if all works correctly.', destination_folder=project_output_dir)
                list_to_filter = [list_to_filter]
            # Filter out unwanted values.
            for item in list_to_filter:
                df = df[df[key] != item]
                write_into_report(f'Filtering column {key} with value {item} as {type(item)}. After filter, the shape of the DataFrame is {df.shape}', destination_folder=project_output_dir)
        write_into_report(f'{dict_to_filter} were filtered out from the DataFrame.\nNew shape of the DataFrame is: {df.shape}', destination_folder=project_output_dir)
        return df
    except KeyError as e:
        write_into_report(f'Used dictionary for filterin was: {dict_to_filter}', destination_folder=project_output_dir)
        write_into_report(f'Error: {e}', destination_folder=project_output_dir)

################################################################################
#          FUNCTIONS TO PROCESS THE CRISPRESSO2 ALLES_FREQ_TABLES              #
################################################################################

def clasify_INDELSUB_from_columns(row):
    """
    Classify sequences based on the counts of deletions, insertions, and mutations, including handling cases where gaps occur around the cut site.

    Parameters:
    row (dict): A dictionary representing a row of data with keys 'n_deleted', 'n_inserted', 'n_mutated', and 'Gap_around_cut'.

    Returns:
    str: Classification of the sequence.
    """
    # CRISPResso2 has a mistake in the algorithm and some of the 'Alignes_sequence' with '-' are classified as n_deleted=0
    # in the column 'No_of_gaps' I analysed the sequence by my script to identify DEL around the cut
    n_del, n_ins, n_mut, gap_around_cut = row['n_deleted'], row['n_inserted'], row['n_mutated'], row['Gap_around_cut']
    if ((n_del != 0 and n_ins == 0 and n_mut == 0) or (n_del == 0 and n_ins == 0 and n_mut == 0 and gap_around_cut)):
        return 'DEL'
    elif ((n_del != 0 and n_ins != 0 and n_mut == 0) or (n_del == 0 and n_ins != 0 and n_mut == 0 and gap_around_cut)):
        return 'INDEL'
    elif ((n_del != 0 and n_ins == 0 and n_mut != 0) or (n_del == 0 and n_ins == 0 and n_mut != 0 and gap_around_cut)):
        return 'DELSUB'
    elif ((n_del != 0 and n_ins != 0 and n_mut != 0) or (n_del == 0 and n_ins != 0 and n_mut != 0 and gap_around_cut)):
        return 'INDELSUB'
    elif (n_del == 0 and n_ins != 0 and n_mut == 0):
        return 'IN'
    elif (n_del == 0 and n_ins != 0 and n_mut != 0):
        return 'INSUB'
    elif (n_del == 0 and n_ins == 0 and n_mut != 0):
        return 'SUB'
    elif (n_del == 0 and n_ins == 0 and n_mut == 0):
        return 'no_mutation'
    else:
        raise ValueError('Uncovered_category')

def change_Event_based_on_mut_position(event, mutation_around_cut):
    """
    Check if the mutation event is around the cut site. If not, change it to the correct one.

    Parameters:
    event (str): The current mutation event classification.
    mutation_around_cut (str or bool): Indicates whether the mutation is around the cut site.

    Returns:
    str: Updated mutation event classification.
    """
    if mutation_around_cut == 'cut_out_of_range':
        return mutation_around_cut
    if mutation_around_cut:
        return event
    else:
        new_event = event.replace('SUB','')
        if not new_event:
            new_event = 'no_mutation'
        return new_event

def find_ranges(iterable):
    """
    Yield ranges of consecutive numbers as tuples or single values.

    Parameters:
    iterable (iterable): An iterable containing numbers.

    Yields:
    int or tuple: Single numbers or tuples representing ranges of consecutive numbers.
    """
    for group in mit.consecutive_groups(iterable):
        group = list(group)
        if len(group) == 1:
            yield group[0]
        else:
            yield group[0], group[-1] # Yield a tuple representing the start and end of the range

def count_no_of_gaps(aligned):
    """
    Returns the number of gaps and a list of ranges where gaps occur in the `aligned` string.

    Parameters:
    aligned (str): A string in which gaps are represented by '-'.

    Returns:
    tuple: A tuple containing the number of gaps and a list of ranges.
    """
    deletion_at = [i + 1 for i, char in enumerate(aligned) if char == "-"]

    # Generate ranges from the positions list
    ranges = []
    for r in find_ranges(deletion_at):
        if isinstance(r, int):
            ranges.append((r, r))
        else:
            ranges.append(r)

    return len(ranges), ranges

def find_DEL_closest_to_cut(sequence, cut_position):
    """
    Returns the range of the closest gap to the cut site, whether the gap is within a deletion or not.

    Parameters:
    sequence (str): The sequence to analyze for gaps.
    cut_position (int): The position of the cut site.

    Returns:
    list: A list containing the start and end positions of the closest gap.
    """
    no_of_gaps, ranges_updated = count_no_of_gaps(sequence)

    if no_of_gaps == 0:
        return ["no_gap", "no_gap"]

    def calculate_distance(gap_position):
        distance = abs(cut_position - gap_position)
        return distance + 1 if gap_position <= cut_position else distance

    closest_gap = min(ranges_updated, key=lambda gap: min(calculate_distance(gap[0]), calculate_distance(gap[1])))

    return [closest_gap[0], closest_gap[1]]

def find_min_and_max_gap_ranges(reference, aligned, original_gap_start, original_gap_end, cut_position, cut_maximum, cut_minimum):
    """
    Shifts nucleotides around the deletion to the left and right to find the possible most left and right range of the DEL.
    Determines if the deletion is around the cut, based on the specified cut_position, cut_maximum and cut_minimum.

    Parameters:
    reference (str): The reference sequence.
    aligned (str): The aligned sequence.
    original_gap_start (int or str): Start of the gap determined by CRISPResso2.
    original_gap_end (int or str): End of the gap determined by CRISPResso2.
    cut_position (int): Position of the optimal cut, which was set for CRIPSResso2
    cut_maximum (int): No of nucleotides upstream from cut_position - the most 5' end cut site. Example -2 for SpCas9.
    cut_minimum (int): No of nucleotides upstream from cut_position - the most 3' end cut site. Example 0 for SpCas9.

    Returns:
    tuple: ([min_range], [max_range], bool indicating if DEL is around the specified cut)
    """
    if original_gap_start == 'no_gap':
        return 'no_gap', 'no_gap', False

    if original_gap_start == 1:
        if original_gap_end <(cut_position + cut_maximum):
            return 'gap_span_analysed_window', 'gap_span_analysed_window', False
        return 'gap_span_analysed_window', 'gap_span_analysed_window', True

    if original_gap_end == len(reference):
        if original_gap_start > (cut_position + cut_minimum + 1):
            return 'gap_span_analysed_window', 'gap_span_analysed_window', False
        return 'gap_span_analysed_window', 'gap_span_analysed_window', True

    def shift_to_right(start, end, seq1, seq2):
        counter = 0
        while end + counter < len(seq2):
            if seq2[end + counter] == seq1[start - 1 + counter]:
                counter += 1
            else:
                break
        return [start + counter, end + counter]

    def shift_to_left(start, end, seq1, seq2):
        counter = 0
        while start - 2 - counter >= 0:
            if seq2[start - 2 - counter] == seq1[end - 1 - counter]:
                counter += 1
            else:
                break
        return [start - counter, end - counter]

    max_range = shift_to_right(original_gap_start, original_gap_end, reference, aligned)
    min_range = shift_to_left(original_gap_start, original_gap_end, reference, aligned)

    # Analyse if the DEL is around the cut
    if ((min_range[1] < cut_position + cut_maximum) or (min_range[0] > cut_position + cut_minimum + 1)) and ((max_range[1] < cut_position + cut_maximum) or (max_range[0] > cut_position + cut_minimum + 1)):
        return min_range, max_range, False
    return min_range, max_range, True

def find_DEL_around_cut_range(Gap_closest_to_cut_range, Gap_around_cut):
    """
    If the gap is around the cut site, returns the start and end of the deletion.

    Parameters:
    Gap_closest_to_cut_range (list): The range of the closest gap to the cut site.
    Gap_around_cut (bool): Whether the gap is around the cut site.

    Returns:
    tuple: Start and end of the deletion if around cut, otherwise "gap_out_of_cut".
    """
    return (Gap_closest_to_cut_range[0], Gap_closest_to_cut_range[1]) if Gap_around_cut else (False, False)


def gap_span_analysed_window(gap_range, aligned_seq):
    """
    Check if the gap spans the analyzed window. If the gap starts or ends with '-', those sequences cannot be analyzed for micro-homologies.

    Parameters:
    gap_range (list): A list containing the start and end positions of the gap.
    aligned_seq (str): The aligned sequence to analyze.

    Returns:
    bool: True if the gap spans the analyzed window, False otherwise.
    """
    return gap_range[0] == 1 or gap_range[1] == len(aligned_seq)


def find_microhomology(gap_start, gap_end, gap_span_analysed_window, aligned, reference):
    """ 
    Searches the maximum microhomology around the deletion.

    Parameters:
        gap_start (int): position of the start of the gap
        gap_end (int): position of the end of the gap
        gap_span_analysed_window (boolean): information if the gap span the analysed window
        aligned (str): sequence of the Aligned_Sequence from CRISPResso2 output
        reference (str): sequence of the Reference_Sequence from CRISPResso2 output

    Returns: 
        microhomology_to_right_lenght (int): how many nt on the left side of the gap has microhomology with the right side od the gap
        microhomology_to_left_lenght (int): how many nt on the left side of the gap has microhomology with the right side od the gap
        microhomology_lenght (int): lenght of the maximal found microhomology, joined from the left and right
        microhomology (str): sequence of the maximal found microhomology, joined from the left and right"""
    if gap_span_analysed_window:
        statement = 'gap_out_of_cut'
        return statement, statement, statement, statement
    if gap_start == 'gap_out_of_cut':
        statement = 'gap_out_of_cut'
        return statement, statement, statement, statement
    if gap_start == 'no_gap':
        statement = 'no_gap'
        return statement, statement, statement, statement

    microhomology_to_right_lenght = 0
    microhomology_to_left_lenght = 0
    size_of_gap = gap_end-gap_start+1

    ##################################################
    ##### HOW MANY nt CAN BE SHIFT TO THE RIGHT? #####
    # for 1bp deletions
    if size_of_gap == 1:
        #      0123456789
        # A -> NNC-NNNNNN (gap_start = 4, gap_end = 4)
        # R -> NNCCNNNNNN
        if aligned[gap_start-2] == reference[gap_end-1]:
            microhomology_to_right_lenght += 1
    # for >1bp deletions
    else:
        if size_of_gap != 0:
            #      0123456789
            # A -> NGC---NNNN (gap_start = 4, gap_end = 6, size_of_gap = 3)
            # R -> NGCAGCNNNN
            for s in range(size_of_gap):
                if aligned[gap_start-2-s] == reference[gap_end-1-s]:
                    microhomology_to_right_lenght += 1
                else:
                    break

    ##################################################
    ##### HOW MANY nt CAN BE SHIFT TO THE LEFT?  #####
    # for 1bp deletions
    if size_of_gap == 1:
        #      0123456789
        # A -> NN-CNNNNNN (gap_start = 3, gap_end = 3)
        # R -> NNCCNNNNNN
        if aligned[gap_end] == reference[gap_end-1]:
            microhomology_to_left_lenght += 1
    # for >1bp deletions
    else:
        if size_of_gap != 0:
            #      0123456789
            # A -> NNN---AGAN (gap_start = 4, gap_end = 6, size_of_gap = 3)
            # R -> NNNAGCNNNN
            for s in range(size_of_gap):
                if gap_end+s < len(aligned):
                    if aligned[gap_end+s] == reference[gap_start-1+s]:
                        microhomology_to_left_lenght += 1
                    else:
                        break

    ##################################################
    ##### DEFINE MICROHOMOLOGY FROM BOTH SIDES   #####
    microhomology = reference[(gap_start-1-microhomology_to_right_lenght):(gap_start-1+microhomology_to_left_lenght)]
    microhomology_lenght = microhomology_to_right_lenght + microhomology_to_left_lenght
    # if the microhomology would be from both sides and together forms longer sequence than size of the gap, than:
    if (microhomology_to_right_lenght + microhomology_to_left_lenght) >= size_of_gap:
        microhomology = microhomology[0:size_of_gap]
        microhomology_lenght = size_of_gap

    #####################################################################
    ##### EDIT MICROHOMOLOGY TO LENGHT FOR ANALYSED TYPES OF CUTS   #####
    return microhomology_to_right_lenght, microhomology_to_left_lenght, microhomology_lenght, microhomology

def get_insertion_event(row, cut_position):
    """
    Extract the insertion event sequence if there is a PAM-proximal gap next to the cut site.

    Parameters:
    row (dict): A dictionary containing sequence information.
    cut_position (int): The position of the cut site.

    Returns:
    str or None: The insertion event sequence if present, otherwise None.
    """
    if row['Ref_seq_has_PAM_proximal_gap_next_to_cut_site'] == 1:
        n = row['n_inserted']
        insertion_event = row['Aligned_Sequence'][cut_position:cut_position+n]
    else:
        insertion_event = None
    return insertion_event

def gap_in_aligned_seq(row, position):
    """
    Determine if there's a gap in the aligned sequence at the specified position.

    Parameters:
    row (dict): A dictionary containing sequence information.
    position (int): The position to check for a gap in the aligned sequence.

    Returns:
    str or bool: A message indicating an inability to determine if there's a gap if the reference sequence contains a gap,
                 otherwise a boolean indicating if there's a gap at the specified position.
    """
    if '-' in row['Reference_Sequence']:
        return 'cannot determine - Ref.seq with gap'
    else:
        return row['Aligned_Sequence'][position] == '-'

def gap_in_range(row, column, lower_range, higher_range):
    """
    Determine if there's a gap in the specified range of the aligned sequence.

    Parameters:
    row (dict): A dictionary containing sequence information.
    column (str): The key for accessing the relevant sequence in the row.
    lower_range (int): The lower bound of the range to check (1-based index).
    higher_range (int): The higher bound of the range to check.

    Returns:
    str or bool: A message indicating an inability to determine if there's a gap if the reference sequence contains a gap, otherwise a boolean indicating if there's a gap in the specified range.
    """
    if '-' in row['Reference_Sequence']:
        return 'cannot determine - Ref.seq with gap'
    else:
        return '-' in row[column][lower_range-1:higher_range]

def find_difference(aligned, reference):
    """
    Identify differences between aligned and reference sequences.

    Parameters:
    aligned (str): The aligned sequence.
    reference (str): The reference sequence.

    Returns:
    tuple: A list of substitution bases and their positions.
    """
    SUB_list = []
    SUB_pos_list = []

    for i, (a, r) in enumerate(zip(aligned, reference), start=1):
        if a != r and a != '-' and r != '-':
            SUB_list.append(a)
            SUB_pos_list.append(i)

    return SUB_list, SUB_pos_list

def mutation_in_distance(cut_position, list_of_mutation_position, distance_from_cut):
    """
    Determine if any mutations are within a specified distance from the cut site.

    Parameters:
    cut_position (int): The position of the cut site.
    list_of_mutation_position (list): A list of positions where mutations occur.
    distance_from_cut (int): The maximum distance from the cut site to check for mutations.

    Returns:
    bool: True if any mutation is within the specified distance from the cut site, False otherwise.
    """
    for position in list_of_mutation_position:
        if (cut_position - distance_from_cut) <= position <= (cut_position + distance_from_cut):
            return True
    return False

def make_original_seq_dic(Reference_Sequence, Unedited, gRNA_seq):
    """
    Create a dictionary mapping gRNA sequences to their reference sequences from the given data.

    Parameters:
    Reference_Sequence (list): List of reference sequences.
    Unedited (list of bool): List indicating whether each sequence is unedited.
    gRNA_seq (list): List of gRNA sequences.

    Returns:
    dict: Dictionary mapping gRNA sequences to reference sequences.
    """
    df = pd.DataFrame({'Reference_Sequence': Reference_Sequence, 'Unedited': Unedited, 'gRNA_seq': gRNA_seq})
    # !! DO NOT use df[df['Unedited'] is True] as the code will not work correctly
    df = df[df['Unedited'] == True].drop_duplicates(subset=['gRNA_seq'])
    original_seq_dic = df.set_index('gRNA_seq')['Reference_Sequence'].to_dict()
    return original_seq_dic

def make_seq_around_optimum_cut(Original_Sequence, cut_position, distance):
    """
    Extracts sequences proximal and distal to the blunt cut from the original sequence.

    Parameters:
    Original_Sequence (str): The original nucleotide sequence.
    cut_position (int): The position of the cut site.
    distance (int): Number of nucleotides to extract on either side of the cut.

    Returns:
    str: The sequence around the blunt cut with a '|' separator at the cut site.
    """
    seq = Original_Sequence[cut_position - distance:cut_position] + '|' + Original_Sequence[cut_position:cut_position + distance]
    return seq

def size_of_del_around_the_cut (gap_around_cut_starts, gap_around_cut_ends):
    """
    Calculate the size of a deletion around the cut.

    Parameters:
    gap_around_cut_starts (int or str): Start position of the gap, or a string indicating no gap.
    gap_around_cut_ends (int): End position of the gap.

    Returns:
    int or False: Size of the deletion. Returns False if there is no gap or it is not around cut.
    """
    if gap_around_cut_starts is False:
        return False
    else:
        return gap_around_cut_ends - gap_around_cut_starts + 1

def no_of_gRNA_samples_txt(df, order_list):
    """
    Generate a summary string of the number of gRNA samples and biological replicates.

    Parameters:
    df (DataFrame): DataFrame containing the data, including a "replicate" column.
    order_list (list): List of gRNA samples to be analyzed.

    Returns:
    Nothing
    """
    no_of_gRNA_samples = f'Analysed {len(order_list)} gRNAs in {len(df["replicate"].unique())} biological replicates.'
    write_into_report(no_of_gRNA_samples)
    return no_of_gRNA_samples

def return_shifted_DEL_range(cut_position, DEL_range):
    """Follow-up functions as return_category_for_DEL and return_cut_type_based_on_cut_categories needs for correct classifications cut_optimum==0. For that, this function shifts the ranges.
    Parameters:
    cut_position (int): The position of the cut site.
    DEL_range (list): is the x, y coordinates of DELetion
    
    RETURNS:
    int: shifted x coordinates
    int: shifted y coordinates"""
    # assign x and y coordinate
    x=DEL_range[0]
    y=DEL_range[1]

    # Check if the values in x and y are strings (likely 'no_gap')
    if isinstance(x, str) or isinstance(x, str):
        return x, x

    # Solve walking around 0
    if x - cut_position > 0:
        x_new = x - cut_position
    else:
        x_new = x - cut_position - 1
    if y - cut_position > 0:
        y_new = y - cut_position
    else:
        y_new = y - cut_position - 1

    print(f'New DEL_range is ({x_new},{y_new}).')
    return x_new, y_new

def return_category_of_DEL(cut_position, cut_maximum, cut_minimum, DEL_range):
    """This function returns the cut's category based on the allowed maximum, optimum, and minimum position of the cut and based on the DELetion range.
    Parameters:
    cut_position (int): The position of the cut site.
    cut_maximum (int): is the most 5' position where the Cas can cut
    cut_minimum (int): is the most 3' postition where the Cas can cut
    DEL_range (list): is the x, y coordinates of DELetion

    Returns:
    int: category of the cut
    str: comment of the category"""
    # assign x and y coordinate
    x=DEL_range[0]
    y=DEL_range[1]

    # Check if the values in x and y are strings (likely 'no_gap')
    if isinstance(x, str) or isinstance(y, str):
        return 97, 'x or y coordinates is string, likely no DELetion.'

    # Categories bellow works for cut_optimum = 0, shift everything so cut_optimum is 0, cut_maximum < 0 and cut_min >= 0
    x, y = return_shifted_DEL_range(cut_position, DEL_range)
    if x==0 or y==0:
        return 98, 'x or y coordinate is 0. Problem with determining the range of DELetion.'

    cut_optimum = 0

    # General comment
    comment = 'Range is not classified yet.'
    category = 0

    # Set category based on the logic decision.
    if x>y:
        comment = 'x<y which should not exists.'
        category = 99
        return category, comment

    if y<cut_maximum-1:
        comment = 'DEL is outside the cut upstream.'
        category = 1
        return category, comment

    if y==cut_maximum-1:
        comment = 'DEL is upstream and next to the cut_maximum.'
        category = 2
        return category, comment

    if cut_maximum<=y<cut_optimum-1:
        comment = 'DEL is between cut_maximum and cut_optimum-2 -> not a blunt cut.'
        category = 3
        return category, comment

    if y==cut_optimum-1:
        comment = 'DEL is upstream and next to the cut_optimum -> blunt cut.'
        category = 4
        return category, comment

    if x<cut_maximum and cut_optimum<y<= cut_minimum:
        comment = 'DEL starts upstream from cut_maximum end ends between cut_optimum and cut_maximum.'
        category = 5
        return category, comment

    if x<cut_maximum and cut_minimum<y:
        comment = 'DEL span all possible cut events -> both side DEL.'
        category = 6
        return category, comment

    if cut_maximum<=x<cut_optimum and cut_optimum<y:
        comment = 'DEL starts betweem cut_maximum and cut_optimum and ends after cut_optimum.'
        category = 7
        return category, comment

    if x==cut_optimum+1:
        comment = 'DEL starts downstream from the cut_optimum.'
        category = 8
        return category, comment

    if cut_optimum+1<x<cut_minimum:
        comment = 'DEL starts between Optimum and cut_minimum, but not adjacent to them.'
        category = 9
        return category, comment

    if x==cut_minimum:
        comment = 'DEL strats downstream adjacent to the cut_minimum.'
        category = 10
        return category, comment

    if x==cut_minimum+1:
        comment = 'DEL starts upstream adjacent to the cut_minimum.'
        category = 11
        return category, comment

    if cut_minimum+1<x:
        comment = 'DEL is outside the cut downstream.'
        category = 12
        return category, comment

    return category, comment

def return_cut_type_based_on_cut_categories(category_1, range_1, category_2, range_2, cut_maximum, cut_minimum):
    """Returns the predicted cleavage type based on the default cleavage category for DELetion range 1 and 2, where all nucleotides were 'stacked' to the right and then to the left, respectively.
    
    Parameters:
    category_1 (int): category for the first range based on the return_category_of_DEL, where the shiftable-nucleotides were moved to the right in function return_category_of_DEL
    range_1 (list): range of the DELetion of category_1 processed by return_shifted_DEL_range
    category_2 (int): category for the second range based on the return_category_of_DEL, where the shiftable-nucleotides were moved to the left in function return_category_of_DEL
    range_2 (list): range of the DELetion of category_2 processed by return_shifted_DEL_range
    cut_maximum (int): is the most 5' position where the Cas can cut
    cut_minimum (int): is the most 3' postition where the Cas can cut

    Returns:
    str: Cut type category
    str: End participating on the microhomology orientoriented according to the cleavage site
    str: Orientation of the DELetion according to the cleavage site
    int: DNA end participating in repair (likely) - 3', 5', 0
    bool: If the Cut type categori is almost sure (0) or likely to be based on the homology lenght (1)
    """
    cut_optimum = 0
    # Check if the values in x1 and x2 are strings (likely 'no_gap')
    if isinstance(range_1[0], str) or isinstance(range_2[0], str):
        return 'X', 'X', 'X', 0, 0

    # Define list of outputs
    list_1 = ['O','X','U',0,0]
    list_2 = [range_1[1]+1, 'X', 'U', 0, 0]
    list_3 = [range_1[1], 'X', 'U', 0, 0]
    list_4 = ['C', 'X', 'B', 0, 0]
    list_5 = ['X', 'X', 'D', 0, 0]
    list_6 = [range_1[0]-1, 'X', 'D', 0, 0]
    list_7 = ['O','X','D',0,0]
    list_8 = ['O','C','U',0,0]
    list_9 = [cut_maximum,'C','U',0,0]
    list_10 = [cut_maximum,'D','U',5,1]
    list_11 = ['C','C','C',0,0]
    list_12 = [range_2[0]-1,'U','D',3,1]
    list_13 = [cut_minimum,'U','D',3,1]
    list_14 = [range_1[1]+1,'D','U',5,1]
    list_15 = [cut_optimum,'C','C',0,0]
    list_16 = [range_1[1],'D','U',5,1]
    list_17 = [range_1[1],'U','D',3,1]
    list_18 = ['C','C','C',0,1]
    list_19 = [range_2[0],'U','D',3,1]
    list_20 = ['O','C','D',0,0]

    ########################################################################################
    # Based on the position of the cut_maximum and cut_minimum decide which list is going to be used
    if abs(cut_maximum)>abs(range_2[0])-1:
        list_a = list_10
    elif abs(cut_maximum)==abs(range_2[0])-1:
        list_a = list_11
    else:
        list_a = list_12

    if abs(cut_maximum)>abs(cut_minimum):
        list_b = list_10
    elif abs(cut_maximum)==abs(cut_minimum):
        list_b = list_11
    else:
        list_b = list_13

    if abs(range_1[1])-1>abs(range_2[0])-1:
        list_c = list_14
        list_e = []
    elif abs(range_1[1])-1==abs(range_2[0])-1:
        list_c = list_11
        list_e = list_15
    else:
        list_c = list_12
        list_e = list_12

    if abs(range_1[1])-1>abs(cut_minimum):
        list_d = list_14
    elif abs(range_1[1])-1==abs(cut_minimum):
        list_d = list_11
    else:
        list_d = list_13

    if abs(range_2[0])-cut_optimum>abs(range_1[1])-cut_optimum:
        list_f = list_16
    elif abs(range_2[0])-cut_optimum==abs(range_1[1])-cut_optimum:
        list_f = list_11
    else:
        list_f = list_19

    if range_1[1]>=range_2[0]:
        list_g = list_12
    else:
        list_g = list_17

    if range_1[1]<range_2[0]:
        list_h = list_17
    else:
        list_h = list_12
    ########################################################################################

    # Determine for DELetions without shiftable nucleotides - no homology around the cut site
    if range_1 == range_2:
        dict_NoMH={1:list_1,
                    2:list_2,
                    3:list_2,
                    4:list_2,
                    5:list_3,
                    6:list_4,
                    7:list_5,
                    8:list_6,
                    9:list_6,
                    10:list_6,
                    11:list_6,
                    12:list_7}
        cut_type = dict_NoMH[category_1][0]
        end_participating = dict_NoMH[category_1][1]
        DEL_orientation = dict_NoMH[category_1][2]
        strand_participating = dict_NoMH[category_1][3]
        likely_type = dict_NoMH[category_1][4]
        return cut_type, end_participating, DEL_orientation, strand_participating, likely_type
    if category_1 == 1:
        dict_NoMH={1:list_8,
                   2:list_9,
                   3:list_10,
                   4:list_10,
                   5:list_10,
                   6:list_10,
                   7:list_a,
                   8:list_a,
                   9:list_a,
                   10:list_a,
                   11:list_a,
                   12:list_b}
    elif category_1 == 2:
        dict_NoMH={3:list_14,
                   4:list_14,
                   5:list_14,
                   6:list_14,
                   7:list_c,
                   8:list_c,
                   9:list_c,
                   10:list_c,
                   11:list_c,
                   12:list_d}
    elif category_1 == 3:
        dict_NoMH={3:list_14,
                   4:list_14,
                   5:list_14,
                   6:list_14,
                   7:list_14,
                   8:list_c,
                   9:list_c,
                   10:list_c,
                   11:list_c,
                   12:list_d}
    elif category_1 == 4:
        dict_NoMH={5:list_15,
                   6:list_15,
                   7:list_15,
                   8:list_e,
                   9:list_e,
                   10:list_e,
                   11:list_e,
                   12:list_13}
    elif category_1 == 5:
        dict_NoMH={5:list_16,
                   6:list_16,
                   7:list_11,
                   8:list_15,
                   9:list_12,
                   10:list_17,
                   11:list_17,
                   12:list_17}
    elif category_1 == 6:
        dict_NoMH={6:list_18,
                   7:list_19,
                   8:list_15,
                   9:list_12,
                   10:list_12,
                   11:list_12,
                   12:list_13}
    elif category_1 == 7:
        dict_NoMH={7:list_f,
                   8:list_15,
                   9:list_g,
                   10:list_g,
                   11:list_g,
                   12:list_13}
    elif category_1 in (8, 9, 10, 11):
        dict_NoMH={9:list_h,
                   10:list_h,
                   11:list_h,
                   12:list_13}
    elif category_1 == 12:
        dict_NoMH={12:list_20}
    else:
        raise ValueError(f'Category_1 from return_category_of_DEL is not in function return_cut_type_based_on_cut_categories. Value of category_1 is {category_1}')

    cut_type = dict_NoMH[category_2][0]
    end_participating = dict_NoMH[category_2][1]
    DEL_orientation = dict_NoMH[category_2][2]
    strand_participating = dict_NoMH[category_2][3]
    likely_type = dict_NoMH[category_2][4]

    return cut_type, end_participating, DEL_orientation, strand_participating, likely_type

# Add information about the sample, guide and its selected nucleoids into the table as new columns
def process_table(table, gRNA_seq, sample, parameters_dic):
    """Adds to the given table many columns with stuff which are needed for CRISPRlysis."""
    df = pd.DataFrame(table)
    mutation_diff = df.apply(lambda row: find_difference(row["Aligned_Sequence"], row["Reference_Sequence"]), axis=1)
    df['sample'] = sample
    df['gRNA_seq'] = gRNA_seq
    df['Cut_position'] = parameters_dic['cut_position']
    df['Nuclease'] = parameters_dic['nuclease']
    # make dictionary with original sequences according to the gRNA_seq
    original_seq_dic = make_original_seq_dic(df['Reference_Sequence'], df['Unedited'], df['gRNA_seq'])
    print(original_seq_dic)
    # check if key is present in the dictionary, if no, use 'Reference_sequence' instead
    df['Original_Sequence'] = df.apply(lambda row: original_seq_dic[gRNA_seq] if (gRNA_seq in original_seq_dic.keys()) else row['Reference_Sequence'], axis=1)
    df['Original_Sequence_around_cut'] = df.apply(lambda row: make_seq_around_optimum_cut(row['Original_Sequence'], parameters_dic['cut_position'], 20), axis=1)
    df['Aligned_Sequence_around_cut'] = df.apply(lambda row: make_seq_around_optimum_cut(row['Aligned_Sequence'], parameters_dic['cut_position'], 20), axis=1)
    df['nt_around_cut'] = df.apply(lambda row: row['Original_Sequence'][parameters_dic['cut_position']-1:parameters_dic['cut_position']+1], axis=1)
    if parameters_dic['nuclease'] == 'SpCas9':
        df['PAM'] = df.apply(lambda row: row['Original_Sequence'][parameters_dic['cut_position']+3:parameters_dic['cut_position']+6], axis=1)
    elif parameters_dic['nuclease'] == 'LbCas12a':
        df['PAM'] = df.apply(lambda row: row['Original_Sequence'][parameters_dic['cut_position']-27:parameters_dic['cut_position']-23], axis=1)
    else:
        df['PAM'] = f'Unknown nuclease - {parameters_dic["nuclease"]}'
    df['3rd_nt'] = df.apply(lambda row: row['Original_Sequence'][parameters_dic['cut_position']], axis=1)
    df['4th_nt'] = df.apply(lambda row: row['Original_Sequence'][parameters_dic['cut_position']-1], axis=1)
    df['5th_and_4th'] = df.apply(lambda row: row['Original_Sequence'][parameters_dic['cut_position']-2:parameters_dic['cut_position']], axis=1)
    df['No_of_gaps'] = df.apply(lambda row: count_no_of_gaps(row['Aligned_Sequence'])[0], axis=1)
    df['Gap_ranges'] = df.apply(lambda row: count_no_of_gaps(row['Aligned_Sequence'])[1], axis=1)
    df['Gap_closest_to_cut_range'] = df.apply(lambda row: find_DEL_closest_to_cut(row['Aligned_Sequence'], parameters_dic['cut_position']), axis=1)
    df['Gap_min_max_ranges'] = df.apply(lambda row: find_min_and_max_gap_ranges(\
                                                    reference = row['Reference_Sequence']\
                                                    ,aligned=row['Aligned_Sequence']\
                                                    ,original_gap_start=row['Gap_closest_to_cut_range'][0]\
                                                    ,original_gap_end=row['Gap_closest_to_cut_range'][1], cut_position=parameters_dic['cut_position'], cut_maximum=parameters_dic['cut_maximum'], cut_minimum=parameters_dic['cut_minimum'])[0:2], axis=1)
    df['Gap_around_cut'] = df.apply(lambda row: find_min_and_max_gap_ranges(\
                                                    reference = row['Reference_Sequence']\
                                                    ,aligned=row['Aligned_Sequence']\
                                                    ,original_gap_start=row['Gap_closest_to_cut_range'][0]\
                                                    ,original_gap_end=row['Gap_closest_to_cut_range'][1], cut_position=parameters_dic['cut_position'], cut_maximum=parameters_dic['cut_maximum'], cut_minimum=parameters_dic['cut_minimum'])[2], axis=1)
    df['Gap_around_cut_starts'] = df.apply(lambda row: find_DEL_around_cut_range(row['Gap_closest_to_cut_range'], row['Gap_around_cut'])[0], axis=1)
    df['Gap_around_cut_ends'] = df.apply(lambda row: find_DEL_around_cut_range(row['Gap_closest_to_cut_range'], row['Gap_around_cut'])[1], axis=1)
    df['n_deleted_around_cut'] = df.apply(lambda row: size_of_del_around_the_cut(row['Gap_around_cut_starts'], row['Gap_around_cut_ends']), axis=1)
    df['Gap_span_analysed_window'] = df.apply(lambda row: gap_span_analysed_window(\
                                                    row['Gap_closest_to_cut_range']\
                                                    , row['Aligned_Sequence']), axis=1)

    def apply_find_microhomology(df, func, new_column_names):
        """Adds new columns into table based on the function find_microhomology."""
        results = df.apply(lambda row: func(row['Gap_closest_to_cut_range'][0],
                                            row['Gap_closest_to_cut_range'][1],
                                            row['Gap_span_analysed_window'],
                                            row['Aligned_Sequence'],
                                            row['Reference_Sequence']), axis=1)
        # Unpack results into separate columns
        for i, col_name in enumerate(new_column_names):
            df[col_name] = [result[i] for result in results]
    # Call the above function to add columns based on the function find_microhomology
    apply_find_microhomology(df, find_microhomology, ['MH_to_proximal_len', 'MH_to_distal_len', 'Maximal_MH_len', 'Maximal_MH'])

    def apply_return_cut_type_based_on_cut_categories(df, func, new_column_names):
        """Adds new columns into table based on the function return_cut_type_based_on_cut_categories."""
        results = df.apply(lambda row: func(\
        return_category_of_DEL(parameters_dic['cut_position'], parameters_dic['cut_maximum'], parameters_dic['cut_minimum'], row['Gap_min_max_ranges'][0])[0],\
        return_shifted_DEL_range(parameters_dic['cut_position'], row['Gap_min_max_ranges'][0]),\
        return_category_of_DEL(parameters_dic['cut_position'], parameters_dic['cut_maximum'], parameters_dic['cut_minimum'], row['Gap_min_max_ranges'][1])[0],\
        return_shifted_DEL_range(parameters_dic['cut_position'], row['Gap_min_max_ranges'][1]),\
        parameters_dic['cut_maximum'],\
        parameters_dic['cut_minimum']), axis=1)
        # Unpack results into separate columns
        for i, col_name in enumerate(new_column_names):
            df[col_name] = [result[i] for result in results]

    df['Gap_category_1'] = df.apply(lambda row: return_category_of_DEL(parameters_dic['cut_position'], parameters_dic['cut_maximum'], parameters_dic['cut_minimum'], row['Gap_min_max_ranges'][0])[0], axis=1)
    df['Gap_category_2'] = df.apply(lambda row: return_category_of_DEL(parameters_dic['cut_position'], parameters_dic['cut_maximum'], parameters_dic['cut_minimum'], row['Gap_min_max_ranges'][1])[0], axis=1)

    # Call the above function to add columns based on the function find_microhomology
    apply_return_cut_type_based_on_cut_categories(df, return_cut_type_based_on_cut_categories, ['Cut_type', 'End_participating_on_MH', 'DEL_orientation', 'End_3_5', 'Likely_cut_type'])

    # TODO once the code is working, make the function for the following lines
    # df['MH_around_cut'] = df.apply(lambda row: MH_around_cut(row['DEL_orientation_min_max_range'], row['Maximal_MH'], row['Aligned_Sequence'], row['End_participating_on_MH'], parameters_dic['cut_position']), axis=1)
    # df['MH_around_cut_len'] = df.apply(lambda row: len(row['MH_around_cut']), axis=1)

    df['Ref_seq_has_PAM_proximal_gap_next_to_cut_site'] = df.apply(lambda row: row['Reference_Sequence'][parameters_dic['cut_position']] == '-', axis=1)
    df['Insertion_event'] = df.apply(lambda row:  get_insertion_event(row, parameters_dic['cut_position']), axis=1)
    df['4th_nt_inserted'] = df.apply(lambda row: row['Insertion_event'] == row['4th_nt'], axis=1)
    df['5th_and_4th_nt_inserted'] = df.apply(lambda row: row['Insertion_event'] == row['5th_and_4th'], axis=1)
    df['Aligned_nt_cut_upstream'] = df.apply(lambda row:  row['Aligned_Sequence'][parameters_dic['cut_position']-1], axis=1)

    # it is around middle, which does not means that it is around the cut, as for sequences with insertions there is nt shift
    df['Aligned_2nt_mid'] = df.apply(lambda row:  row['Aligned_Sequence'][parameters_dic['cut_position']-1:parameters_dic['cut_position']+1], axis=1)
    df['Reference_2nt_mid'] = df.apply(lambda row:  row['Reference_Sequence'][parameters_dic['cut_position']-1:parameters_dic['cut_position']+1], axis=1)
    df['Same_2nt_mid'] = df.apply(lambda row:  row['Aligned_2nt_mid'] == row['Reference_2nt_mid'], axis=1)

    # TODO Do I need the followitn four columns?
    # # this works only if there is no '-' in Reference Sequence, thus no IN event around predicted cut, else 'cannot determine'
    # df['Aligned_sequence_has_gap_next_to_blunt_cut_in_PAM_proximal'] = df.apply(lambda row:  gap_in_aligned_seq(row, parameters_dic['cut_position']), axis=1)
    # df['Aligned_sequence_has_gap_next_to_blunt_cut_in_PAM_distal'] = df.apply(lambda row:  gap_in_aligned_seq(row, parameters_dic['cut_position']-1), axis=1)

    # # This does not count with staggered end formation and microhomology, just DEL in first half or second half of Aligned sequence
    # df['Aligned_sequence_contain_gap_distal_from_blunt_cut'] = df.apply(lambda row:  gap_in_range(row, 'Aligned_Sequence', 0, parameters_dic['cut_position']), axis=1)
    # df['Aligned_sequence_contain_gap_proximal_from_blunt_cut'] = df.apply(lambda row:  gap_in_range(row, 'Aligned_Sequence', parameters_dic['cut_position'], (len(row['Aligned_Sequence'])+1)), axis=1)

    df['Referene_sequence_contain_gap'] = df.apply(lambda row: '-' in row['Reference_Sequence'], axis=1)

    # # to know the letter of the mutation
    df['mutation'] = mutation_diff.apply(lambda row: row[0])
    # # to know the position of the mutation
    df['mutation_pos'] = mutation_diff.apply(lambda row: row[1])
    # # to know, if the mutation is in given distance from the cut -> now it is 3, but play with the number to see if it works
    df['mutation_around_cut'] = df.apply(lambda row:  mutation_in_distance(parameters_dic['cut_position'], row['mutation_pos'], 3), axis=1)
    # Add column with edited event format -> if the mutation in not around the cut. Use for example IN instead of INSUB.
    df['Event'] = df.apply(lambda row: clasify_INDELSUB_from_columns(row), axis=1)
    df['Event_around_cut'] = df.apply(lambda row: change_Event_based_on_mut_position(row['Event'], row['mutation_around_cut']), axis=1 )
    return df

def create_mastertable (project_dir, project_output_dir, parameters_dic):
    """Finds matching files 'Alleles_frequency_table_around' from CRISPResso2 output, and makes a joined table from them. Plus adds columns from the names."""
    now = datetime.now()
    date = now.strftime('%y%m%d')

    # Get list of tables, generated by CRISPResso2
    list_of_tables = find_files_paths_with_criteria(project_dir, 'Alleles_frequency_table_around', '.txt')
    write_into_report(f'The code will process the followng list of files: {list_of_tables}')
    list_of_edited_tables = []
    for filepath in list_of_tables:
        # get the sequence of the gRNA from the file name
        gRNA_name = split(basename(filepath))[-1]
        gRNA_name = gRNA_name.split("_")[-2]

        folder_name = split(dirname(filepath))[-1]

        sample_name = gRNA_name + '_' + folder_name[14:]
        # nuclease = filepath.split('/')[-2].split('_')[3]
        table = pd.read_csv(filepath, sep='\t')
        print(table)
        table = process_table(table, filepath[-24:-4], sample_name, parameters_dic)
        print(table.shape)
        # the name of the Allels_frequency_table_* must have sequence of the gRNA at the end of the name
        # it is there automatically if the name of the gRNA is not given into CRISPResso2 tool
        # example, where I added the sequence post analysis via search in totalcommander -> Alleles_frequency_table_around_gRNA10_GTAAAGTTTTCTGGTTAATT.txt
        table_folder_path = dirname(os.path.normpath(filepath))
        table_name = f'00_mastertable_{folder_name[14:]}.tsv'
        table.to_csv(Path(str(table_folder_path) + '/' + table_name), sep='\t', index=False)
        list_of_edited_tables.append(table_folder_path + '/' + table_name)
        write_into_report(f'{table_name} with shape {table.shape} was created, no filter was used.')

    write_into_report('Start to make mastertable from individual sample tables.')
    # Create a mastertable CSV file to write the merged data from individual tables
    mastertable_dir = Path(str(project_output_dir) + f'/00_mastertable_{date}.tsv')
    counter = 1
    no_of_tables = len(list_of_edited_tables)
    skip_header = False
    with open(mastertable_dir, 'w', newline='', encoding='utf-8') as merged_file:
        # Initialize the CSV writer object
        writer = csv.writer(merged_file)
        # Loop through each CSV file
        for file_name in list_of_edited_tables:
            # Open the CSV file and read the data
            with open(file_name, 'r', encoding='utf-8') as csv_file:
                reader = csv.reader(csv_file)
                if skip_header is True:
                    next(reader)     # skip the header row
                else:
                    skip_header = True     # keep the header row from the first file
                data = [row for row in reader]
            # Write the data to the merged CSV file
            writer.writerows(data)
            write_into_report(f'Table {counter} out of {no_of_tables} was appended to mastertable.')
            counter+=1
    write_into_report(f'00_mastertable_{date}.tsv was created.')
    return

# Add columns to the 01_mastertable
def add_columns_to_df(df, project_output_dir, parameters_dic, date):
    """Adds columns to the mastertable based on the given parameters and filters for #Reads column higher than minimum_Reads. The initial names of the source files from CRISPREsso2 must follow exact orders: Order_CasVariant_Gene_gRNAnum_Mutant_NoOfTargets_Overhangs_Replicate. Saves the updated mastertable as a TSV file."""
    # filter out #Reads<minimum_Reads
    df = df[df['#Reads'] >= int(parameters_dic['minimum_Reads'])].copy()
    df['order'] = df.apply(lambda row: row['sample'].split('_')[1], axis=1)
    df['Cas_variant'] = df.apply(lambda row: row['sample'].split('_')[2], axis=1)
    df['gene'] = df.apply(lambda row: row['sample'].split('_')[3], axis=1)
    df['gRNA_num'] = df.apply(lambda row: row['sample'].split('_')[4], axis=1)
    df['mutant'] = df.apply(lambda row: row['sample'].split('_')[5], axis=1)
    df['no_of_targets'] = df.apply(lambda row: row['sample'].split('_')[6], axis=1)
    df['overhang'] = df.apply(lambda row: row['sample'].split('_')[7], axis=1)
    df['replicate'] = df.apply(lambda row: row['sample'].split('_')[8], axis=1)
    for key in parameters_dic['add_columns']:
        df[key] = parameters_dic['add_columns'][key]
    df['name'] = df['sample'].astype(str) + '_' + df['gRNA_seq'].astype(str)

    df.sort_values(by='name', inplace=True)
    df.to_csv(Path(project_output_dir, '01_mastertable_' + date + '.tsv'), sep='\t', mode='w', header=True, index=None)
    write_into_report(f'01_mastertable_{date}.tsv with shape {df.shape} was created.\n')
    return df

def make_subtables_df(df, project_output_dir, subtables_param, date):
    """
    Creates subtables from the main dataframe based on the given parameters and saves them as TSV files.
    
    Parameters:
    - df (pd.DataFrame): The main dataframe to process.
    - project_output_dir (str or Path): The directory where the subtable TSV files will be saved.
    - subtables_param (list): A list of column names to create subtables for. The function will create combinations of subtables based on unique values in these columns.
    - date (str): The date to be included in the filenames of the output TSV files.
    
    Returns:
    - tables_dic (dict): A dictionary where keys are the file paths of the created subtables and values are the lists of parameter values used to create each subtable.
    """
    tables_dic = {}
    print(subtables_param)

    # Recursive function to handle any number of parameters
    def create_subtables(df, params, prefix):
        if len(params) == 0:
            # Convert all elements in prefix to strings before joining
            prefix_str = [str(item) for item in prefix]
            csv_dir = Path(str(project_output_dir) + f"/02_mastertable_{'_'.join(prefix_str)}_{str(date)}.tsv")
            df.to_csv(csv_dir, sep="\t", mode="w", header=True, index=None)
            write_into_report(f"Subtable 02_mastertable_{'_'.join(prefix_str)}_{date}.tsv with shape {df.shape} was created.")
            tables_dic[csv_dir] = prefix
        else:
            param = params[0]
            remaining_params = params[1:]
            for value in df[param].unique():
                df_filtered = df[df[param] == value]
                create_subtables(df_filtered, remaining_params, prefix + [value])

    # Start the recursive subtable creation
    create_subtables(df, subtables_param, [])

    print("Subtables were created and saved. Output is table_dic.")
    return tables_dic
