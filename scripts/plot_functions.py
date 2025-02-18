# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""File contains functions for plotting charts, used in the CRISPRlysis project."""

from datetime import datetime
from glob import glob
import os
from os.path import dirname, split, basename
from pathlib import Path
import shutil
import time
import sys

import matplotlib
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.patches import PathPatch
import matplotlib.transforms as transforms
import more_itertools as mit
import numpy as np
import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation

from general_functions import write_into_report

def save_stat(stat_info, name_to_save, project_output_dir):
    """
    Save information as statistic into output *_stat.txt file

    Parameters:
    stat_info (str): The statistical information to be saved.
    name_to_save (str): The base name for the output file. The function will append '_stat.txt' to this name.
    project_output_dir (str): The directory where the output file will be saved.
    
    Returns:
    no value, the *stat.txt file will be saved"""
    with open(Path(project_output_dir, f'{name_to_save}_stat.txt'), 'w', encoding='utf-8') as f:
        f.write(str(stat_info))
    write_into_report(f'Statistic was saved as {name_to_save}_stat.txt', destination_folder=project_output_dir)

def make_Ubuntu_mnt_path(path):
    """
    When using WSL in Windows, the paths from format r'D:\\path\\to\\file' are converted to /mnt/d/path/to/file
    
    Parameters:
        path (str): A Windows-style path, e.g., r'D:\\path\\to\\file'

    Returns:
        path (str): The converted WSL path, e.g., /mnt/d/path/to/file"""
    # Check if the letter of the disk drive is capitalised, if yes make it lower
    if path[0].isupper():
        path = path[0].lower() + path[1:]
    path = path.replace(':','').replace('\\', '/')
    FinalPath = '/mnt/' + path
    return FinalPath

def make_reverse_complement(DNA_seq):
    """
    Makes reverse complement of given DNA sequence.

    Parameters:
        seq (str): A,C,G,T,- sequence.

    Returns:
        reverse_complement (str): The reverse complement of the input DNA sequence.
    """
    # Valid characters in a DNA sequence
    valid_chars = {'A', 'C', 'G', 'T', '-'}
    # Check if all characters in the input sequence are valid
    if not set(DNA_seq).issubset(valid_chars):
        raise ValueError('Invalid characters found in DNA sequence. Only A, C, G, T, and - are allowed.')
    # Dictionary mapping each nucleotide to its complement
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
    # Generate the reverse complement by reversing the sequence and replacing each nucleotide with its complement
    reverse_complement = ''.join(complement.get(base, base) for base in reversed(DNA_seq))
    return reverse_complement

def find_element_in_list(element, list_element):
    """
    Finds the 1-based index of an element in a list.

    Parameters:
    - element: The element to search for in the list.
    - list_element: The list in which to search for the element.

    Returns:
    - The 1-based index of the element if found, otherwise None.
    """
    try:
        # Find the index of the element (0-based index)
        index_element = list_element.index(element)
        # Return the 1-based index
        return index_element + 1
    except ValueError:
        # Return None if the element is not found
        return None

def adjust_box_widths(g, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.

    Parameters:
    - g (seaborn.axisgrid.FacetGrid): The seaborn FacetGrid object containing the boxplot.
    - fac (float): The factor by which to adjust the width of the boxes.

    Requirements:
    - from matplotlib.patches import PathPatch
    - import numpy as np

    Returns:
    - None
    """
    # Import PathPatch if not already imported
    if 'PathPatch' not in sys.modules:
        import matplotlib.patches as PathPatch
    if 'numpy' not in sys.modules:
        import numpy as np

    # Iterating through Axes instances
    for ax in g.axes:
        # Iterating through axes artists
        for c in ax.get_children():
            # Searching for PathPatches (boxplot boxes)
            if isinstance(c, PathPatch):
                # Getting current width of box
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]  # All but the last vertex
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)
                # Setting new width of box
                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new
                # Setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
    return
