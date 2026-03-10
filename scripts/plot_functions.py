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
# TODO remove statannot once all works on statannotations bellow
# from statannot import add_stat_annotation
from statannotations.Annotator import Annotator

from .general_functions import write_into_report

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

def adjust_box_widths(g, fac, no_of_columns):
    """
    Adjust the widths of a seaborn-generated boxplot.

    Parameters:
    - g (seaborn.axisgrid.FacetGrid): The seaborn FacetGrid object containing the boxplot.
    - fac (float): The factor by which to adjust the width of the boxes.
    - no_of_columns (integer): The script did not work when there were only 2 categories to be plotted, this is workaround.

    Requirements:
    - from matplotlib.patches import PathPatch
    - import numpy as np

    Returns:
    - None
    """
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
                    if no_of_columns == 2:
                        xdata = l.get_xdata()
                        if np.isclose(np.min(xdata), xmin) and np.isclose(np.max(xdata), xmax):
                            l.set_xdata([xmin_new, xmax_new])
                    else:
                        if np.all(l.get_xdata() == [xmin, xmax]):
                            l.set_xdata([xmin_new, xmax_new])
                
    return

# TODO offset of datapoints does not work
def plot_sns_violin_and_print(data, x, y, number_of_columns, project_output_dir\
        , hue=None, hue_order=None, hue_only_for_data_points=True\
        , boxplot_distances_with=[0.15, 0.9]\
        , ylabel=None, xlabel=None, title=None, xtickt_rotation=0\
        , order=None\
        , integrated_properties_dic=False, color_dic=None\
        , hatches=False, hatches_list=None\
        , change_edge_color=False, edge_color_dic=None\
        , show_datapoints=True, jitter=True, stripplot_offset=-10, dodge=True\
        , show_hue_medians=True\
        , transparency = 0.8 \
        , clip_right=False\
        , yaxis_formater=1, name_to_save='Boxplot', data_source = None\
        , y_min=None, y_max=None\
        , x_min=None, x_max=None\
        , plot_vertical_line_middle=False, middle_line_position=1\
        , plot_legend=False\
        , add_stat=False, stat_pairs=None, stat_test='t-test_ind', stat_format='star', stat_correction=None\
        , text_bellow_chart=None, y_legend_coordinate=0\
        , output_format='png', dpi=200, show_plots=False):
    """Plots violinplot and saves it.
    data: DataFrame
    x: x to plot
    y: y to plot
    number_of_columns: int() to calculate the chart width
    project_output_dir: dir for output
    hue: hue or None
    hue_order: list or None, default None
    hue_only_for_data_points: True/False, default True
    boxplot_distances_with: [distance, width], default [0.15, 0.9]
    ylabel: ylabel, default None
    xlabel: xlabel, default None
    title: title, default None
    xtickt_rotation: int(), default 0
    order: order of x, default None
    integrated_properties_dic (bool):True/False use integrated dictionary for color/hatch/edge_color
    color_dic: dic{} if integrated_properties_dic==False, default None
    hatches (bool): True/False, if you want to make hatches, default False
    hatches_list (list, dict, None): is used if integrated_properties_dic== False, deault None
    change_edge_color (bool): True/False, default False
    edge_color_dic: dic{} if integrated_properties_dic==False, default None
    show_datapoints (bool): True/False, default True
    jitter: True/False/int() for stripplot, default True
    stripplot_offset (int): to move stripplot to the left or right, default -10
    dodge: True/False separated stripplot hues - needed when hues are used, default True
    show_hue_medians: True/False for stripplot, default True
    transparency = int() 0-1, default 0.8
    clip_right (bool): True/False removes right part of the violine plot, default False
    yaxis_formater: None or usually 1 of 100, to make yaxis in %, default 1
    name_to_save: hae to contain full path with name of the boxplot (date and format is added automatialy)
    data_source: like 'Přibylová et al., (2021)'
    y_min: int(), default None
    y_max: int(), default None
    x_min: int(), default None
    x_max: int(), default None
    plot_vertical_line_middle: True/False, default False
    middle_line_position: int(), default 1
    plot_legend: True/False, default False
    add_stat: add t-test_ind statistic annotation, default False
    stat_pairs: [(('pairs'), ('to compare')), (('pairs'), ('to compare')), ...]
    stat_test: t-test_ind, t-test_welch, t-test_paired, anova, pearsonr, Mann-Whitney, Wilcoxon, Kruskal, Friedman, Spermanr, default t-test_ind
    stat_correction: bonferroni (the most strict, for 1-3 comparisions), holm (2-20 comparisions, good one), fdr_bh (omice/multi-group experiments, more powerfull than holm), sidak (as bonferroni but les onservative), none, default None
    stat_format: star, simple, full - format of statistic, default star
    text_bellow_chart: str() to comment chart bellow, default None
    y_legend_coordinate: int() coordinates of figure are 0-1, default 0
    output_format: png, svg,..., default 'png'
    dpi: int(), default 200
    report_dir: where should be the report, default os.getcwd()
    show_plots: True/False, default False"""
    # Define the custom exception class within the function
    class PlotSnsViolinPlotError(Exception):
        """Rise a function specific error."""
    
    # Check if the integrated_properties_dic is a dictionary, if not use the default dictionary
    if isinstance(integrated_properties_dic, dict) is True:
        print('The integrated_properties_dic is a dictionary, and it is going to be used.')
        properties_dic = integrated_properties_dic
        print(properties_dic)
    elif integrated_properties_dic is True:
        # TODO warning, markers does not work, yet...
        properties_dic = {'stagg2nt':{'color':'#E890B8','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'stagg2nt_likely':{'color':'#E890B8','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'stagg1nt':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'stagg1nt_likely':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'////', 'marker':'o'}
            ,'blunt':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'blunt_likely':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'cannot_be_determined':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'can_not_be_determined':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'can not be determined':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'cannot_be_determined_ambiguous':{'color':'white','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'cannot_be_determined_NO_ambi':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'distal_stagg2nt_likely':{'color':'#E890B8','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'distal_stagg2nt':{'color':'#E890B8','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'distal_stagg1nt_likely':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'distal_stagg1nt':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'distal_blunt_likely':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'distal_blunt':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_blunt':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_blunt_likely':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'proximal_stagg1nt':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_stagg1nt_likely':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'proximal_stagg2nt':{'color':'#E890B8','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_stagg2nt_likely':{'color':'#E890B8','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'gap_out_of_cut':{'color':'gray','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'IN':{'color':'#091f92','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'INs':{'color':'#091f92','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'DEL':{'color':'#ba8759','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'DELs':{'color':'#ba8759','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'SUB':{'color':'#ffcc33','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'OTHERS':{'color':'#808080','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'no_mutation':{'color':'#111111','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'U':{'color':'#f0d271','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'B':{'color':'#808080','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'D':{'color':'#f5ad62','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'X':{'color':'#d3d3d3','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'X_gap_out_of_cut':{'color':'#d3d3d3','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'upstream':{'color':'#f0d271','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'upstream_ambiguous':{'color':'#f0d271','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'upstream_NO_ambi':{'color':'#f0d271','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'upstream_both':{'color':'#C4C4C4','edge_color':'#f0d271', 'hatch':'\\\\\\', 'marker':'o'}
            ,'upstream_both_downstream':{'color':'#979797','edge_color':'black', 'hatch':'xx', 'marker':'o'}
            ,'both':{'color':'#808080','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'both_ambiguous':{'color':'#808080','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'both_NO_ambi':{'color':'#808080','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'both_downstream':{'color':'#C4C4C4','edge_color':'#f5ad62', 'hatch':'////', 'marker':'o'}
            ,'downstream':{'color':'#f5ad62','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'downstream_ambiguous':{'color':'#f5ad62','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'downstream_NO_ambi':{'color':'#f5ad62','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'deletion_spans_analysed_window':{'color':'#a9a9a9','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'deletion spans analysed window':{'color':'#a9a9a9','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'PAM_distal_nt_in_MH':{'color':'#52D8FF','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'PAM_proximal_nt_in_MH':{'color':'#FADA52','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'no_MH_likely':{'color':'#A9A9A9','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'no_MH':{'color':'#D3D3D3','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'A':{'color':'#70c936','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'C':{'color':'#ff0000','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'G':{'color':'#ffd42a','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'T':{'color':'#9955ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,True:{'color':'#D3D3D3','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,False:{'color':'#B8BEDB','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'True':{'color':'#D3D3D3','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'False':{'color':'#B8BEDB','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'unmethylated':{'color':'#000000','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'CG_only':{'color':'#4D4D4D','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'heterochromatin':{'color':'#999999','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'RdDM':{'color':'#fafafa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'Homo sapiens - HCT116 - 4': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Homo sapiens - HCT116 - 8': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Homo sapiens - HCT116 - 16': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': '^'}
            ,'Homo sapiens - HCT116 - 24': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': '<'}
            ,'Homo sapiens - HCT116 - 48': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo sapiens - HEK293 - 4': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Homo sapiens - HEK293 - 8': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Homo sapiens - HEK293 - 16': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': '^'}
            ,'Homo sapiens - HEK293 - 24': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': '<'}
            ,'Homo sapiens - HEK293 - 48': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo sapiens - K562 - 4': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'*'}
            ,'Homo sapiens - K562 - 8': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'Homo sapiens - K562 - 16': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'^'}
            ,'Homo sapiens - K562 - 24': {'color': "#A55D84", 'edge_color':'gray', 'hatch':'', 'marker':'<'}
            ,'Homo sapiens - K562 - 48': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'s'}
            ,'Homo sapiens - HAP1 - 48to72': {'color': '#7D3E60', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo sapiens - Jurkat - 48to72': {'color': '#810F72', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo sapiens - Hela - 72': {'color': '#520047', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Arabidopsis thaliana - xxx - 504': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Arabidopsis thaliana - root culture - 504': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Arabidopsis thaliana - xxx - unknown': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'A. thaliana - embryo - unknown': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Nicotiana benthamiana - xxx - 120': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'v'}
            ,'Petunia hybrida - xxx - 24': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '>'}
            ,'Silene latifoliaa - xxx - xx': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '<'}
            ,'Solanum lycopersicum - xxx - 24': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'plant': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'animal': {'color': '#cc79a6', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            }
    else:
        print('The properties_dic is not a dictionary, and also the integrated dictionary is not going to be used.')

    # Check if there are data in the data[x], if not skip making the plot.
    if data[x].notnull().any() is False:
        write_into_report(f'data["x"] are empty, plotting of the plot {name_to_save} was terminated.', destination_folder=project_output_dir)
        raise PlotSnsViolinPlotError(f'data["x"] are empty, plotting of the plot {name_to_save} was terminated.')

    now = datetime.now()
    date = now.strftime('%y%m%d')
    # Change matplotlib style and define figure size
    sns.set_theme(style='whitegrid')
    fig, ax = plt.subplots(figsize=(number_of_columns, 4))

    # This script can not make hatches when usin hue for violinplot
    if hue is not None:
        if hue_only_for_data_points is False:
            hatches=False
            print('This script can not make hatches when usin hue for violinplot. Hatches was set to False.')
    # Make color,edge and hatch dic from nested dictionary
    if integrated_properties_dic is not False:
        color_dic = {}
        edge_color_list = []
        hatches_list = []
        test = data.groupby(x).value.mean().to_dict()
        for item in list(data[x].unique()):
            if item in properties_dic.keys():
                color_dic[item] = properties_dic[item]['color']
                if test[item] != 0:
                    edge_color_list.append(properties_dic[item]['edge_color'])
                    hatches_list.append(properties_dic[item]['hatch'])
            if item not in properties_dic.keys():
                print('The ' + str(item) + ' is not in properties_dic of plot_sns_violin_and_print function. The dic contains:' + str(color_dic))
        if hue is not None:
            for item in list(data[hue].unique()):
                if item in properties_dic.keys():
                    color_dic[item] = properties_dic[item]['color']
                if item not in properties_dic.keys():
                    print('The ' + str(item) + ' is not in properties_dic of plot_sns_violin_and_print function.')
    else:
        color_dic = color_dic
        hatches_list = hatches_list
        edge_color_dic = edge_color_dic

    # to make hue color only for stripplot and not for violin
    if hue_only_for_data_points is False:
        hue_violin = hue
    else:
        hue_violin=None

    # Set withts for plotted boxplots
    width_list = boxplot_distances_with
    # TODO solve in the main CRISPRlysis.py code and set up the new parameter boxplot_distances_with likely as [0.8, 0.3]
    # if hue_only_for_data_points is False:
    #     if hue is not None:
    #         if len(data[hue].unique()) == 4:
    #             width_list = [0.8, 0.3]
    #         else:
    #             width_list = [0.8, 0.4]
    #             print('Problem with boxplot wisht could be in settings width_list. Change values in the with_list = [0.8, 0.4]')

    # Plot Violinplot in order (True/False) with stripplot (True/False)
    if order is None:
        # Make Violinplot
        # Change color_dic keys to 'txt' as the violinplot bellow makes from True -> 'True'
        # and the plot does not work without the conversion -> ideas how to solve it?
        color_dic_txt = {str(key): value for key, value in color_dic.items()}
        sns.violinplot(data=data, x=x, y=y, hue=hue_violin, hue_order=hue_order, density_norm='width', palette=color_dic_txt, cut=0, inner=None, zorder=0)
        ## Set hatches if wanted
        if hatches is True:
            # changing hatches for if hatches_list is not None to be able to use predefined list
            if hatches_list is not None:
                if isinstance(hatches_list, list):
                    hatches_list = hatches_list
                elif isinstance(hatches_list, dict):
                    hatches_list_new = []
                    for key in order:
                        if key in hatches_list:
                            hatches_list_new.append(hatches_list[key])
                    hatches_list = hatches_list_new
                else:
                    raise PlotSnsViolinPlotError('Problem with hatches_list variable, it is not list or dict. The hatches is set to True.')
            else:
                hatches_list = []
                for item in order:
                    if test[item] != 0:
                        hatches_list.append(properties_dic[item]['hatch'])
            ihatch = iter(hatches_list)
            hatch = [i.set_hatch(next(ihatch)) for i in ax.get_children() if isinstance(i, matplotlib.collections.PolyCollection)]
        # Clip the right half of each violin.
        if clip_right is True:
            for item in ax.collections:
                x0, y0, width, height = item.get_paths()[0].get_extents().bounds
                item.set_clip_path(plt.Rectangle((x0-0.1, y0), (width/2)+0.1, height,
                                angle=0, transform=ax.transData))
        # Make stripplot to show individual data
        if show_datapoints is True:
            # move all of the datapoints a bit to the left or right
            def offset(p):
                return transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
            trans = plt.gca().transData
            # Note: be careful when using dodge as it can mess up the coordinates when combined with jitter
            stripplot = sns.stripplot(x=x, y=y, hue=hue, hue_order=hue_order, data=data, jitter=jitter\
                , order=None, dodge=dodge\
                , palette=color_dic\
                , size=3, linewidth=0.5\
                , transform=trans+offset(stripplot_offset)\
                , zorder=20)
            # Add median for stripplot
            if show_hue_medians is True:
                sns.boxplot(showmeans=False,
                    meanline=False,
                    meanprops={'marker':'x'},
                    medianprops={'visible': True},
                    whiskerprops={'visible': False},
                    boxprops={'visible':False},
                    zorder=20,
                    x=x, y=y, hue=hue, hue_order=hue_order, data=data,
                    palette=color_dic,
                    showfliers=False,
                    showbox=True, # must be True for Change color of meanprops later
                    showcaps=False,
                    ax=stripplot)
                # Change color of the mean according to the used palette in color_dic
                
                # box_patches = [patch for patch in ax.patches if type(patch) == matplotlib.patches.PathPatch]
                # if len(box_patches) == 0:  # in matplotlib older than 3.5, the boxes are stored in ax.artists
                #     box_patches = ax.artists
                # num_patches = len(box_patches)
                # lines_per_boxplot = len(ax.lines) // num_patches
                # for i, patch in enumerate(box_patches):
                #     # Set the linecolor on the patch to the facecolor, and set the facecolor to None
                #     col = patch.get_facecolor()
                #     patch.set_edgecolor(col)
                #     patch.set_facecolor('None')
                #     # Each box has asociated Line2D objects (to make the whiskers, fliers, etc.)
                #     # Loop over them here, and use the same color as above
                #     for line in ax.lines[i * lines_per_boxplot: (i + 1) * lines_per_boxplot]:
                #         line.set_color(col)
                #         line.set_mfc(col)  # facecolor of fliers
                #         line.set_mec(col)  # edgecolor of fliers
        # Make inner boxplot with mean as 'x'
        sns.boxplot(x=x, y=y, hue=hue_violin, hue_order=hue_order, data=data, width=width_list[0], showfliers=False, showmeans=True, meanprops=dict(marker='x', markeredgecolor='black', markerfacecolor='black', markersize=5, linewidth=1.0, zorder=14), boxprops=dict(facecolor='gray', edgecolor='black', linewidth=1.5, alpha=1, zorder=10), whiskerprops=dict(color='black', linewidth=1.5, zorder=11), capprops=dict(color='black', linewidth=1.5, zorder=12), medianprops=dict(color='black', linewidth=1.5, zorder=13))
    else:
        # Make Violinplot
        sns.violinplot(x=x, y=y, hue=hue_violin, hue_order=hue_order, data=data, density_norm='width', palette=color_dic, cut=0, inner=None, order=order, zorder=1)
        ## Set hatches if wanted
        if hatches is True:
            # changing hatches for if hatches_list is not None to be able to use predefined list
            if hatches_list is not None:
                if isinstance(hatches_list, list):
                    hatches_list = hatches_list
                elif isinstance(hatches_list, dict):
                    hatches_list_new = []
                    for key in order:
                        if key in hatches_list:
                            hatches_list_new.append(hatches_list[key])
                    hatches_list = hatches_list_new
                else:
                    raise PlotSnsViolinPlotError('Problem with hatches_list variable, it is not list or dict. The hatches is set to True.')
            else:
                hatches_list = []
                for item in order:
                    if test[item] != 0:
                        hatches_list.append(properties_dic[item]['hatch'])
            ihatch = iter(hatches_list)
            hatch = [i.set_hatch(next(ihatch)) for i in ax.get_children() if isinstance(i, matplotlib.collections.PolyCollection)]

        # Clip the right half of each violin.
        if clip_right is True:
            for item in ax.collections:
                x0, y0, width, height = item.get_paths()[0].get_extents().bounds
                item.set_clip_path(plt.Rectangle((x0-0.1, y0), (width/2)+0.1, height,
                                angle=0, transform=ax.transData))
        # Make stripplot to show individual data
        if show_datapoints is True:
            # move all of the datapoints a bit to the left or right
            def offset(p):
                return transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
            trans = plt.gca().transData
            stripplot = sns.stripplot(x=x, y=y, hue=hue, hue_order=hue_order, data=data, jitter=jitter\
                , order=order, dodge=dodge\
                , palette=color_dic\
                , size=3, linewidth=0.5\
                , transform=trans+offset(stripplot_offset)\
                , zorder=20)
            # Add mean for stripplot
            if show_hue_medians is True:
                sns.boxplot(showmeans=False,
                    meanline=False,
                    # meanprops={'marker':'x'},
                    medianprops={'visible': True},
                    whiskerprops={'visible': False},
                    boxprops={'visible':False},
                    zorder=20,
                    x=x, y=y, hue=hue, hue_order=hue_order, data=data,
                    palette=color_dic, order=order,
                    showfliers=False,
                    showbox=True, # must be True for Change color of meanprops later
                    showcaps=False,
                    ax=stripplot)
                # Change color of the mean according to the used palette
                box_patches = [patch for patch in ax.patches if type(patch) == matplotlib.patches.PathPatch]
                if len(box_patches) == 0:  # in matplotlib older than 3.5, the boxes are stored in ax.artists
                    box_patches = ax.artists
                num_patches = len(box_patches)
                lines_per_boxplot = len(ax.lines) // num_patches
                for i, patch in enumerate(box_patches):
                    # Set the linecolor on the patch to the facecolor, and set the facecolor to None
                    col = patch.get_facecolor()
                    patch.set_edgecolor(col)
                    patch.set_facecolor('None')
                    # Each box has asociated Line2D objects (to make the whiskers, fliers, etc.)
                    # Loop over them here, and use the same color as above
                    for line in ax.lines[i * lines_per_boxplot: (i + 1) * lines_per_boxplot]:
                        line.set_color(col)
                        line.set_mfc(col)  # facecolor of fliers
                        line.set_mec(col)  # edgecolor of fliers
        # Make inner boxplot with mean as 'x'
        sns.boxplot(x=x, y=y, hue=hue_violin, hue_order=hue_order, data=data\
            , width=width_list[0]\
            , order=order, showfliers=False, showmeans=True\
            , meanprops=dict(marker='x', markeredgecolor='black', markerfacecolor='black'\
                , markersize=5, linewidth=1.0, zorder=14)\
            , boxprops=dict(facecolor='gray', edgecolor='black', linewidth=1.5, alpha=1, zorder=10)\
            , whiskerprops=dict(color='black', linewidth=1.5, zorder=11)\
            , capprops=dict(color='black', linewidth=1.5, zorder=12)\
            , medianprops=dict(color='black', linewidth=1.5, zorder=13))

    # Change edge color if wanted
    if change_edge_color is True:
        if order is not None:
            edge_color_list = []
            for item in order:
                if test[item] != 0:
                    edge_color_list.append(properties_dic[item]['edge_color'])
        iedge_color = iter(edge_color_list)
        for i in range(0,len(edge_color_list)):
            ax.collections[i].set_edgecolor(next(iedge_color))

    # Statistic
    if add_stat is True:
        # add_stat_annotation(ax, data=data, x=x, y=y\
        #     , box_pairs=stat_pairs,\
        #     test=stat_test, text_format='star', loc='outside', verbose=1, fontsize='large')
        annot = Annotator(ax, pairs=stat_pairs, data=data, x=x, y=y)
        if stat_correction in [None, "none", "no_correction"]:
            annot.configure(
                test=stat_test,
                text_format=stat_format
            )
        else:
            annot.configure(
                test=stat_test,
                text_format=stat_format,
                pvalue_correction_method=stat_correction
            )

        # After configuration:
        test_used = annot.test
        annot.apply_and_annotate()
        ax.text( 1.02, 0.00, f"Test: {test_used}\nCorrection: {stat_correction}\nα = 0.05", transform=ax.transAxes, ha='left', va='bottom', fontsize=10 )

    # Plot vertical middle line
    if plot_vertical_line_middle is True:
        plt.axvline(x=middle_line_position-0.5, color='black', ls='--', label='blunt cut site')

    # Legend
    if plot_legend is True:
        if show_datapoints is True:
            if show_hue_medians is True:
                # Plot legend only from the stripplot and not boxplot (they are identical)
                handles, labels = ax.get_legend_handles_labels()
                plt.legend(handles[int(len(handles)/2):int(len(handles))],\
                    labels[int(len(labels)/2):int(len(labels))],\
                    loc='center left', bbox_to_anchor=(1, 0.5))    
            else:
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            if show_hue_medians is True:
                # Plot legend only from the stripplot and not boxplot (they are identical)
                handles, labels = ax.get_legend_handles_labels()
                plt.legend(handles[0:int(len(handles)/2)],\
                    labels[int(len(labels)/2):int(len(labels))],\
                    loc='center left', bbox_to_anchor=(1, 0.5))    
            else:
                # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                # Plot legend only from the stripplot and not boxplot (they are identical)
                handles, labels = ax.get_legend_handles_labels()
                plt.legend(handles[0:int(len(handles)/2)],\
                    labels[int(len(labels)/2):int(len(labels))],\
                    loc='center left', bbox_to_anchor=(1, 0.5))    
    # Plot labels, set lims, etc.
    plt.title(title + '\n(' + date + ')\n' + data_source)
    plt.ylabel(ylabel)
    if xlabel is not None:
        plt.xlabel(xlabel)
    plt.ylim(y_min, y_max)
    plt.xlim(x_min, x_max)
    plt.xticks(rotation=xtickt_rotation, ha='left', rotation_mode='anchor')

    # Add transparency to colors
    plt.setp(ax.collections, alpha=transparency)

    if yaxis_formater is not None:
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(yaxis_formater, decimals=0))
    # makes gaps between boxes when hue is used
    # If order is None, get order from data
    if order is None:
        order = list(data[x].unique())
    adjust_box_widths(fig, width_list[1], len(order))
    # Add text below chart
    if text_bellow_chart is not None:
        props = {'backgroundcolor':'lightgray', 'size':12}  # bbox features
        at = AnchoredText(text_bellow_chart, frameon=True, loc='upper left'\
            , bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=fig.transFigure, borderpad=0., prop=props)
        at.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
        ax.add_artist(at)

    plt.savefig(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.' + output_format), dpi=dpi, bbox_inches='tight')
    write_into_report(f'Violin plot {name_to_save} was created.', destination_folder=project_output_dir)
    if show_plots is True:
        return plt.show()
    return plt.close()

def check_if_all_events_are_present_after_pivot_reset_index(df, list_of_events, cut_maximum, cut_minimum, return_only_predef_list=False):
    """
    Parameters:
    df (DataFrame): DataFrame
    list_of_events (str or list): INDELSUBOTHERS or BluntStag or ProximalDistalBars or [individual, list]
    cut_maximum (int): is the most 5' position where the Cas can cut
    cut_minimum (int): is the most 3' postition where the Cas can cut
    return_only_predef_list (bool): default False
    
    Returns:
    Edited DataFrame
    list_of_events
    """
    # TODO solve list of events fo bluntstag
    # Make list based on the cut_maximum and cut_optimum
    default_events = ['C_0', 'C_1', 'O_0', 'O_1', 'X_0', 'X_1']
    if list_of_events == 'BluntStag':
        list_of_events = []
        for i in range(cut_maximum, cut_minimum+1):
            list_of_events.append(f'{i}_0')
            list_of_events.append(f'{i}_1')
        list_of_events.extend(default_events)
        # concert int() to str()
        # for i in range(len(list_of_events)):
        #     list_of_events[i]=str(list_of_events[i])
        # TODO delete onec it works
        # list_of_events = ['stagg2nt', 'stagg2nt_likely', 'stagg1nt', 'stagg1nt_likely', 'blunt_likely', 'blunt', 'cannot_be_determined']
    if list_of_events == 'INDELSUBOTHERS':
        list_of_events = ['IN', 'DEL', 'SUB', 'OTHERS']
    # TODO solve list of events for prosimal distal
    if list_of_events == 'ProximalDistalBars':
        list_of_events = ['distal_blunt', 'distal_blunt_likely', 'distal_stagg1nt', 'distal_stagg1nt_likely', 'distal_stagg2nt_likely', 'proximal_blunt', 'proximal_blunt_likely', 'proximal_stagg1nt', 'proximal_stagg2nt']
    if return_only_predef_list is True:
        return list_of_events
    # Check if the list is empty. If yes, return unchanged dataframe
    if not list_of_events:
        return df
    for item in list_of_events:
        if item not in df.columns.tolist():
            df[item] = 0
    return df, list_of_events

def plot_general_stacked_horizontal_barchart(df, x, y, values, project_output_dir\
    , check_columns=False, cut_maximum=-2, cut_minimum=0, check_columns_list=[]\
    , data_source=None, title=None, name_to_save="stacked_horizontal_barchart"\
    , output_format='png', dpi=200, xlabel=None, sns_palette=True, palette='hls'\
    , hatches_dic=None, x_as_perc=False, xaxis_formater=1, decimals=1, xlim=None\
    , text_bellow_chart=None, y_legend_coordinate='Default'\
    , y_order=None, bars_order=None, show_plots=False):
    """Plots horizontal stacked histogram. 
    df: DataFrame
    x: single column
    y: single column
    values: single column
    project_output_dir: path to the output folder
    check_columns: check if the data have all neede columns to plot the plot, default False
    cut_maximum: int(), default -2
    cut_minimum: int(), default 0
    check_columns_list: 'INDELSUBOTHERS' or 'BluntStag' or 'ProximalDistalBars' or [individual, list], default []
    data_source: usually name of the study as 'Přibylová et al., (2022)', default None
    title: Any, default None
    name_to_save: str(), default 'stacked_horizontal_barchart'
    output_format: png, svg..., default 'png'
    dpi: any, default 200
    xlabel: Anly, default None
    sns_palette: True or False, if True set parameter palette from sns palettes, if False enter your own dictionary, default True
    palette: choose from sns.color_palette, default 'hls'
    hatches_dic: {'key':'value'}, default None
    x_as_perc: True makes from xaxis axes with %, anything else is False, default=False
    xaxis_formater: int(), multiplies x axis
    decimals: if X_as_perc is True, set number of decimals, default 1
    xlim: (min, max), default None
    text_bellow_chart: str() to comment chart bellow, default None
    y_legend_coordinate: int() coordinates of figure are 0-1, default 'Default' calculates y_legend_coordinate=(-(1/y_len)*3)
    y_order: paste list of Ys with desired order, default None
    bars_order: paste list of bars with desired order, default None
    show_plots: True/False, default False"""

    # Define the custom exception class within the function
    # class PlotGeneralStackedHorizonError(Exception):
    #     """Rise a function specific error."""
    # Check if there are data in the data[x], if not skip making the plot.
    if df[x].notnull().any() is False:
        write_into_report(f'df["x"] are empty, plotting of the plot {name_to_save} was terminated.', destination_folder=project_output_dir)
        return
        # raise PlotGeneralStackedHorizonError(f'df["x"] are empty, plotting of the plot {name_to_save} was terminated.')

    now = datetime.now()
    date = now.strftime('%y%m%d')

    y_len = len(df[y].unique())
    # skip the plot if there are more than 20 unique ys
    if y_len > 40:
        write_into_report(f'There are more than 40 unique values in df["{y}"], plotting of the plot {name_to_save} was terminated.', destination_folder=project_output_dir)
        return
        # raise PlotGeneralStackedHorizonError(f'There are more than 20 unique values in df["{y}"], plotting of the plot {name_to_save} was terminated.')

    # Try to make the pivot table, if not possible, skip the plot
    try:
        pivot = pd.pivot_table(df, index=y, columns=x, values=values)
    except Exception as e:
        write_into_report(f'Pivot table could not be created for {name_to_save} plot. Error: {e}', destination_folder=project_output_dir)
        return
    # make dic with updated values -> value + sum from previous values nd from it df again
    # pivot = pd.pivot_table(df, index=y, columns=x, values=values)
    if bars_order != 0:
        pivot = pivot.reindex(bars_order, axis=1)
    bars = pivot.columns.tolist()
    no_of_bars = len(bars)
    df_pivot= pivot.reset_index()
    if check_columns is not False:
        df_pivot = check_if_all_events_are_present_after_pivot_reset_index(df_pivot, check_columns_list, cut_maximum, cut_minimum)
    y = df[y].unique()
    i = 0
    to_plot_dic = {}
    for row in y:
        list_of_row_values = []
        value=0
        for column in range(len(df_pivot.iloc[i,1:])):
            value=df_pivot.iloc[i,column+1] + value
            list_of_row_values.append(value)
        to_plot_dic[row] = list_of_row_values
        i+=1
    df = pd.DataFrame.from_dict(to_plot_dic, orient='index', columns=bars)

    # change order of rows if wanted
    if y_order is not None:
        df = df.rename_axis('y').reset_index()
        df.y = df.y.astype('category')
        df.y = df.y.cat.set_categories(y_order)
        df = df.sort_values(['y'])
        df = df.set_index('y')

    # plot
    sns.set_theme(style='whitegrid')                      # makes gid behind the plot
    f, ax = plt.subplots(figsize=(9,y_len/4))       # set the figure size

    # make color palette from sns.color_palette or own dictionary
    if sns_palette is True:
        pal =sns.color_palette(palette, no_of_bars)
    else:
        # if some of the column does not exist as a key in the dictionary -> use sns palette
        pal = []
        for column in bars:
            if column in palette:
                pal.append(palette[column])
            else:
                pal = sns.color_palette('hls', no_of_bars)
                print('Key ' + str(column) + ' is not present in color dictionary, thus default hls palette was used.')
                break

    if hatches_dic is not None:
        for column in range(no_of_bars,0,-1):              # plot in reverse
            sns.barplot(x=df.iloc[:,column-1], y=df.index, color=pal[column-1], edgecolor='#000000', \
                    alpha=1, errorbar=None, linewidth=0.2, hatch=hatches_dic[list(df.columns)[column-1]])
    else:
        for column in range(no_of_bars,0,-1):              # plot in reverse
            sns.barplot(x=df.iloc[:,column-1], y=df.index, color=pal[column-1], edgecolor='#000000', \
                    alpha=1, errorbar=None, linewidth=0.2)        
    # get the current ticks and use it to set the new ticks by sampling every other tick
    # ax.set_xticks(ax.get_xticks()[::2])

    plt.title(f'{title} ({date}) \n {data_source}')
    plt.xlabel(xlabel, size=15)
    plt.ylabel('')
    plt.yticks(fontname='monospace')
    bars = np.flip(bars)                            # change the bars order as they were plotted 
    plt.legend(labels=bars, loc='center left', bbox_to_anchor=(1.1, 0.5))
    if xlim is not None:
        plt.xlim(xlim)
    if x_as_perc is True:
        ax.xaxis.set_major_formatter(mtick.PercentFormatter(xaxis_formater, decimals=decimals))
    # Add text below chart
    if text_bellow_chart is not None:
        if y_legend_coordinate == 'Default':
            if y_len == 0:
                y_len = 3
            y_legend_coordinate=-(1/y_len)*3
        props = {'backgroundcolor':'lightgray', 'size':12}  # bbox features
        at = AnchoredText(text_bellow_chart, frameon=True, loc='upper left', bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=f.transFigure, borderpad=0., prop=props)
        at.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
        ax.add_artist(at)

    plt.savefig(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.' + output_format), dpi=dpi, bbox_inches='tight')
    write_into_report(f'Stacked horizontal barchart {name_to_save} was created.', destination_folder=project_output_dir)
    if show_plots is True:
        return plt.show()
    return plt.close()

# TODO plot_barplot when fill_empty_x=True changes the original data which are
# imput for the plot, solve why
def plot_barplot(data, x, y, project_output_dir, name_to_save, data_source, hue=None, color=None\
                , title='Bar plot', x_lim=None, y_lim=None, xlabel=None\
                , yaxis_formater=100, decimals=1, xtickt_rotation=0, text_bellow_chart=None\
                , errorbar=None, fill_empty_x=False, output_format='png', dpi=200, show_plots=False):
    ''' Return plot
    data: source DataFrame
    project_output_dir: path to the output folder
    name_to_save: name for the figure
    data_source: usually name of the study as 'Přibylová et al., (2022)'
    hue: default None
    color: color of bars as '#091f92', default None
    title: str(), default 'Bar plot'
    x_lim: (x_min, x_max), default None
    y_lim: (y_min, y_max), default None
    xlabel = str(), default None 
    yaxis_formater: to format y axis as %
    decimals: if yaxis_formater is set, default 1 
    xtickt_rotation: int(), default 0
    text_bellow_chart: str() to comment chart bellow, default None
    ci: sd, default None
    fill_empty_x: True/False -> True fills missing valuese on X axis when they are int(), to make it continuous, deault False
    output_format: png, svg...
    dpi: any, usually 200
    show_plots: True/False, default False'''
    # Define the custom exception class within the function
    class PlotBarplotError(Exception):
        """Rise a function specific error."""
    # Check if there are data in the data[x], if not skip making the plot.
    if data[x].notnull().any() == False:
        write_into_report(f'data["x"] are empty, plotting of the plot {name_to_save} was terminated.', destination_folder=project_output_dir)
        return
        # raise PlotBarplotError(f'data["x"] are empty, plotting of the plot {name_to_save} was terminated.')
    
    now = datetime.now()
    date = now.strftime('%y%m%d')
    sns.set_theme(style='whitegrid')

    # if you want to have continuous x-axis like 1,2,3,4,5,6,... -. works only for numbers
    if fill_empty_x is True:
        # Generate a range of continuous x values
        min_x = data[x].min()
        max_x = data[x].max()
        all_x_values = pd.DataFrame({x: range(min_x, max_x + 1)})

        # Merge DataFrames to fill in missing rows with 0 in the y column
        data_filled = pd.merge(all_x_values, data, on=x, how='left').fillna(0)
    else:
        data_filled = data

    f,ax = plt.subplots()
    sns.barplot(x=x, y=y, data=data_filled, hue=hue, capsize=0.2, errorbar='sd', alpha=0.8, color=color)
    plt.title(title + ' (' + date + ') \n' + data_source)
    plt.xlabel(xlabel)
    plt.ylabel('Average of %Reads')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    # for container in ax.containers:
    #     ax.bar_label(container)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(yaxis_formater, decimals=decimals))
    if xtickt_rotation != 0:
        plt.xticks(rotation=xtickt_rotation, ha='left', rotation_mode='anchor')
    plt.text(0, -0.02, text_bellow_chart, fontsize=14, transform=plt.gcf().transFigure)
    plt.savefig(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.' + output_format), dpi=dpi, bbox_inches='tight')
    write_into_report('Barplot {0} was created.'.format(name_to_save))
    if show_plots is True:
        return plt.show()
    return plt.close()

# TODO shift the chart when no positive values are present
def plot_general_stacked_vertical_barchart(df, x, y, values, project_output_dir, order_y=None\
    , data_source=None, title=None, name_to_save='stacked_vertical_barchart'\
    , output_format='png', dpi=200, xlabel=None, ylabel=None, sns_palette=True\
    , hatches_dic=None, palette='hls', y_as_perc=False, yaxis_formater=1, xlim=None\
    , ylim=None, plot_vertical_line_middle=False\
    , text_bellow_chart=None, y_legend_coordinate=0, fill_empty_x=False, show_plots=False):
    """Plots horizontal stacked histogram. 
    df: DataFrame
    x: single column
    y: single column
    values: single column
    project_output_dir: path to the output folder
    check_columns: check if the data have all neede columns to plot the plot, default False
    check_columns_list: 'INDELSUBOTHERS' or 'BluntStag' or 'ProximalDistalBars' or [individual, list], default []
    order_y: list['order', 'of', 'y', 'categories'], default None
    data_source: usually name of the study as 'Přibylová et al., (2022)', default None
    title: Any, default None
    name_to_save: str(), default 'stacked_vertical_barchart'
    output_format: png, svg..., default 'png'
    dpi: any, default 200
    xlabel: Any, default None
    ylabel: Any, default None
    sns_palette: True or False, if True set parameter palette from sns palettes, if False enter your own dictionary, default True
    palette: if sns_palette is True -> choose from sns.color_palette, if False enter your dictionary, default 'hls'
    hatches_dic: {'key':'value'}, default None
    y_as_perc: True makes from xaxis axes with %, anything else is False, default=False
    xaxis_formater: int(), multiplies x axis
    xlim: (min, max), default None
    ylim: (min, max), defaule None
    text_bellow_chart: str() to comment chart bellow, default None
    plot_vertical_line_middle: True/False, default False, makes dashed line between -1 and 1
    fill_empty_x: True/False -> True fills missing valuese on X axis when they are int(), to make it continuous, deault False
    show_plots: True/False, default False"""
    # Check if there are data in the data[x], if not skip making the plot.
    if df[x].notnull().any() == False:
        print(f'df["x"] are empty, plotting of the plot {name_to_save} was terminated.')
        write_into_report(f'df["x"] are empty, plotting of the plot {name_to_save} was terminated.')
        return
    # Get actual date and time for plot annotation
    now = datetime.now()
    date = now.strftime('%y%m%d')

    # Set the limit of X axis.
    if xlim == None:
        minimum = 0
        maximum = len(df[x].unique())
    else:
        minimum = xlim[0]
        maximum = xlim[1]

    # Make dic with updated values -> value + sum from previous values and from it df again.
    pivot = pd.pivot_table(df, index=x, columns=y, values=values, fill_value=0)
    # Reindex the pivot if order_y is specified.
    if order_y == None:
        pivot2=pivot
    else:
        pivot2 = pivot.reindex(order_y, axis=1)
    # Convert pivot table to DataFrame.
    df_pivot= pivot2.reset_index()

    # if you want to have continuous x-axis like 1,2,3,4,5,6,... -. works only for numbers
    if fill_empty_x:
        # If min(df_pivot[x]) is empty skip the plot
        if df_pivot[x].empty:
            print(f'df["x"] are empty after pivot, plotting of the plot {name_to_save} was terminated.')
            write_into_report(f'df["x"] are empty after pivot, plotting of the plot {name_to_save} was terminated.')
            return
        if min(df_pivot[x]) > 0:
            list_of_x = list(range(max(df_pivot[x]) + 1))
        else:
            if max(df_pivot[x]) < 0:
                list_of_x = list(range(min(df_pivot[x]), 2))
            else:
                list_of_x = list(range(min(df_pivot[x]), max(df_pivot[x]) + 1))

        if 0 in list_of_x:
            list_of_x.remove(0)

        # Track the index of 1 in the cleaned-up list
        try:
            middle = list_of_x.index(1)
        except ValueError:
            middle = 0  # fallback just in case

        no_of_columns = len(df_pivot.columns)

        for lx in list_of_x:
            if lx not in df_pivot[x].values:
                empty_row = [lx] + [0] * (no_of_columns - 1)
                df_pivot.loc[len(df_pivot)] = empty_row

        df_pivot.sort_values(x, inplace=True)

    # if fill_empty_x == True:
    #     # check if all columns exists, if not make an empty one
    #     if min(df_pivot[x]) > 0:
    #         list_of_x = list(range(max(df_pivot[x])+1))
    #     else:
    #         list_of_x = list(range(min(df_pivot[x]), max(df_pivot[x])+1))
    #     if 0 in list_of_x:
    #         list_of_x.remove(0)
    #     no_of_columns = len(df_pivot.columns)

    #     if 1 not in list_of_x:
    #         empty_row = [1] + [0] * (no_of_columns - 1)
    #         df_pivot.loc[len(df_pivot)] = empty_row
    #         # append 1 to the list of x values

    #     # try:
    #     middle = list_of_x.index(1)
    #     # except:
    #     #     middle=0
        
    #     for lx in list_of_x:
    #         if lx not in list(df_pivot[x]):
    #             empty_row = [lx]
    #             for i in range(no_of_columns-1):
    #                 empty_row.append(0)
    #             df_pivot.loc[len(df_pivot)] = empty_row
    #     df_pivot.sort_values(x, inplace=True)

    i = 0
    to_plot_dic = {}
    for row in df_pivot[x]:
        # print('Column is ' + str(row))
        list_of_row_values = []
        value=0
        for column in range(len(df_pivot.iloc[i,1:])):
            value=df_pivot.iloc[i,column+1] + value
            list_of_row_values.append(value)
        to_plot_dic[row] = list_of_row_values
        i+=1
    
    # Sort the bars
    if order_y == None:
        bars = sorted(df_pivot.columns[1:])
    else:
        bars = order_y
    no_of_bars = len(bars)

    df = pd.DataFrame.from_dict(to_plot_dic, orient='index', columns=bars)

    # plot
    sns.set_style(style='whitegrid')                      # makes gid behind the plot
    f, ax = plt.subplots(figsize=(maximum/2, 3))       # set the figure size

    # make color palette from sns.color_palette or own dictionary
    if sns_palette == True:
        if isinstance(palette, str):
            pal = sns.color_palette(palette, no_of_bars)
        else:
            print('Set palette to valid sns.palete. Default hls was used.')
            pal = sns.color_palette('hls', no_of_bars)
    else:
        # if some of the bar does not exist as a key in the dictionary -> use sns palette
        pal = []
        for column in bars:
            if column in palette:
                pal.append(palette[column])
            else:
                pal = sns.color_palette('hls', no_of_bars)
                print('Key ' + str(column) + ' is not present in color dictionary, thus default hls palette was used.') 
                break
    if hatches_dic != None:
        for column in range(no_of_bars,0,-1):              # plot in reverse
            sns.barplot(x=df.index, y=df.iloc[:,column-1], color=pal[column-1], edgecolor='#000000', \
                alpha=1, errorbar=None, linewidth=0.2, hatch=hatches_dic[list(df.columns)[column-1]]) 
    else:
        for column in range(no_of_bars,0,-1):              # plot in reverse
            sns.barplot(x=df.index, y=df.iloc[:,column-1], color=pal[column-1], edgecolor='#000000', \
                alpha=1, errorbar=None, linewidth=0.2)        

    plt.title(f'{title} \n ({date}) \n {data_source}')
    plt.xlabel(xlabel, size=10)
    plt.ylabel(ylabel, size=15)
    bars = np.flip(bars)                            # change the bars order as they were plotted 
    plt.legend(labels=bars, loc='center left', bbox_to_anchor=(1, 0.5))
    if xlim != None:
        if maximum <= 0:
            maximum = 1
        plt.xlim(middle+minimum-0.5,middle+maximum-0.5)
    if ylim != None:
        plt.ylim(ylim[0], ylim[1])
    if plot_vertical_line_middle == True:
        plt.axvline(x=middle-0.5, color='black', ls='--', label='blunt cut site')
    if y_as_perc == True:
        # set decimal points dependent on the ylim maxima
        ymax = ax.get_ylim()[1]
        if ymax >= 10:
            decimals_value = 0
        elif ymax >= 1:
            decimals_value = 1
        else:
            decimals_value = 2
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(yaxis_formater, decimals=decimals_value))
        # Add text below chart
    if text_bellow_chart != None:
        # at = AnchoredText(text_bellow_chart, prop=dict(size=12), frameon=True, loc='upper left'\
        #     , bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=f.transFigure, borderpad=0.)
        # at.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
        # ax.add_artist(at)
        props = {'backgroundcolor':'lightgray', 'size':12}  # bbox features
        at = AnchoredText(text_bellow_chart, frameon=True, loc='upper left'\
            , bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=f.transFigure, borderpad=0., prop=props)
        at.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
        ax.add_artist(at)
            # plt.text(0, -0.35, text_bellow_chart, fontsize=14, transform=plt.gcf().transFigure)
    plt.savefig(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.' + output_format), dpi=dpi, bbox_inches='tight')
    write_into_report(f'Stacked vertical barchart {name_to_save} was created.')
    if show_plots == True:
        return plt.show()
    return plt.close()

# def plot_sns_boxplot_and_print(data, x, y, number_of_columns, project_output_dir\
#                                , hue=None, title=None, xtickt_rotation=0\
#                                , ylabel=None, order=None, integrated_palette=False\
#                                , transparency = 0.8 \
#                                , hatches=False, hatches_dic=None, hatches_color_dic=None\
#                                , yaxis_formater=1, name_to_save='Boxplot', data_source = None\
#                                , y_min=None, y_max=None, y_legend_coordinate=0, plot_legend=False\
#                                , add_stat=False, stat_pairs=None, stat_test='t-test_ind'\
#                                , text_bellow_chart=None, output_format='png', dpi=200, show_plots=False):
#     """Plots boxplot and saves it.
#     data: DataFrame
#     x: x to plot
#     y: y to plot
#     number_of_columns: int() to calculate the chart width
#     hue: hue or None
#     title: title, default None
#     xtickt_rotation: int(), default 0
#     ylabel: ylabel, default None
#     order: order of x, default None
#     integrated_palette:True use integrated dictionary; False -> makes default colors; enter own palette dic plain for color or nasted {'key':{'color':'gray', 'hatch_color':'', 'hatch':'////'}}
#     transparency = int() 0-1, default 0.8
#     hatches: True/False, if you want to make hatches, default False
#     hatches_dic: use when not nasted dic as integrated_palette is used
#                 for example {'key':'gray'}, deault None
#     hatches_color_dic: fuse when not nasted dic as integrated_palette is used
#                     for example {'key':'gray'}, deault None
#     yaxis_formater: usually 1 of 100, to make yaxis in %, default 1
#     project_output_dir: dir for output
#     name_to_save: hae to contain full path with name of the boxplot (date and format is added automatialy)
#     data_source: like 'Přibylová et al., (2021)'
#     y_min: int(), default None
#     y_max: int(), default None
#     y_legend_coordinate: int() coordinates of figure are 0-1, default 0
#     plot_legend: True/False, default False
#     add_stat: add t-test_ind statistic annotation, default False
#     stat_pairs: [(('pairs'), ('to compare')), (('pairs'), ('to compare')), ...]
#     stat_test: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, default t-test_ind
#     text_bellow_chart: str() to comment chart bellow, default None
#     output_format: png, svg,..., default 'png'
#     dpi: int(), default 200
#     show_plots: True/False, default False"""
    
#     # TODO move the dict as a universal parameter for all charts
#     color_dic = {'stagg2nt':{'color':'#E890B8','hatch_color':'gray', 'hatch':''}\
#             ,'stagg2nt_likely':{'color':'#E890B8','hatch_color':'gray', 'hatch':'//'}\
#             ,'stagg1nt':{'color':'#ffaaaa','hatch_color':'gray', 'hatch':''}\
#             ,'stagg1nt_likely':{'color':'#ffaaaa','hatch_color':'gray', 'hatch':'////'}\
#             ,'blunt':{'color':'#6ea7ff','hatch_color':'gray', 'hatch':''}\
#             ,'blunt_likely':{'color':'#6ea7ff','hatch_color':'gray', 'hatch':'\\\\'}\
#             ,'cannot_be_determined':{'color':'white','hatch_color':'gray', 'hatch':''}\
#             ,'can_not_be_determined':{'color':'white','hatch_color':'gray', 'hatch':''}\
#             ,'gap_out_of_cut':{'color':'gray','hatch_color':'', 'hatch':''}\
#             ,'IN':{'color':'#091f92','hatch_color':'', 'hatch':''}\
#             ,'INs':{'color':'#091f92','hatch_color':'', 'hatch':''}\
#             ,'DEL':{'color':'#ba8759','hatch_color':'', 'hatch':''}\
#             ,'DELs':{'color':'#ba8759','hatch_color':'', 'hatch':''}\
#             ,'SUB':{'color':'#ffcc33','hatch_color':'', 'hatch':''}\
#             ,'OTHERS':{'color':'#808080','hatch_color':'', 'hatch':''}\
#             ,'no_mutation':{'color':'#111111','hatch_color':'', 'hatch':''}\
#             ,'distal':{'color':'#f0d271','hatch_color':'', 'hatch':''}\
#             ,'distal_both':{'color':'#C4C4C4','hatch_color':'', 'hatch':''}\
#             ,'distal_both_proximal':{'color':'#979797','hatch_color':'', 'hatch':''}\
#             ,'both':{'color':'#808080','hatch_color':'', 'hatch':''}\
#             ,'both_proximal':{'color':'#C4C4C4','hatch_color':'', 'hatch':''}\
#             ,'proximal':{'color':'#f5ad62','hatch_color':'', 'hatch':''}\
#             ,'PAM_distal_nt_in_MH':{'color':'#52D8FF','hatch_color':'', 'hatch':''}\
#             ,'PAM_proximal_nt_in_MH':{'color':'#FADA52','hatch_color':'', 'hatch':''}\
#             ,'no_MH_likely':{'color':'#A9A9A9','hatch_color':'', 'hatch':''}\
#             ,'no_MH':{'color':'#D3D3D3','hatch_color':'', 'hatch':''}\
#             ,'A':{'color':'#70c936','hatch_color':'', 'hatch':''}\
#             ,'C':{'color':'#ff0000','hatch_color':'', 'hatch':''}\
#             ,'G':{'color':'#ffd42a','hatch_color':'', 'hatch':''}\
#             ,'T':{'color':'#9955ff','hatch_color':'', 'hatch':''}}
    
#     # Check if there are data in the data[x], if not skip making the plot.
#     if data[x].notnull().any() == False:
#         print(f'data["x"] are empty, plotting of the plot {name_to_save} was terminated.')
#         write_into_report(f'data["x"] are empty, plotting of the plot {name_to_save} was terminated.', destination_folder=project_output_dir)
#         return
    
#     now = datetime.now()
#     date = now.strftime('%y%m%d')
#     # Change matplotlib style and define figure size
#     sns.set(style='whitegrid')
#     f, ax = plt.subplots(figsize=(number_of_columns, 4))
#     # Define which palette will be used
#     if integrated_palette == True:
#         color_dic = color_dic
#     else:
#         color_dic = integrated_palette
#     # Make color list from dictionary
#     color_list=[]
#     hatches_list=[]
#     edge_color_list=[]
#     if order==None:
#         order = list(data[x].unique())    
#     for order_key in order:
#         if ((any(isinstance(i, dict) for i in color_dic.values()))) == True:
#             color_list.append(color_dic[order_key]['color'])
#             hatches_list.append(color_dic[order_key]['hatch'])
#             edge_color_list.append(color_dic[order_key]['hatch_color'])
#         else:
#             color_list = color_dic
#             hatches_list = hatches_dic
#             edge_color_list = hatches_color_dic

#     #  Plot Boxplot in order (True/False)
#     if order == None:
#         boxplot = sns.boxplot(x=x, y=y, hue=hue, data=data, palette=color_list, \
#                 showmeans=True, meanprops={'marker':'x'
#                                         ,'markerfacecolor':'black'
#                                         , 'markeredgecolor':'black'
#                                         , 'markersize':'8'})

#     else:
#         boxplot = sns.boxplot(x=x, y=y, hue=hue, data=data, palette=color_list, \
#             showmeans=True, meanprops={'marker':'x'
#                                     ,'markerfacecolor':'black'
#                                     , 'markeredgecolor':'black'
#                                     , 'markersize':'8'}
#                                     , order=order)
#     # Statistic
#     if add_stat == True:
#         add_stat_annotation(ax, data=data, x=x, y=y\
#             , box_pairs=stat_pairs,\
#             test=stat_test, text_format='star', loc='outside', verbose=1, fontsize='large')
#     # Legend
#     if plot_legend == True:
#         plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     plt.title(title + '\n(' + date + ')\n' + data_source)
#     plt.ylabel(ylabel)
#     plt.ylim(y_min, y_max)
#     plt.xticks(rotation=xtickt_rotation, ha='left', rotation_mode='anchor')
#     # Add transparency to colors 
#     for patch in boxplot.artists:
#         r, g, b, a = patch.get_facecolor()
#         patch.set_facecolor((r, g, b, transparency))
#     # set hatches if wanted
#     if hatches == True:
#         # Loop over the bars
#         for i,thisbar in enumerate(boxplot.artists):
#             # Set a different hatch for each column
#             thisbar.set_hatch(hatches_list[i])
#             thisbar.set_edgecolor(edge_color_list[i])
    
#     ax.yaxis.set_major_formatter(mtick.PercentFormatter(yaxis_formater, decimals=0))
#     # makes gaps between boxes when hue is used
#     adjust_box_widths(f, 0.9, len(order)) 
#     # Add text below chart
#     if text_bellow_chart != None:
#         # at = AnchoredText(text_bellow_chart, prop=dict(size=12), frameon=True, loc='upper left'\
#         #     , bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=f.transFigure, borderpad=0.)
#         # at.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
#         # ax.add_artist(at)
#         props = {'backgroundcolor':'lightgray', 'size':12}  # bbox features
#         at = AnchoredText(text_bellow_chart, frameon=True, loc='upper left'\
#             , bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=f.transFigure, borderpad=0., prop=props)
#         at.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
#         ax.add_artist(at)

#     plt.savefig(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.' + output_format), dpi=dpi, bbox_inches='tight')
#     write_into_report(f'Boxplot plot {name_to_save} was created.', destination_folder=project_output_dir)
#     if show_plots == True:
#         return plt.show()
#     return plt.close()