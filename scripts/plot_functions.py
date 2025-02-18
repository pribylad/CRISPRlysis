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

def plot_sns_violin_and_print(data, x, y, number_of_columns, project_output_dir\
        , hue=None, hue_order=None, hue_only_for_data_points=True\
        , ylabel=None, title=None, xtickt_rotation=0\
        , order=None\
        , integrated_properties_dic=False, color_dic=None\
        , hatches=False, hatches_list=None\
        , change_edge_color=False, edge_color_dic=None\
        , show_datapoints=True, jitter=True, stripplot_offset=-10, dodge=True\
        , show_hue_medians=True\
        , transparency = 0.8 \
        , clip_right=False\
        , yaxis_formater=1, name_to_save="Boxplot", data_source = None\
        , y_min=None, y_max=None\
        , x_min=None, x_max=None\
        , plot_vertical_line_middle=False, middle_line_position=1\
        , plot_legend=False\
        , add_stat=False, stat_pairs=None, stat_test="t-test_ind"\
        , text_bellow_chart=None, y_legend_coordinate=0\
        , output_format="png", dpi=200, show_plots=False):
    """Plots violinplot and saves it.
    data: DataFrame
    x: x to plot
    y: y to plot
    number_of_columns: int() to calculate the chart width
    project_output_dir: dir for output
    hue: hue or None
    hue_order: list or None, default None
    hue_only_for_data_points: True/False, default True
    ylabel: ylabel, default None
    title: title, default None
    xtickt_rotation: int(), default 0
    order: order of x, default None
    integrated_properties_dic:True/False use integrated dictionary for color/hatch/edge_color
    color_dic: dic{} if integrated_properties_dic==False, default None
    hatches: True/False, if you want to make hatches, default False
    hatches_list: list[] if integrated_properties_dic== False, deault None
    change_edge_color: True/False, default False
    edge_color_dic: dic{} if integrated_properties_dic==False, default None
    show_datapoints: True/False, default True
    jitter: True/False/int() for stripplot, default True
    stripplot_offset: int() to move stripplot to the left or right, default -10
    dodge: True/False separated stripplot hues, default True
    show_hue_medians: True/False for stripplot, default True
    transparency = int() 0-1, default 0.8
    clip_right: True/False removes right part of the violine plot, default False
    yaxis_formater: usually 1 of 100, to make yaxis in %, default 1
    name_to_save: hae to contain full path with name of the boxplot (date and format is added automatialy)
    data_source: like "Přibylová et al., (2021)"
    y_min: int(), default None
    y_max: int(), default None
    x_min: int(), default None
    x_max: int(), default None
    plot_vertical_line_middle: True/False, default False
    middle_line_position: int(), default 1
    plot_legend: True/False, default False
    add_stat: add t-test_ind statistic annotation, default False
    stat_pairs: [(("pairs"), ("to compare")), (("pairs"), ("to compare")), ...]
    stat_test: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal, default t-test_ind
    text_bellow_chart: str() to comment chart bellow, default None
    y_legend_coordinate: int() coordinates of figure are 0-1, default 0
    output_format: png, svg,..., default "png"
    dpi: int(), default 200
    report_dir: where should be the report, default os.getcwd()
    show_plots: True/False, default False"""

    properties_dic = {"stagg2nt":{"color":"#E890B8","edge_color":"gray", "hatch":""}\
        ,"stagg2nt_likely":{"color":"#E890B8","edge_color":"gray", "hatch":"//"}\
        ,"stagg1nt":{"color":"#ffaaaa","edge_color":"gray", "hatch":""}\
        ,"stagg1nt_likely":{"color":"#ffaaaa","edge_color":"gray", "hatch":"////"}\
        ,"blunt":{"color":"#6ea7ff","edge_color":"gray", "hatch":""}\
        ,"blunt_likely":{"color":"#6ea7ff","edge_color":"gray", "hatch":"\\\\"}\
        ,"cannot_be_determined":{"color":"white","edge_color":"gray", "hatch":""}\
        ,"can_not_be_determined":{"color":"white","edge_color":"gray", "hatch":""}\
        ,'distal_stagg2nt_likely':{"color":"#E890B8","edge_color":"gray", "hatch":"\\\\"}\
        ,'distal_stagg2nt':{"color":"#E890B8","edge_color":"gray", "hatch":""}\
        ,'distal_stagg1nt_likely':{"color":"#ffaaaa","edge_color":"gray", "hatch":"\\\\"}\
        ,'distal_stagg1nt':{"color":"#ffaaaa","edge_color":"gray", "hatch":""}\
        ,'distal_blunt_likely':{"color":"#6ea7ff","edge_color":"gray", "hatch":"\\\\"}\
        ,'distal_blunt':{"color":"#6ea7ff","edge_color":"gray", "hatch":""}\
        ,'proximal_blunt':{"color":"#6ea7ff","edge_color":"gray", "hatch":""}\
        ,'proximal_blunt_likely':{"color":"#6ea7ff","edge_color":"gray", "hatch":"//"}\
        ,'proximal_stagg1nt':{"color":"#ffaaaa","edge_color":"gray", "hatch":""}\
        ,'proximal_stagg1nt_likely':{"color":"#ffaaaa","edge_color":"gray", "hatch":"//"}\
        ,'proximal_stagg2nt':{"color":"#E890B8","edge_color":"gray", "hatch":""}\
        ,'proximal_stagg2nt_likely':{"color":"#E890B8","edge_color":"gray", "hatch":"//"}\
        ,"gap_out_of_cut":{"color":"gray","edge_color":"gray", "hatch":""}\
        ,"IN":{"color":"#091f92","edge_color":"gray", "hatch":""}\
        ,"INs":{"color":"#091f92","edge_color":"gray", "hatch":""}\
        ,"DEL":{"color":"#ba8759","edge_color":"gray", "hatch":""}\
        ,"DELs":{"color":"#ba8759","edge_color":"gray", "hatch":""}\
        ,"SUB":{"color":"#ffcc33","edge_color":"gray", "hatch":""}\
        ,"OTHERS":{"color":"#808080","edge_color":"gray", "hatch":""}\
        ,"no_mutation":{"color":"#111111","edge_color":"gray", "hatch":""}\
        ,"distal":{"color":"#f0d271","edge_color":"gray", "hatch":""}\
        ,"distal_both":{"color":"#C4C4C4","edge_color":"#f0d271", "hatch":"\\\\\\"}\
        ,"distal_both_proximal":{"color":"#979797","edge_color":"black", "hatch":"xx"}\
        ,"both":{"color":"#808080","edge_color":"black", "hatch":""}\
        ,"both_proximal":{"color":"#C4C4C4","edge_color":"#f5ad62", "hatch":"////"}\
        ,"proximal":{"color":"#f5ad62","edge_color":"gray", "hatch":""}\
        ,"PAM_distal_nt_in_MH":{"color":"#52D8FF","edge_color":"gray", "hatch":""}\
        ,"PAM_proximal_nt_in_MH":{"color":"#FADA52","edge_color":"gray", "hatch":""}\
        ,"no_MH_likely":{"color":"#A9A9A9","edge_color":"gray", "hatch":""}\
        ,"no_MH":{"color":"#D3D3D3","edge_color":"gray", "hatch":""}\
        ,"A":{"color":"#70c936","edge_color":"gray", "hatch":""}\
        ,"C":{"color":"#ff0000","edge_color":"gray", "hatch":""}\
        ,"G":{"color":"#ffd42a","edge_color":"gray", "hatch":""}\
        ,"T":{"color":"#9955ff","edge_color":"gray", "hatch":""}\
        ,True:{"color":"#D3D3D3","edge_color":"gray", "hatch":""}\
        ,False:{"color":"#B8BEDB","edge_color":"gray", "hatch":""}\
        ,"True":{"color":"#D3D3D3","edge_color":"gray", "hatch":""}\
        ,"False":{"color":"#B8BEDB","edge_color":"gray", "hatch":""}\
        ,"unmethylated":{"color":"#000000","edge_color":"gray", "hatch":""}\
        ,"CG_only":{"color":"#4D4D4D","edge_color":"gray", "hatch":""}\
        ,"heterochromatin":{"color":"#999999","edge_color":"gray", "hatch":""}\
        ,"RdDM":{"color":"#fafafa","edge_color":"gray", "hatch":""}}

    # Check if there are data in the data[x], if not skip making the plot.
    if data[x].notnull().any() == False:
        print("data['x'] are empty, plotting of the plot {} was terminated.".format(name_to_save))
        write_into_report("data['x'] are empty, plotting of the plot {} was terminated.".format(name_to_save), destination_folder=project_output_dir)
        return

    now = datetime.now()
    date = now.strftime("%y%m%d")
    # Change matplotlib style and define figure size
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(number_of_columns, 4))

    # This script can not make hatches when usin hue for violinplot
    if hue != None:
        if hue_only_for_data_points == False:
            hatches=False
            print("This script can not make hatches when usin hue for violinplot. Hatches was set to False.")
    # Make color,edge and hatch dic from nested dictionary
    if integrated_properties_dic == True:
        color_dic = {}
        edge_color_list = []
        hatches_list = []
        test = data.groupby(x).value.mean().to_dict()
        for item in list(data[x].unique()):
            if item in properties_dic.keys():
                color_dic[item] = properties_dic[item]["color"]
                if test[item] != 0:
                    edge_color_list.append(properties_dic[item]["edge_color"])
                    hatches_list.append(properties_dic[item]["hatch"])
            if item not in properties_dic.keys():
                print("The " + str(item) + " is not in properties_dic of plot_sns_violin_and_print function. The dic contains:" + str(color_dic))
        if hue != None:
            for item in list(data[hue].unique()):
                if item in properties_dic.keys():
                    color_dic[item] = properties_dic[item]["color"]
                if item not in properties_dic.keys():
                    print("The " + str(item) + " is not in properties_dic of plot_sns_violin_and_print function.")
    else:   
        color_dic = color_dic
        hatches_list = hatches_list
        edge_color_dic = edge_color_dic 

    # to make hue color only for stripplot and not for violin
    if hue_only_for_data_points == False:
        hue_violin = hue
    else:
        hue_violin=None

    # Set withts for plotted boxplots
    width_list = [0.15, 0.9]
    if hue_only_for_data_points == False:
        if hue != None:
            if len(data[hue].unique()) == 4:
                width_list = [0.8, 0.4]
            else:
                width_list = [0.8, 0.4]
                print("Problem with boxplot wisht could be in settings width_list.")
    
    # Plot Violinplot in order (True/False) with stripplot (True/False)
    if order == None:
        # Make Violinplot
        # Change color_dic keys to 'txt' as the violinplot bellow makes from True -> 'True'
        # and the plot does not work without the conversion -> ideas how to solve it?
        color_dic_txt = {str(key): value for key, value in color_dic.items()}
        violinplot = sns.violinplot(data=data, x=x, y=y, hue=hue_violin, hue_order=hue_order, density_norm="width", palette=color_dic_txt, cut=0, inner=None, zorder=0)
        ## Set hatches if wanted
        if hatches == True:
            ihatch = iter(hatches_list)
            hatch = [i.set_hatch(next(ihatch)) for i in ax.get_children() if isinstance(i, matplotlib.collections.PolyCollection)]
        # Clip the right half of each violin.
        if clip_right == True:
            for item in ax.collections:
                x0, y0, width, height = item.get_paths()[0].get_extents().bounds
                item.set_clip_path(plt.Rectangle((x0-0.1, y0), (width/2)+0.1, height,
                                angle=0, transform=ax.transData))
        # Make stripplot to show individual data
        if show_datapoints == True:
            # move all of the datapoints a bit to the left or right
            offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
            trans = plt.gca().transData
            stripplot = sns.stripplot(x=x, y=y, hue=hue, hue_order=hue_order, data=data, jitter=jitter\
                , order=None, dodge=dodge\
                , palette=color_dic\
                , size=3, linewidth=0.5\
                , transform=trans+offset(stripplot_offset)\
                , zorder=20)
            # Add median for stripplot
            if show_hue_medians == True:
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
        # Make inner boxplot with mean as "x"
        boxplot = sns.boxplot(x=x, y=y, hue=hue_violin, hue_order=hue_order, data=data, width=width_list[0]\
            , showfliers=False, showmeans=True\
            , meanprops=dict(marker='x', markeredgecolor='black', markerfacecolor='black'\
                , markersize=5, linewidth=1.0, zorder=14)\
            , boxprops=dict(facecolor="gray", edgecolor="black", linewidth=1.5, alpha=1, zorder=10)\
            , whiskerprops=dict(color="black", linewidth=1.5, zorder=11)\
            , capprops=dict(color="black", linewidth=1.5, zorder=12)\
            , medianprops=dict(color="black", linewidth=1.5, zorder=13))
    else:
        # Make Violinplot
        violinplot = sns.violinplot(x=x, y=y, hue=hue_violin, hue_order=hue_order, data=data, density_norm="width"\
            , palette=color_dic\
            , cut=0, inner=None, order=order, zorder=0)
        ## Set hatches if wanted
        if hatches == True:
            hatches_list = []
            for item in order:
                if test[item] != 0:
                    hatches_list.append(properties_dic[item]["hatch"])
            ihatch = iter(hatches_list)
            hatch = [i.set_hatch(next(ihatch)) for i in ax.get_children() if isinstance(i, matplotlib.collections.PolyCollection)]
        # Clip the right half of each violin.
        if clip_right == True:
            for item in ax.collections:
                x0, y0, width, height = item.get_paths()[0].get_extents().bounds
                item.set_clip_path(plt.Rectangle((x0-0.1, y0), (width/2)+0.1, height,
                                angle=0, transform=ax.transData))
        # Make stripplot to show individual data
        if show_datapoints == True:
            # move all of the datapoints a bit to the left or right
            offset = lambda p: transforms.ScaledTranslation(p/72.,0, plt.gcf().dpi_scale_trans)
            trans = plt.gca().transData
            stripplot = sns.stripplot(x=x, y=y, hue=hue, hue_order=hue_order, data=data, jitter=jitter\
                , order=order, dodge=dodge\
                , palette=color_dic\
                , size=3, linewidth=0.5\
                , transform=trans+offset(stripplot_offset)\
                , zorder=20)
            # Add mean for stripplot
            if show_hue_medians == True:
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
        # Make inner boxplot with mean as "x"
        boxplot = sns.boxplot(x=x, y=y, hue=hue_violin, hue_order=hue_order, data=data\
            , width=width_list[0]\
            , order=order, showfliers=False, showmeans=True\
            , meanprops=dict(marker='x', markeredgecolor='black', markerfacecolor='black'\
                , markersize=5, linewidth=1.0, zorder=14)\
            , boxprops=dict(facecolor="gray", edgecolor="black", linewidth=1.5, alpha=1, zorder=10)\
            , whiskerprops=dict(color="black", linewidth=1.5, zorder=11)\
            , capprops=dict(color="black", linewidth=1.5, zorder=12)\
            , medianprops=dict(color="black", linewidth=1.5, zorder=13))
    
    # Change edge color if wanted
    if change_edge_color == True:
        if order!=None:
            edge_color_list = []
            for item in order:
                if test[item] != 0:
                    edge_color_list.append(properties_dic[item]["edge_color"])
        iedge_color = iter(edge_color_list)
        for i in range(0,len(edge_color_list)):
            ax.collections[i].set_edgecolor(next(iedge_color))
    
    # Statistic
    if add_stat == True:
        add_stat_annotation(ax, data=data, x=x, y=y\
            , box_pairs=stat_pairs,\
            test=stat_test, text_format='star', loc='outside', verbose=1, fontsize='large')
    
    # Plot vertical middle line
    if plot_vertical_line_middle == True:
        plt.axvline(x=middle_line_position-0.5, color="black", ls="--", label='blunt cut site')
    
    # Legend
    if plot_legend == True:
        if show_datapoints == True:
            if show_hue_medians == True:
                # Plot legend only from the stripplot and not boxplot (they are identical)
                handles, labels = ax.get_legend_handles_labels()
                plt.legend(handles[int(len(handles)/2):int(len(handles))],\
                    labels[int(len(labels)/2):int(len(labels))],\
                    loc="center left", bbox_to_anchor=(1, 0.5))    
            else:
                plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        else:
            if show_hue_medians == True:
                # Plot legend only from the stripplot and not boxplot (they are identical)
                handles, labels = ax.get_legend_handles_labels()
                plt.legend(handles[0:int(len(handles)/2)],\
                    labels[int(len(labels)/2):int(len(labels))],\
                    loc="center left", bbox_to_anchor=(1, 0.5))    
            else:
                # plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
                # Plot legend only from the stripplot and not boxplot (they are identical)
                handles, labels = ax.get_legend_handles_labels()
                plt.legend(handles[0:int(len(handles)/2)],\
                    labels[int(len(labels)/2):int(len(labels))],\
                    loc="center left", bbox_to_anchor=(1, 0.5))    
    # Plot labels, set lims, etc.
    plt.title(title + "\n(" + date + ")\n" + data_source)
    plt.ylabel(ylabel)
    plt.ylim(y_min, y_max)
    plt.xlim(x_min, x_max)
    plt.xticks(rotation=xtickt_rotation, ha="left", rotation_mode="anchor")

    # Add transparency to colors
    plt.setp(ax.collections, alpha=transparency)

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(yaxis_formater, decimals=0))
    # makes gaps between boxes when hue is used
    adjust_box_widths(fig, width_list[1])
    # Add text below chart
    if text_bellow_chart is not None:
        props = {"backgroundcolor":"lightgray", "size":12}  # bbox features
        at = AnchoredText(text_bellow_chart, frameon=True, loc='upper left'\
            , bbox_to_anchor=(0, y_legend_coordinate), bbox_transform=fig.transFigure, borderpad=0., prop=props)
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(at)

    plt.savefig(Path(str(project_output_dir) + "/" + name_to_save + "_" + date + "." + output_format), dpi=dpi, bbox_inches="tight")
    write_into_report("Violin plot {0} was created.".format(name_to_save), destination_folder=project_output_dir)
    if show_plots is True:
        return plt.show()
    return plt.close()

# def plot_sns_boxplot_and_print(data, x, y, number_of_columns, project_output_dir\
#                                , hue=None, title=None, xtickt_rotation=0\
#                                , ylabel=None, order=None, integrated_palette=False\
#                                , transparency = 0.8 \
#                                , hatches=False, hatches_dic=None, hatches_color_dic=None\
#                                , yaxis_formater=1, name_to_save="Boxplot", data_source = None\
#                                , y_min=None, y_max=None, y_legend_coordinate=0, plot_legend=False\
#                                , add_stat=False, stat_pairs=None, stat_test="t-test_ind"\
#                                , text_bellow_chart=None, output_format="png", dpi=200, show_plots=False):
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
#             # Set a different hatch for each bar
#             thisbar.set_hatch(hatches_list[i])
#             thisbar.set_edgecolor(edge_color_list[i])
    
#     ax.yaxis.set_major_formatter(mtick.PercentFormatter(yaxis_formater, decimals=0))
#     # makes gaps between boxes when hue is used
#     adjust_box_widths(f, 0.9) 
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