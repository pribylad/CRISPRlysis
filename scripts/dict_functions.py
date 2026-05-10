# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
This code is for making dictionaries for plotting colour palettes and hatches
"""


def make_color_dic_CutType(list_of_events):
    ''' Makes color dictionary for plotting cut types.
    Parameters:
    list_of_events (list): based on function check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, 'BluntStag', cut_maximum, cut_minimum)

    Returns:
    color_dic
    '''
    from random import randint
    color_dic = {'-2nt staggered cut':'#64b96a',
                '-2nt staggered cut with complementarity':'#64b96a',
                '-1nt staggered cut':'#ff6a2a',
                '-1nt staggered cut with complementarity':'#ff6a2a',
                'blunt cut':'#ffe46a',
                'blunt cut with complementarity':'#ffe46a',
                'type of cut cannot be determined':'#CCCCCC',
                'type of cut cannot be determined with complementarity':'#CCCCCC',
                'DEL span analysed window':'#464646'
            }
    # If list_of_events not in color_dic, add it with random color
    for l in list_of_events:
        if l not in color_dic:
            # Generate random HEX color
            color_dic[l] = f'#{randint(0, 0xFFFFFF):06x}'
            print(f'Note: in make_color_dic_CutType function. Parameter: {l} is not in the colour decision list. It has been added with random color {color_dic[l]}.')

    return color_dic

def make_hatches_dic_CutType(list_of_events):
    """
    Makes hatches dictionary for plotting cut types.
    Parameters:
    list_of_events (list): based on function check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, 'BluntStag', cut_maximum, cut_minimum)

    Returns:
    hatches_dic
    """
    hatches_dic = {}
    for l in list_of_events:
        if 'complementarity' in l:
            hatches_dic[l] = '\\\\'
        else:
            hatches_dic[l] = ''
    return hatches_dic

def make_color_dic_blunt_stagg_orientation(list_of_events):
    """
    Makes color dictionary for plotting blunt stagg orientation.
    Parameters:
    list_of_events (list): ['U_-2_0', 'U_-2_1', 'U_-1_0', ...]

    Returns:
    color_dic
    """
    from random import randint
    color_dic={'upstream -2nt staggered cut':{'color':'#64b96a','edge_color':'gray', 'hatch':''},
                'upstream -2nt staggered cut with complementarity':{'color':'#64b96a','edge_color':'gray', 'hatch':'\\\\'},
                'upstream -1nt staggered cut':{'color':'#ff6a2a','edge_color':'gray', 'hatch':''},
                'upstream -1nt staggered cut with complementarity':{'color':'#ff6a2a','edge_color':'gray', 'hatch':'\\\\'},
                'upstream blunt cut':{'color':'#ffe46a','edge_color':'gray', 'hatch':''},
                'upstream blunt cut with complementarity':{'color':'#ffe46a','edge_color':'gray', 'hatch':'\\\\'},
                'downstream -2nt staggered cut':{'color':'#64b96a','edge_color':'gray', 'hatch':''},
                'downstream -2nt staggered cut with complementarity':{'color':'#64b96a','edge_color':'gray', 'hatch':'\\\\'},
                'downstream -1nt staggered cut':{'color':'#ff6a2a','edge_color':'gray', 'hatch':''},
                'downstream -1nt staggered cut with complementarity':{'color':'#ff6a2a','edge_color':'gray', 'hatch':'\\\\'},
                'downstream blunt cut':{'color':'#ffe46a','edge_color':'gray', 'hatch':''},
                'downstream blunt cut with complementarity':{'color':'#ffe46a','edge_color':'gray', 'hatch':'\\\\'},
                'orientation cannot be determined blunt cut with complementarity':{'color':'#ffe46a','edge_color':'gray', 'hatch':'**'},
                'bidirectional':{'color':'#CCCCCC','edge_color':'gray', 'hatch':''},
                'bidirectional with complementarity':{'color':'#CCCCCC','edge_color':'gray', 'hatch':'\\\\'},
                'DEL span analysed window':{'color':'#464646','edge_color':'gray', 'hatch':''}
                }
    # If list_of_events not in color_dic, add it with random color
    for l in list_of_events:
        if l not in color_dic:
            # Generate random HEX color
            color_dic[l] = f'#{randint(0, 0xFFFFFF):06x}'
            print(f'Note: in make_color_dic_blunt_stagg_orientation function. Parameter: {l} is not in the colour decision list. It has been added with random color {color_dic[l]}.')
    return color_dic

def make_color_dic_Dino(list_of_events):
    ''' Makes color dictionary for plotting cut types in Dino-charts.
    Parameters:
    list_of_events (list): 
        example 1: ['U_-1_1', 'U_0_0', 'D_0_1', 'D_0_0', 'U_-1_0', 'U_-2_0', 'U_-2_1', 'D_X_0']
        example 2: ['-1_1', 'O_1', '-2_0', '-2_1', '0_0', 'C_0', 'C_1', 'X_0', '0_1', 'X_1', 'O_0', '-1_0']

    Returns:
    color_dic
    '''
    # Initialize an empty dictionary to store the results
    color_dic = {}
    # Iterate over each event in the list
    for l in list_of_events:
        split = l.split('_')
        if len(split) == 2:
            if split[0] == str(0):
                color_dic[l] = '#ffe46a'
            elif split[0] == str(-1):
                color_dic[l] = '#ff6a2a'
            elif split[0].startswith('-'):
                color_dic[l] = '#64b96a'
        elif len(split) == 3:
            if split[0].startswith('s'):
                color_dic[l] = '#CCCCCC'
            elif split[1].startswith('X'):
                color_dic[l] = '#CCCCCC'
            elif split[1] == str(0):
                color_dic[l] = '#ffe46a'
            elif split[1] == str(-1):
                color_dic[l] = '#ff6a2a'
            elif split[1].startswith('-'):
                color_dic[l] = '#64b96a'
        else:
            color_dic[l] = 'red'  # Default case if it doesn't match any criteria
            print(f'Error in make_color_dic_Dnino function. Parameter: {l} is not in the colour decision list.')
    return color_dic

def make_hatches_dic_Dino(list_of_events):
    """
    Makes hatches dictionary for plotting cut types in Dino-charts.
    Parameters:
    list_of_events (list): example: ['U_-1_1', 'U_0_0', 'D_0_1', 'D_0_0', 'U_-1_0', 'U_-2_0', 'U_-2_1']
    
    Returns:
    "hatches_dic
    """
    hatches_dic = {}
    # Iterate over each event in the list
    for l in list_of_events:
        split = l.split('_')
        if split[2] == str(0):
            hatches_dic[l] = ''
        elif split[2] == str(1):
            if split[0] == 'U':
                hatches_dic[l] = '\\\\'
            elif split[0] == 'D':
                hatches_dic[l] = '//'
            else:
                hatches_dic[l] = '*'
                print(f'Error in make_hatches_dic_Dino function. Parameter: {l} is not in the hatch decision list.')
    return hatches_dic

def make_color_dic_Dino_end(list_of_events):
    ''' Makes color dictionary for plotting cut types in Dino-charts.
    Parameters:
    list_of_events (list): 
        example: ['U_5_True', 'U_5_False', 'U_3_True', 'U_3_False', 'D_5_True', 'D_5_False', 'D_3_True', 'D_3_False']

    Returns:
    color_dic
    '''
    # Initialize an empty dictionary to store the results
    color_dic = {}
    # Iterate over each event in the list
    for l in list_of_events:
        split = l.split('_')
        if split[1] == str(5):
            if split[2] == 'False':
                color_dic[l] = '#d22d40'
            elif split[2] == 'True':
                color_dic[l] = '#d3838c'
            else:
                color_dic[l] = 'red'  # Default case if it doesn't match any criteria
                print(f'Error in make_color_dic_Dnino_end function. Parameter: {l} is not in the colour decision list.')
        elif split[1] == str(3):
            if split[2] == 'False':
                color_dic[l] = '#003657'
            elif split[2] == 'True':
                color_dic[l] = "#668699"
            else:
                color_dic[l] = 'red'  # Default case if it doesn't match any criteria
                print(f'Error in make_color_dic_Dnino_end function. Parameter: {l} is not in the colour decision list.')
        elif split[1] == str(0):
            if split[2] == 'False':
                color_dic[l] = "#ffe46a"
            elif split[2] == 'True':
                color_dic[l] = "#CCCCCC"
            else:
                color_dic[l] = 'red'  # Default case if it doesn't match any criteria
                print(f'Error in make_color_dic_Dnino_end function. Parameter: {l} is not in the colour decision list.')
        else:
            color_dic[l] = 'red'  # Default case if it doesn't match any criteria
            print(f'Error in make_color_dic_Dnino function. Parameter: {l} is not in the colour decision list.')
    return color_dic

def get_color_dic_35_violin():
    """
    Returns:
    color_dic
    """
    color_dic = {'3\' possible blunt':"#668699",
                 '5\' possible blunt':'#d3838c',
                 'eider 3\' or 5\' possible blunt':'#FFFFFF',
                 '3\' NOT blunt':'#003657',
                 '5\' NOT blunt':'#d22d40',
                 'eider 3\' or 5\' NOT blunt':"#C7C7C7"}
    return color_dic