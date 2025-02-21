# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
This code is for making dictionaries for plotting colour palettes and hatches"""
def make_color_dic_CutType(list_of_events):
    ''' Makes color dictionary for plotting cut types.
    Parameters:
    list_of_events (list): based on function check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, 'BluntStag', cut_maximum, cut_minimum)

    Returns:
    color_dic
    '''
    # Initialize an empty dictionary to store the results
    color_dic = {}
    # Iterate over each event in the list
    for l in list_of_events:
        if l.startswith('0'):
            color_dic[l] = '#6ea7ff'
        elif l.startswith('-1'):
            color_dic[l] = '#ffaaaa'
        elif l.startswith('-'):
            color_dic[l] = '#E890B8'
        elif l[0].isdigit():
            color_dic[l] = 'yellow'
        elif l.startswith(('C', 'O', 'X')):
            color_dic[l] = 'gray'
        else:
            color_dic[l] = 'red'  # Default case if it doesn't match any criteria
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
        if l.startswith('C'):
            hatches_dic[l] = ''
        elif l.startswith('O'):
            hatches_dic[l] = 'oo'
        elif l.startswith('X'):
            hatches_dic[l] = 'xxx'
        elif l.endswith('0'):
            hatches_dic[l] = ''
        elif l.endswith('1'):
            if l.startswith('0'):
                hatches_dic[l] = '--'
            elif l.startswith('-1'):
                hatches_dic[l] = '//'
            elif l.startswith('-'):
                hatches_dic[l] = '////'
            else:
                hatches_dic[l] = '\\'
        else:
            hatches_dic[l] = '*'
    return hatches_dic