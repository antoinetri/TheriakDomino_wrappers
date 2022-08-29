#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################################
# Read plt files from Theriak-Domino pseudosection
# Plotting phase diagram with variance ready for vector data editing
# Author: Antoine Triantafyllou 2022-08 - Antoine.Triantafyllou@univ-lyon1.fr
# Associated paper: Triantafyllou et al. (2022) - Geology - Add title and vol.
#############################################################################


# Import useful libraries
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib


def import_param(filename):
    """read the plt file produced by domino/guzzler and extract phases diagram parameters"""
    f = open(filename, "r")
    lines = list(f)

    # gather P-T limits
    cond = lines[2]
    p_t_cond = cond.split('   ')
    p_cond = [float(p_t_cond[3].replace(' ', '')), float(p_t_cond[4].replace(' ', ''))]
    t_cond = [float(p_t_cond[1].replace(' ', '')), float(p_t_cond[2].replace(' ', ''))]

    # gather list and number of cations
    chm = lines[4]
    chm_eq = chm.split('    ')[-1]
    chm_line = chm_eq.split('= ')[1]
    el_list = chm_line.split('(')
    count = 0
    for el in el_list:
        if 'O' in el:
            count = count
        else:
            count = count+1
    num_of_cat = count - 1

    # number of row to skip
    num_key = 0
    county = 0
    for line_in in lines:
        county = county + 1
        if "    0    0    0    0" in line_in:
            num_key = num_key + 1
            if num_key == 2:
                break
            else:
                pass
        else:
            pass
    row_to_skip = county - 1

    return p_cond, t_cond, num_of_cat, row_to_skip


def import_reactions(filename, skip_row=25, phase_list=[], el=0):
    """read the plt file produced by domino/guzzler and create a df out of it"""
    f = open(filename, "r")
    lines = list(f)[skip_row::]
    it_temp = []
    all_list = []
    # iterate line by line in txt file and make a list of list including each reaction
    for line_in in lines:
        if "    0    0    0    0" in line_in:
            all_list.append(it_temp)
            it_temp = []
        else:
            it_temp.append(line_in)
    all_list.append(it_temp)

    # iterate through each list of list to sort assemblage vs pt values
    list_assemblage, list_p, list_t, reaction_number, list_num_phases = [], [], [], [], []
    number = 0
    for reaction in all_list[1::]:
        number = number + 1
        counter = 0
        # gather assemblage names from line and make a counter to know how many lines of PT values we have for
        # this reaction
        for it in range(len(reaction)):
            if any(phase in reaction[it] for phase in phase_list):
                assemblage = reaction[it].split('  0')[-1]
                num_phases = len(assemblage[:-1].split(' '))
                list_assemblage.append(assemblage[:-1])
                list_num_phases.append(num_phases)
            else:
                counter = counter + 1
        pt_last2 = []

        # just create this exception for reaction having only one line of PT values
        if counter == 1:
            pt_1 = reaction[0].split('   ')
            pt_2 = [w.replace(' ', '') for w in pt_1]
            pt_3 = [y.rstrip('\n') for y in pt_2]
            pt_last2 = [float(i) for i in pt_3[1::]]
            pressure = pt_last2[1::2]
            temperature = pt_last2[::2]
            list_p.append(pressure)
            list_p.append(pressure)
            list_t.append(temperature)
            list_t.append(temperature)

        # otherwise, just iterate through these lines only and concatenate PT value in one list
        else:
            for it2 in range(counter-1):
                pt_1 = reaction[it2].split('   ')
                pt_2 = [w.replace(' ', '') for w in pt_1]
                pt_3 = [y.rstrip('\n') for y in pt_2]
                pt_last = [float(i) for i in pt_3[1::]]
                pt_last2.extend(pt_last)
            # distinguish P and T values in list
            pressure = pt_last2[1::2]
            temperature = pt_last2[::2]
            list_p.append(pressure)
            list_p.append(pressure)
            list_t.append(temperature)
            list_t.append(temperature)
        reaction_number.append(number)
        reaction_number.append(number)

    # Put all these lists in a single dataframe
    df_all = pd.DataFrame(columns=['number', 'assemblage', 'p_list', 't_list', 'num_phases', 'variance'])
    df_all['number'] = reaction_number
    df_all['assemblage'] = list_assemblage
    df_all['p_list'] = list_p
    df_all['t_list'] = list_t
    df_all['t_list'] = list_t
    df_all['num_phases'] = list_num_phases
    df_all['variance'] = el - df_all['num_phases'] + 2
    return df_all


def plot_polygons(df_data, p_range, t_range, color_v='Greens', label_mode=False):
    """plot all these data"""
    import alphashape
    from descartes import PolygonPatch

    variance_range = [df_data['variance'].min(), df_data['variance'].max()]
    grouped = df_data.groupby(["assemblage"], as_index=False)
    it = 0
    # TODO : add limit conditions to close the ps !!! if it touche p_range, add local min max t_range
    fig = plt.figure(figsize=(10, 10), dpi=80)
    ax = fig.add_subplot(111)
    for name, group in grouped:
        it = it + 1
        x_data, y_data, z_data, poly = [], [], [], []
        for _, row in group.iterrows():
            p_data = row['p_list']
            y_data.extend(p_data)
            s_p = set(p_data)
            if p_range[0] in s_p:
                print("at lower limit of P")
            elif p_range[1] in s_p:
                print("at upper limit of P")
            else:
                pass
            t_data = row['t_list']
            s_t = set(t_data)
            if t_range[0] in s_t:
                print("at lower limit of T")
            elif t_range[1] in s_t:
                print("at upper limit of T")
            else:
                pass
            x_data.extend(t_data)
            var = row['variance']
        print('Reaction # ' + str(it) + ' in progress')

        merged_xy = [(x_data[i], y_data[i]) for i in range(0, len(x_data))]
        # condition to ignore two points list (no polygon)
        print(len(merged_xy))
        if len(merged_xy) > 3:
            alpha_conv = 0.9 * alphashape.optimizealpha(merged_xy, max_iterations=100000)
            if alpha_conv:
                hull = alphashape.alphashape(merged_xy, alpha_conv)
                print(alpha_conv)
            else:
                alpha_conv = 2.0
                hull = alphashape.alphashape(merged_xy, alpha_conv)
                print(alpha_conv)
                print('ok')

            if label_mode:
                y_center = min(y_data) + ((max(y_data) - min(y_data))/2)
                x_center = min(x_data) + ((max(x_data) - min(x_data))/2)
                list_of_phase_in = name.replace(' ', '\n')
                plt.text(x_center, y_center, list_of_phase_in, fontsize=3, ma='center', horizontalalignment='center', verticalalignment='center')
            else:
                pass

            norm = matplotlib.colors.Normalize(vmin=variance_range[0], vmax=float(variance_range[1]))
            if color_v == 'Blues':
                rgba_color = cm.Blues(norm(var), bytes=True)
            elif color_v == 'Greens':
                rgba_color = cm.Greens(norm(var), bytes=True)
            elif color_v == 'Reds':
                rgba_color = cm.Reds(norm(var), bytes=True)
            elif color_v == 'Greys':
                rgba_color = cm.Greys(norm(var), bytes=True)
            elif color_v == 'Oranges':
                rgba_color = cm.Oranges(norm(var), bytes=True)
            else:
                rgba_color = cm.Greens(norm(var), bytes=True)
            rgba_color2 = [float(i) / 255.0 for i in rgba_color]
            ax.add_patch(PolygonPatch(hull, fill=True, color=rgba_color2, alpha=0.8, label=str(name)))
        else:
            pass
    ax.set_ylim(p_range[0], p_range[1])
    ax.set_xlim(t_range[0], t_range[1])
    plt.show()

###########################################################################
# this is where users should modify parameters for treatment and plotting
###########################################################################
if __name__ == '__main__':
    # Add the directory of the domplt.plt file
    filename_plt = r'your\directory\Theriak-Domino_plots\domplt.plt'
    ps_param = import_param(filename_plt)

    # Add in this list most (not all) of the phases that are in your ps (previous ones can still be in the list,no worries)
    phase_list = ['OPXW14', 'CAMPG16', 'AUGG16', '(2)AUGG16', 'GRTW14', 'PLc03', '(2)PLc03', 'sph', 'ru', 'ky', 'q', 'ab', 'H2O', 'LIQMB16', 'ILM00', 'OL11', 'MTSP02', 'EP11', 'BI14', 'MU14', 'CHL14']

    # Set number of cations in system for variance calc
    el = ps_param[2]
    # el = 10 # if user defined, change here

    # Row to skip (change manually if needed)
    skip_row = ps_param[3]

    # precise PT (or TX) ranges for the plot
    p_range = ps_param[0]
    print(p_range)
    t_range = ps_param[1]
    print(t_range)
    # p_range = [5000, 15000] # if user defined, change here [min, max]
    # t_range = [500, 1000]

    # pick a color cmap: Greens, Blues, Reds, Purples for variance
    color_v = 'Greens'

    # choose if you want to plot the phases assemblage or not (False or True) (still bad at the moment but handy)
    label_mode = True

    df_data = import_reactions(filename_plt, skip_row, phase_list, el)
    plot_polygons(df_data, p_range, t_range, color_v, label_mode)

