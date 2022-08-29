#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################################
# Run Theriak ietartively based on chemical system composition from a datatable file
# Author: Antoine Triantafyllou 2022-08 - Antoine.Triantafyllou@univ-lyon1.fr
# Associated paper: Triantafyllou et al. (2022) - Geology - Add title and vol.
#############################################################################

# import some libraries
import pandas as pd
import subprocess
import shutil
import os
import datetime
import matplotlib.pyplot as plt


def import_data(filename, sheet_n):
    """import xls files to df"""
    data_in = pd.read_excel(filename, sheet_name=str(sheet_n))
    return data_in


def theriak_conv(df_conc, rat_Fe2O3FeOt=0.32):
    """convert from oxides composition to theriak domino expression"""
    # list_of_el without Mn
    list_of_el = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'MgO', 'CaO', 'FeOT', 'Fe2O3', 'BaO', 'Na2O', 'K2O', 'SO3', 'O2', 'H2O/LOI', 'CO2', 'B2O3', 'P2O5', 'Li2O']
    mm_file = r'C:\Users\Utilisateur\PycharmProjects\Theriak-Domino_plots\mm.xlsx'

    df_conc["sum_ox"] = df_conc[list_of_el].sum()
    sum_ox = float(df_conc["sum_ox"])
    df_mm = pd.read_excel(mm_file, sheet_name='mm')
    df_cpfu = pd.read_excel(mm_file, sheet_name='cpfu')
    df_opfu = pd.read_excel(mm_file, sheet_name='opfu')
    mole_sum = 0.0
    for el in list_of_el:
        df_conc[el+'_norm'] = 100 * df_conc[el] / sum_ox
        df_conc[el+'_mole'] = df_conc[el+'_norm'] / df_mm[el]
        mole_sum = mole_sum + float(df_conc[el+'_mole'])

    df_conc['FeO_only_norm'] = df_conc['FeOT_norm'] * (1-float(rat_Fe2O3FeOt))
    df_conc['Fe2O3_only_norm'] = df_conc['FeOT_norm'] * float(rat_Fe2O3FeOt) * (df_mm['Fe2O3']/df_mm['FeOT']) * 0.5
    df_conc['FeO_only_mole'] = df_conc['FeO_only_norm'] / float(df_mm['FeOT'])
    df_conc['Fe2O3_only_mole'] = df_conc['Fe2O3_only_norm'] / float(df_mm['Fe2O3'])
    df_conc['Fe2_only_mol'] = df_conc['FeO_only_mole'] * 1
    df_conc['Fe3_only_mol'] = df_conc['Fe2O3_only_mole'] * 2
    df_conc['oxfe2_mol'] = df_conc['FeO_only_mole'] * 1
    df_conc['oxfe3_mol'] = df_conc['Fe2O3_only_mole'] * 3
    df_conc['tot_ox_fe'] = float(df_conc['oxfe2_mol']) + float(df_conc['oxfe3_mol'])
    df_conc['ox_from_fe3'] = float(df_conc['tot_ox_fe']) - float(df_conc['FeOT_mole'])
    df_conc['O_mole'] = df_conc['ox_from_fe3']
    mole_sum = mole_sum + float(df_conc['O_mole'])
    list_of_el.append('O')

    atom_sum = 0.0
    for el in list_of_el:
        df_conc[el + '_mol_pct'] = 100 * df_conc[el +'_mole'] / float(mole_sum)
        df_conc[el+'_atom_prop'] = df_conc[el + '_mol_pct'] * float(df_cpfu[el])
        atom_sum = atom_sum + float(df_conc[el+'_atom_prop'])

    oxygen_sum = 0.0
    for el in list_of_el:
        df_conc[el + '_cation'] = 100 * df_conc[el+'_atom_prop'] / float(atom_sum)
        df_conc[el + '_oxygen'] = df_conc[el + '_cation'] * df_opfu[el] / df_cpfu[el]
        if str(el) != 'H2O/LOI':
            if str(el) != 'O':
                oxygen_sum = oxygen_sum + float(df_conc[el + '_oxygen'])
            else:
                pass
        else:
            pass
    oxygen_no_wat = float(oxygen_sum) + float(df_conc['O_oxygen'])
    # td_exp = str('SI(' + str(round(float(df_conc['SiO2_cation']), 4)) + ')AL(' + str(round(float(df_conc['Al2O3_cation']), 4)) + ')TI(' + str(round(float(df_conc['TiO2_cation']), 4)) + ')MN(' + str(round(float(df_conc['MnO_cation']), 4)) + ')FE(' + str(round(float(df_conc['FeOT_cation']), 4)) + ')MG(' + str(round(float(df_conc['MgO_cation']), 4)) + ')CA(' + str(round(float(df_conc['CaO_cation']), 4)) + ')NA(' + str(round(float(df_conc['Na2O_cation']), 4)) + ')K(' + str(round(float(df_conc['K2O_cation']), 4)) + ')O(' + str(round(oxygen_no_wat, 4)) + ')H(' + str(round(float(df_conc['H2O/LOI_cation']), 4)) + ')O(' + str(round(float(df_conc['H2O/LOI_oxygen']), 4)) + ')')
    #write a theriak expression version herafter is without Mn:
    td_exp = str('SI(' + str(round(float(df_conc['SiO2_cation']), 4)) + ')AL(' + str(round(float(df_conc['Al2O3_cation']), 4)) + ')TI(' + str(round(float(df_conc['TiO2_cation']), 4)) + ')FE(' + str(round(float(df_conc['FeOT_cation']), 4)) + ')MG(' + str(round(float(df_conc['MgO_cation']), 4)) + ')CA(' + str(round(float(df_conc['CaO_cation']), 4)) + ')NA(' + str(round(float(df_conc['Na2O_cation']), 4)) + ')K(' + str(round(float(df_conc['K2O_cation']), 4)) + ')O(' + str(round(oxygen_no_wat, 4)) + ')H(' + str(round(float(df_conc['H2O/LOI_cation']), 4)) + ')O(' + str(round(float(df_conc['H2O/LOI_oxygen']), 4)) + ')')

    comm = str(df_conc['comment'])
    ids = str(df_conc['id'])
    print(str(ids), comm, td_exp)
    return td_exp, ids, comm


def therin_file(td_exp, comment, p, t, path_dir, p_out=9999, t_out=9999, steps=9999):
    """create a txt therin and therend file for theriak"""
    # create therin file
    therin_dir = os.path.join(path_dir, "THERIN.txt")
    f = open(therin_dir, "a")
    f.write("""! -----Version: 05.09.06
! Comments in this file start with ! at position 1.
!----------------------------------------------------------------------------------------""" + '\n')
    f.write(str(t) + '    ' + str(p) + '\n')
    f.write('1    ' + str(td_exp) + '          *   ' + str(comment) + '\n')
    f.close()

    # open and read the file after the appending:
    # f = open(os.path.join(path_dir, "THERIN.txt"), "r")
    # print(f.read())

    # create therend file
    if steps != 9999:
        therend_dir = os.path.join(path_dir, "THEREND.txt")
        f2 = open(therend_dir, "a")
        f2.write("""! -----Version: 05.09.06
! Comments in this file start with ! at position 1.
!----------------------------------------------------------------------------------------""" + '\n')
        f2.write("TP    " + str(t_out) + '    ' + str(p_out) + '    ' + str(steps) + '          *   ' + '\n')
        f2.close()
    else:
        pass
    print(therin_dir, therend_dir)
    return therin_dir, therend_dir


def run_theriak(source_trin, source_trend, file_out="THEREND.txt", name_out='output', loop_id='this_run', folder_name='your_results', td_database="td-6axNCKFMASHTOm45.txt", working_dir=r'C:\TheriakDominoWIN\TD-win-15.r187-20180826'):
    """runnning the theriak app via subprocess and a command line PIPE"""
    td_db = str(td_database)
    stream = subprocess.Popen("cd\\ & cd " + working_dir + " & theriak.exe", stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
    stream.communicate(input="td-6axNCKFMASHTOm45.txt\ntherend.txt\npython_output\n")

    outs, errs = stream.communicate()
    print(outs, errs)
    print("THERIAK RUN COMPLETED !!!")

    therin_file = 'THERIN.txt'
    therend_file = 'THEREND.txt'
    thkout_file = 'thkout.out'
    thktab_file = 'thktab.tab'
    folder_path = os.path.join(working_dir, loop_id, folder_name)
    print("folder_path is:" + folder_path)
    os.makedirs(folder_path)

    #copy files therin
    original1 = os.path.join(working_dir, therin_file)
    target1 = os.path.join(folder_path, therin_file)
    shutil.copyfile(original1, target1)

    #copy files therend
    original2 = os.path.join(working_dir, therend_file)
    target2 = os.path.join(folder_path, therend_file)
    shutil.copyfile(original2, target2)

    #copy files thkout_file
    original3 = os.path.join(working_dir, thkout_file)
    target3 = os.path.join(folder_path, thkout_file)
    shutil.copyfile(original3, target3)

    #copy files thktab_file
    original4 = os.path.join(working_dir, thktab_file)
    target4 = os.path.join(folder_path, thktab_file)
    shutil.copyfile(original4, target4)

    print("output files copied ! here: " + folder_path)
    return target4, folder_path


def run_domino(source_trin, source_trend, file_out="THEREND.txt", name_out='output', loop_id='this_run', folder_name='your_results', td_database="td-6axNCKFMASHTOm45.txt", working_dir=r'C:\TheriakDominoWIN\TD-win-15.r187-20180826'):
    """WORK IN PROGRESS"""
    td_db = str(td_database)
    stream = subprocess.Popen("cd\\ & cd " + working_dir + " & theriak.exe", stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
    stream.communicate(input="td-6axNCKFMASHTOm45.txt\ntherend.txt\npython_output\n")

    outs, errs = stream.communicate()
    print(outs, errs)
    print("THERIAK RUN COMPLETED !!!")

    therin_file = 'THERIN.txt'
    therend_file = 'THEREND.txt'
    thkout_file = 'thkout.out'
    thktab_file = 'thktab.tab'
    folder_path = os.path.join(working_dir, loop_id, folder_name)
    print("folder_path is:" + folder_path)
    os.makedirs(folder_path)

    #copy files therin
    original1 = os.path.join(working_dir, therin_file)
    target1 = os.path.join(folder_path, therin_file)
    shutil.copyfile(original1, target1)

    #copy files therend
    original2 = os.path.join(working_dir, therend_file)
    target2 = os.path.join(folder_path, therend_file)
    shutil.copyfile(original2, target2)

    #copy files thkout_file
    original3 = os.path.join(working_dir, thkout_file)
    target3 = os.path.join(folder_path, thkout_file)
    shutil.copyfile(original3, target3)

    #copy files thktab_file
    original4 = os.path.join(working_dir, thktab_file)
    target4 = os.path.join(folder_path, thktab_file)
    shutil.copyfile(original4, target4)

    print("output files copied ! here: " + folder_path)
    return target4, folder_path


def extract_phase(df_in):
    """extract the phase list and put them in ordered list"""
    head = list(df_in.columns)
    head_cl = []
    for col0 in head:
        head_cl.append(col0.replace(' ', ''))
    df_in.columns = head_cl
    itera = 0
    idx = []
    # only pick list of column with main phase (not endmembers)
    for col in head_cl:
        itera = itera + 1
        num_t = col.split('_')
        if len(num_t) == 3:
            pass
        elif len(num_t) == 2:
            if 'V_' in col:
                idx.append(col)
        else:
            pass
    # remove CA column and
    for col2 in idx:
        if 'CA' in col2:
            idx.remove(col2)
        else:
            pass
    phases = idx[2::]
    list_of_phases = [w.replace('V_', '').replace('[', '').replace(']', '') for w in phases]

    if 'melt' in phases:
        phases.remove('melt')
        phases.insert(0, 'melt')
    elif 'V_[h2oL1]' in phases:
        phases.remove('V_[h2oL1]')
        phases.insert(0, 'V_[h2oL1]')
    else:
        pass
    list_of_V = phases
    list_of_rho = [z.replace('V_', 'rho_') for z in phases]

    return list_of_phases, list_of_V, list_of_rho


def merge_col(df_w, list_original, list_to_merge, new_head):
    """merge col of a same phase"""
    df_w[new_head] = df_w[list_to_merge].sum(axis=1)
    list_updated = [e for e in list_original if e not in list_to_merge]
    if new_head == 'melt':
        list_updated.insert(0, new_head)
    elif list_updated[0] == "melt":
        list_updated.insert(1, new_head)
    else:
        list_updated.insert(0, new_head)
    return df_w, list_updated


def read_tkh_out(filename=r'C:\TheriakDominoWIN\TD-win-15.r187-20180826\thktab.tab'):
    """read output theriak file, in df and organize data for plotting"""
    df_o = pd.read_csv(filename)
    df_in = df_o.copy()
    list_of_phases, list_of_V, list_of_rho = extract_phase(df_in)
    print(list_of_V)

    temp = ':Temperature'
    pressure = ':Pressure'

    # hereafter , you'll find a simplification of mineral phases groups for cumulative phase plotting
    # this list of phases may be completed and can be modified on the fly
    # Epidote group : merge epidote group
    set_ep = set(['V_[ep]', 'V_[cz]'])
    int_ep = set_ep.intersection(list_of_V)
    df_newv, list_of_V = merge_col(df_in, list_of_V, int_ep, 'epidote')

    # merge quartz group
    set_qtz = set(['V_q'])
    int_qtz = set_qtz.intersection(list_of_V)
    df_newx, list_of_V = merge_col(df_newv, list_of_V, int_qtz, 'quartz')

    # merge feti min group
    set_sph = set(['V_sph'])
    int_sph = set_sph.intersection(list_of_V)
    df_newy, list_of_V = merge_col(df_newx, list_of_V, int_sph, 'titanite')

    # merge micas group
    set_feld = set(['V_[phl]', 'V_[ann1]', 'V_[mu]', 'V_[tbi]'])
    int_feld = set_feld.intersection(list_of_V)
    df_new2, list_of_V = merge_col(df_newy, list_of_V, int_feld, 'micas')

    # merge feld group
    set_feld = set(['V_[abh1]', 'V_[san2]', 'V_[an2]', 'V_ab'])
    int_feld = set_feld.intersection(list_of_V)
    df_new2, list_of_V = merge_col(df_in, list_of_V, int_feld, 'feld')

    # merge ilm group
    set_ilm = set(['V_[dilm]', 'V_[oilm]'])
    int_ilm = set_ilm.intersection(list_of_V)
    df_new3, list_of_V = merge_col(df_new2, list_of_V, int_ilm, 'ilm')

    # merge opx group
    set_opx = set(['V_[en]', 'V_[fm]'])
    int_opx = set_opx.intersection(list_of_V)
    df_new4, list_of_V = merge_col(df_new3, list_of_V, int_opx, 'opx')

    # merge cpx group
    set_cpx = set(['V_[di1]', 'V_[acm1]', 'V_[ocats]'])
    int_cpx = set_cpx.intersection(list_of_V)
    df_new5, list_of_V = merge_col(df_new4, list_of_V, int_cpx, 'cpx')

    # merge amp group
    set_amp = set(['V_[parg1]', 'V_[b]', 'V_[ts1]', 'V_[gl1]', 'V_[tr]', 'V_[mrb]'])
    int_amp = set_amp.intersection(list_of_V)
    df_new6, list_of_V = merge_col(df_new5, list_of_V, int_amp, 'amp')

    # merge grt group
    set_grt = set(['V_[spss]', 'V_[gr]', 'V_[alm]', 'V_[py]', 'V_[py].1', 'V_[kho]'])
    int_grt = set_grt.intersection(list_of_V)
    df_new7, list_of_V = merge_col(df_new6, list_of_V, int_grt, 'grt')

    # merge melt group
    set_melt = set(['V_[abL1]', 'V_[h2oL1]', 'V_[kspL1]'])
    int_melt = set_melt.intersection(list_of_V)
    df_new8, list_of_V = merge_col(df_new7, list_of_V, int_melt, 'melt')

    return df_new8, temp, pressure, list_of_V


def plot_modes(df_imp, x_val, fixed_val, list_phases, plot_y2='no', id_title=" ", output_dir='new_folder'):
    """Create and display Tx cumulative plot for modelled solid and liquid phases"""
    df_plot = df_imp[list_phases]

    df_imp['vol_sum'] = df_imp[list_phases].sum()
    print(list_phases)
    df_plot = df_plot.divide(df_plot.sum(axis=1), axis=0)
    df_plot[x_val] = df_imp[x_val]
    ax = df_plot.plot.area(x=x_val, stacked=True, legend=False, title="Phase vol. proportions " + str(id_title) + " Ma", color=('firebrick', 'orange', 'mediumseagreen', 'royalblue', 'grey', 'black', 'gold', 'bisque', 'purple', 'olive', "teal", "black", "sienna", "darkgrey"))

    ax.set_ylabel('[vol. %]')
    ax.set_xlabel('Temperature [°C]')

    ax.set_ylim([0, 1])
    min_temp = 600
    max_temp = 1200
    ax.set_xlim([min_temp, max_temp])
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

    run_ide = os.path.basename(os.path.dirname(output_dir))
    plt.savefig(os.path.join(output_dir, "tx_plot_" + id_title + ".png"))
    plt.savefig(os.path.join(output_dir, "tx_plot_" + id_title + ".eps"))

    plt.savefig(os.path.join(theriak_wd, run_ide, 'all_plots', "tx_plot_" + id_title + ".png"))

    # show result plot for 2 sec then close and keep on modelling
    plt.show(block=False)
    plt.pause(0.5)
    plt.close()


def plot_no_melt_mode(df_imp, x_val, fixed_val, list_phases2, id_title=" ", output_dir='new_folder'):
    """Create and display Tx cumulative plot but just solid phases"""
    list_phases2.remove('melt')
    print(list_phases2)
    df_plot = df_imp[list_phases2]
    print(df_plot)
    df_imp['vol_sum'] = df_imp[list_phases2].sum()
    df_plot = df_plot.divide(df_plot.sum(axis=1), axis=0)
    df_plot[x_val] = df_imp[x_val]
    ax2 = df_plot.plot.area(x=x_val, stacked=True, legend=False, title="Phase vol. proportions " + str(id_title) + " Ma", color=('orange', 'mediumseagreen', 'royalblue', 'grey', 'black', 'gold', 'bisque', 'purple', 'olive', "teal", "black", "sienna", "darkgrey"))

    ax2.set_ylabel('[vol. %]')
    ax2.set_xlabel('Temperature [°C]')

    ax2.set_ylim([0, 1])
    #other option here is to set the min and max temp manually
    min_temp = 600
    max_temp = 1200
    ax2.set_xlim([min_temp, max_temp])
    ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

    run_ide = os.path.basename(os.path.dirname(output_dir))
    plt.savefig(os.path.join(output_dir, "tx_plot_no-melt_" + id_title + ".png"))
    plt.savefig(os.path.join(output_dir, "tx_plot_no-melt_" + id_title + ".eps"))

    plt.savefig(os.path.join(theriak_wd, run_ide, 'all_plots', "tx_plot_no-melt_" + id_title + ".png"))

    # show result plot for 2 sec then close and keep on
    plt.show(block=False)
    plt.pause(0.5)
    plt.close()



###########################################################################
# this is where users should modify parameters for treatment and plotting
###########################################################################
if __name__ == '__main__':
    # import xls file as df : each row is one composition to run (see example in repo)
    run_name = "Name_of_this_run"
    file_comp_in = r'your\data_inputs_compo.xlsx'

    # import xls file as df : each row is one composition to run (see example in repo)
    df_comp = import_data(file_comp_in, "input_serie")

    # Set theriak domino working directory
    theriak_wd = r'C:\TheriakDominoWIN\TD-win-15.r187-20180826'
    os.makedirs(os.path.join(theriak_wd, run_name, 'all_plots'))

    #create a folder and dir to transit with therin and therend files
    parent_dir = r"a\transit\folder\directory"
    today = datetime.datetime.now()
    # this allows them to have "unique" IDs
    date_time = today.strftime("%Y%m%d" + "_" + "%H%M%S")
    list_of_VH2O, all_V = [], []

    # here you may iterate for each row in xls files (age, various composition, etc.)
    for index, row in df_comp.iterrows():
        # those are the limit parameters to play with Pressure (P) start and end and Temperature (T) similarly, steps for modeling
        # we may add these parameters into the xls file too so we can iterate easily
        # parameters hereafter are for 500 T/P geotherm and is linearized between each point
        p_in, t_in = 10000, 600
        p_out, t_out, steps = 22000, 1200, 50

        # you can specify here the redox condition if known
        rat_Fe2O3FeO = 0.22

        # run theriak conv to convert oxid prop in mole and taking water content and redox conditions into account
        td_exp, name, comment = theriak_conv(row, rat_Fe2O3FeOt=rat_Fe2O3FeO)
        path = os.path.join(parent_dir, date_time, name)
        os.makedirs(path)
        run_id = str(row['id'])
        print(run_id)

        # create therin and therend files based on content from xls files
        therin_d, therend_d = therin_file(td_exp, comment, p_in, t_in, path, p_out, t_out, steps)

        # copy those files in the right dir so they can be red by theriak-domino app
        target_trin = os.path.join(theriak_wd, 'THERIN.txt')
        print(target_trin)
        try:
            os.remove(target_trin)
            print('it s copying')
            shutil.copyfile(therin_d, target_trin)
        except:
            print('it s copying')
            shutil.copyfile(therin_d, target_trin)

        target_trend = os.path.join(theriak_wd, 'THEREND.txt')
        print(target_trend)
        try:
            os.remove(target_trend)
            shutil.copyfile(therend_d, target_trend)
        except:
            shutil.copyfile(therend_d, target_trend)
        print('files copied in td dir, ready for runs !')

        # run the theriak app now and copy the outputs
        thk_tab_file, output_dir = run_theriak(file_out="THEREND.txt", name_out=run_id, folder_name=run_id, loop_id=run_name, source_trin=therin_d, source_trend=therend_d, working_dir=theriak_wd)




        df_tkh_tab, temp, pressure, list_of_V = read_tkh_out(filename=thk_tab_file)
        all_V.append(list_of_V)
        print(list_of_V)
        if 'V_H2O' in list_of_V:
            list_of_VH2O.append('V_H2O')
        else:
            list_of_VH2O.append('no')
        plot_modes(df_imp=df_tkh_tab, x_val=temp, fixed_val=pressure, list_phases=list_of_V, id_title=run_id, output_dir=output_dir)
        plot_no_melt_mode(df_imp=df_tkh_tab, x_val=temp, fixed_val=pressure, list_phases2=list_of_V, id_title=run_id, output_dir=output_dir)
    print(list_of_VH2O)
    print(all_V)
