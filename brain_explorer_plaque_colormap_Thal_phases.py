# -*- coding: utf-8 -*-
# Created on Tue Sep 18 12:08:12 2018

# @author: jenniferwh

# Changes the colors in the Allen Brain Explorer App to a customized heatmap
# corresponding to plaque density. Can be used to visualize the data from 
# https://www.biorxiv.org/content/early/2018/08/18/395236. Data from the paper
# is included as plaque_densities_per_structure.csv.
# Before using this code, download Brain Explorer from 
# http://mouse.brain-map.org/static/brainexplorer.
# Run Brain Explorer and download the Allen Mouse Common Coordinate Framework
# atlas using the prompts in the Brain Explorer program. The ontology_v2.csv file
# will be installed in C:\Users\*USERNAME*\AppData\Local\Allen Institute\Brain Explorer 2\Atlases\Allen Mouse Brain Common Coordinate Framework
# on a PC, or /Users/*USERNAME*/Library/Application Support/Brain Explorer 2/Atlases/Allen Mouse Brain Common Coordinate Framework
# on a mac. This code modifies the ontology_v2.csv file installed by Brain Explorer 
# to customize the colormap. See http://help.brain-map.org/display/mousebrain/Brain+Explorer 
# forcmore information about Brain Explorer installation.

# Parameters
# ----------
# mouse_line: string. Must exactly match column name in plaque_densities_per_structure.csv.
# age: int
# scale: factor by which to scale density values (1 for no scaling)

# Returns
# -------
# Saves 'ontology_v2.csv' in the user-specified location. This file is used by 
# the Brain Explorer App to assign colors to CCF structures.

import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import os
import argparse
import matplotlib as mpl
from matplotlib import cm
from anatomy.anatomy_api import AnatomyApi
mpl.rcParams['pdf.fonttype'] = 42

def main():
    # Find the path to an existing ontology document that was installed through the Brain Explorer interface
    path = r'C:\Users\jenniferwh\AppData\Local\Allen Institute\Brain Explorer 2\Atlases\Allen Mouse Brain Common Coordinate Framework'
    dat = pd.read_csv(os.path.join(path, 'ontology_v2.csv'))
    aapi = AnatomyApi()
    mcc = MouseConnectivityCache(manifest_file = '/connectivity/mouse_connectivity_manifest.json')
    st = mcc.get_structure_tree()
    ia_map = st.get_id_acronym_map()
    ss = aapi.get_summary_structure_data('id')
    
    if args.Thal_phase == 1:
        thal_structures = ['Isocortex']
    if args.Thal_phase == 2:
        thal_structures = ['CA1', 'ENTl', 'ENTm']
    if args.Thal_phase == 3:
        thal_structures = ['ACA', 'RSP', 'LA', 'BLA', 'BMA', 'PA', 'CEA',
                     'MEA', 'COA', 'DG', 'PRE', 'POST', 'TH', 'CP', 'ACB', 'HY', 'SI', 'MA', 'NDB']
    if args.Thal_phase == 4:
        thal_structures = ['PAG', 
                     'SCs', 'SCm', 'RN', 'IO', 'SNr', 'SNc']
    if args.Thal_phase == 5:
        thal_structures = ['GRN', 'IRN', 'PARN', 'CBX', 'PRNr', 'PRNc', 
                     'RAmb', 'LC', 'PB', 'TRN', 'DTN', 'PG']
    thal_structure_ids = [ia_map[structure] for structure in thal_structures]
    '''
    for structure in thal_structure_ids:
        if structure not in ss:
            children = st.descendant_ids([structure])[0]
            child_ss = [structure_id for structure_id in children if structure_id in ss]
            ss += [structure]
            thal_structure_ids+=child_ss
    '''
    colormap = cm.Reds
    all_structures = dat['abbreviation'].unique()
    all_structures = [structure for structure in all_structures if structure in ia_map.keys()]
    all_structures.remove('root')
    all_structure_ids = [ia_map[structure] for structure in all_structures]
    structure_vals = dict()
    for structure in all_structure_ids:
        if structure in thal_structure_ids:
            structure_vals[structure] = 0.9
        else:
            structure_vals[structure] = 0
    rgb_vals = structure_vals.copy()
    for key in all_structure_ids:
        rgb_vals[key] = tuple([255*i for i in colormap(structure_vals[key])[:3]])
        dat.loc[dat['database_id'] == key, 'red'] = int(rgb_vals[key][0])
        dat.loc[dat['database_id'] == key, 'green'] = int(rgb_vals[key][1])
        dat.loc[dat['database_id'] == key, 'blue'] =int(rgb_vals[key][2])
    rgb_vals[0] = (0, 0, 0)
    dat.loc[dat['abbreviation'] == 'root', 'parent'] = 0.
    dat['parent'] = [int(value) for value in dat['parent']]
    
    dat.to_csv(os.path.join(path, 'ontology_v2.csv'), index=False)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('Thal_phase', type = int, help='Thal Phase (integer')
    args = parser.parse_args()

    main()