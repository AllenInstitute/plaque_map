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
import json
import argparse
import matplotlib as mpl
from matplotlib import cm
mpl.rcParams['pdf.fonttype'] = 42

def main():
    # Find the path to an existing ontology document that was installed through the Brain Explorer interface
    path = r'C:\Users\jenniferwh\AppData\Local\Allen Institute\Brain Explorer 2\Atlases\Allen Mouse Brain Common Coordinate Framework'
    dat = pd.read_csv(os.path.join(path, 'ontology_v2.csv'))

    mcc = MouseConnectivityCache(manifest_file = '/connectivity/mouse_connectivity_manifest.json')
    st = mcc.get_structure_tree()
    ia_map = st.get_id_acronym_map()
    
    jsonfile = r'/Users/jenniferwh/Dropbox (Allen Institute)/Diamond_collaboration/data_files/consensus_structures.json'
    with open(jsonfile, 'r') as file:
        consensus_structures = json.load(file)
    age = str(args.age)+'monthTG'
    print(age)
    
    tau_structures = consensus_structures[age]
    tau_structure_ids = [ia_map[structure] for structure in tau_structures]
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
        if structure in tau_structure_ids:
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
    parser.add_argument('age', type = int, help='age group')
    args = parser.parse_args()

    main()