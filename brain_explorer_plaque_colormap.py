# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:08:12 2018

@author: jenniferwh

Changes the colors in the Allen Brain Explorer App to a customized heatmap
corresponding to plaque density. Can be used to visualize the data from 
https://www.biorxiv.org/content/early/2018/08/18/395236. Data from the paper
is included as plaque_densities_per_structure.csv.
Before using this code, download Brain Explorer from 
http://mouse.brain-map.org/static/brainexplorer.
Run Brain Explorer and download the Allen Mouse Common Coordinate Framework
atlas using the prompts in the Brain Explorer program. The ontology_v2.csv file
will be installed in C:\Users\*USERNAME*\AppData\Local\Allen Institute\Brain Explorer 2\Atlases\Allen Mouse Brain Common Coordinate Framework
on a PC, or /Users/*USERNAME*/Library/Application Support/Brain Explorer 2/Atlases/Allen Mouse Brain Common Coordinate Framework
on a mac. This code modifies the ontology_v2.csv file installed by Brain Explorer 
to customize the colormap. See http://help.brain-map.org/display/mousebrain/Brain+Explorer 
forcmore information about Brain Explorer installation.

Parameters
----------
mouse_line: string. Must exactly match column name in plaque_densities_per_structure.csv.
age: int
scale: factor by which to scale density values (1 for no scaling)

Returns
-------
Saves 'ontology_v2.csv' in the user-specified location. This file is used by 
the Brain Explorer App to assign colors to CCF structures.
"""

import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import os
import argparse
import matplotlib as mpl
from matplotlib import cm
mpl.rcParams['pdf.fonttype'] = 42

# Find the path to an existing ontology document that was installed through the Brain Explorer interface
path = r'/Users/jenniferwh/Library/Application Support/Brain Explorer 2/Atlases/Allen Mouse Brain Common Coordinate Framework'
dat = pd.read_csv(os.path.join(path, 'ontology_v2.csv'))

unionize_dat = pd.read_csv('plaque_densities_per_structure.csv')
unionize_dat = unionize_dat[unionize_dat['control'] == False]
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
unionize_dat['structure_id'] = [ia_map[structure] for structure in unionize_dat['structure_acronym']]

def get_mean_value_per_structure(group, age, structure_ids):
    means = []
    if len(unionize_dat[(unionize_dat['mouse_line'] == group) & (unionize_dat['age_group'] == age)]) > 0:
        isids = unionize_dat[(unionize_dat['mouse_line'] == group) & 
                             (unionize_dat['age_group'] == age)]['image_series_id'].values
        for structure_id in structure_ids:
            density = np.mean(unionize_dat[(unionize_dat['structure_id'] == structure_id) & 
                                  (unionize_dat['image_series_id'].isin(isids))]['plaque_density'])
            means.append(density*100) #To convert from fraction to percent           
        structuredat = dict(zip(structure_ids, means))
    else:
        structuredat = dict(zip(structure_ids, np.zeros(len(structure_ids))))
        isids = []
    return structuredat, len(isids)

def get_cmap(group, age, colormap=cm.gray, scale=1):
    structure_vals, n = get_mean_value_per_structure(group, age, unionize_dat['structure_id'].unique())
    rgb_vals = structure_vals.copy()
    for key in structure_vals:
        rgb_vals[key] = tuple([255*i for i in colormap(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (0, 0, 0)
    return rgb_vals, n

def main():
    rgb_vals, n = get_cmap(args.mouse_line, str(args.age)+' mo', cm.inferno, args.scale)
    for key in rgb_vals.keys():
        dat.loc[dat['database_id'] == key, 'red'] = int(rgb_vals[key][0])
        dat.loc[dat['database_id'] == key, 'green'] = int(rgb_vals[key][1])
        dat.loc[dat['database_id'] == key, 'blue'] =int(rgb_vals[key][2])
    dat.loc[dat['abbreviation'] == 'root', 'parent'] = 0.
    dat['parent'] = [int(value) for value in dat['parent']]
    
    dat.to_csv(os.path.join(path, 'ontology_v2.csv'), index=False)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('mouse_line', type = str, help='mouse line (APP/PS1, hAPP-J20, or Tg2576)')
    parser.add_argument('age', type = int, help='age (5, 7, 9, 13, or 19)')
    parser.add_argument('scale', type = float, help = 'scale values by this factor')
    args = parser.parse_args()

    main()