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


import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import os
import argparse
import matplotlib as mpl
from matplotlib import cm
mpl.rcParams['pdf.fonttype'] = 42

# Find the path to an existing ontology document that was installed through the Brain Explorer interface
path = r'C:\Users\jenniferwh\AppData\Local\Allen Institute\Brain Explorer 2\Atlases\Allen Mouse Brain Common Coordinate Framework'
dat = pd.read_csv(os.path.join(path, 'ontology_v2.csv'))

dataset = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Denali_collaboration\dataset.csv')
dataset = dataset[dataset['hemisphere'] != 'Fail']
dataset = dataset[dataset['age_group'] == '8 mo']
print(len(dataset))

unionize_dat = pd.read_csv(r'C:\Users\jenniferwh\Dropbox (Allen Institute)\Denali_collaboration\plaque_unionizes.csv')
unionize_dat['count_per_mm_3'] = unionize_dat['sum_plaques']/unionize_dat['volume'] #to convert to plaques per mm^3
unionize_dat = unionize_dat[unionize_dat['hemisphere_id'] != 3]
unionize_dat.loc[(unionize_dat['image_series_id'].isin(dataset[dataset['hemisphere'] == 'R']['image_series_id'].unique())) &
        (unionize_dat['hemisphere_id'] == 2), 'include'] = 'yes'
unionize_dat.loc[(unionize_dat['image_series_id'].isin(dataset[dataset['hemisphere'] == 'L']['image_series_id'].unique())) &
        (unionize_dat['hemisphere_id'] == 1), 'include'] = 'yes'
unionize_dat = unionize_dat[unionize_dat['include'] == 'yes']
print(len(unionize_dat['image_series_id'].unique()))

mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()

def get_cmap(colormap=cm.gray, scale=1):
    median_vals = []
    for structure in unionize_dat['structure_id'].unique():
        med = int(unionize_dat[unionize_dat['structure_id'] == structure][
                                            'count_per_mm_3'].median())
        median_vals.append(med)
    structure_vals = dict(zip(unionize_dat['structure_id'].unique(), 
                              median_vals))
    rgb_vals = structure_vals.copy()
    for key in structure_vals.keys():
        rgb_vals[key] = tuple([255*i for i in colormap(structure_vals[key]*scale)[:3]])
    rgb_vals[0] = (0, 0, 0)
    return rgb_vals

def main():
    rgb_vals = get_cmap(cm.gray_r, args.scale)
    for key in rgb_vals.keys():
        dat.loc[dat['database_id'] == key, 'red'] = int(rgb_vals[key][0])
        dat.loc[dat['database_id'] == key, 'green'] = int(rgb_vals[key][1])
        dat.loc[dat['database_id'] == key, 'blue'] =int(rgb_vals[key][2])
    dat.loc[dat['abbreviation'] == 'root', 'parent'] = 0.
    dat['parent'] = [int(value) for value in dat['parent']]
    
    dat.to_csv(os.path.join(path, 'ontology_v2.csv'), index=False)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('scale', type = float, help = 'scale values by this factor')
    args = parser.parse_args()

    main()