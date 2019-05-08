# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:02 2018

@author: JENNIFERWH
Creates and saves a series of coronal image sections from the 3D heatmap of 
plaque density values in the plaque_densities_per_structure.csv file.
This image series can be converted to a movie as in 
https://www.biorxiv.org/content/early/2018/08/18/395236 Movie 1.

Parameters
----------
plaque_movie_params.json file that includes the mouse lines to map (APP/PS1,
Tg2576, hAPP-J20), scale for the heatmap for each mouse line, age to use for
mapping each mouse line, and an output path for saving images.

Returns
-------
Saves a series of .png files containing the images at 25 um intervals.
"""

import os
import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from matplotlib import cm
import json
import matplotlib.pyplot as plt

path = r'/Users/jenniferwh/Dropbox/Diamond_collaboration'
dat = pd.read_csv('/Users/jenniferwh/Dropbox/Diamond_collaboration/transposed_data.csv')
dat = pd.melt(dat, id_vars = ['SampleID', 'Condition'], value_name = 'density', var_name='acronym')
mcc = MouseConnectivityCache(manifest_file = '..connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
dat['acronym'] = [acronym[1:-1] for acronym in dat['acronym'].values]
for structure in dat['acronym']:
    if structure == 'OutsideBrainBoundaries':
        dat.loc[dat['acronym'] == structure, 'structure_id'] = int(0)
    else:
        dat.loc[dat['acronym'] == structure, 'structure_id'] = ia_map[structure]
    
dat.loc[dat['Condition'].isin(['12monthTG', '12monthTG.1', '12monthTG.2', '12monthTG.3', 
        '11monthTG', '11monthTG.1', '11monthTG.2']), 'group'] = '11-12mo'
dat.loc[dat['Condition'].isin(['9monthTG', '9monthTG.1', '9monthTG.2', '8monthTG', '8monthTG.1', 
                        '8monthTG.2']), 'group'] = '8-9mo'
dat.loc[dat['Condition'].isin(['3monthTG', '3monthTG.1', '3monthTG.2']), 'group'] = '3mo'

def get_mean_value_per_structure(group, structure_ids):
    means = []
    structs = []
    isids = dat[(dat['group'] == group)]['SampleID'].values
    for structure_id in structure_ids:
        if len(dat[(dat['structure_id'] == structure_id) & 
                       (dat['SampleID'].isin(isids))]) > 0:
            str_mean = dat[(dat['structure_id'] == structure_id) & 
                           (dat['SampleID'].isin(isids))]['density'].mean()
            means.append(str_mean)
            structs.append(structure_id)
    # Not all structures are represented in this dataset. 
    # Color missing structures by their parent
    structures = pd.DataFrame(st.get_structures_by_set_id([184527634]))
    missing_strs = [structure for structure in structures['id'].values if 
                    structure not in structs]
    print(len(missing_strs))
    while len(missing_strs) > 0:
        for structure in missing_strs:
            ancestors = structures[structures['id'] == structure]['structure_id_path'].values[0]
            for n in range(2, len(ancestors)+1):
                ancestor = ancestors[-n:][0]
                if ancestor in dat['structure_id'].unique():
                    str_mean = dat[(dat['structure_id'] == ancestor) & 
                               (dat['SampleID'].isin(isids))]['density'].mean()
                    means.append(str_mean)
                    structs.append(structure)
                    missing_strs.remove(structure)
                    break;
    structuredat = dict(zip(structs, means))
    return structuredat, len(isids)

def get_cmap(group, scale=1):
    structure_vals, n = get_mean_value_per_structure(group, dat['structure_id'].unique())
    rgb_vals = structure_vals.copy()
    for key in structure_vals:
        rgb_vals[key] = tuple([255*i for i in cm.hot(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (0, 0, 0)
    return rgb_vals, n

def main(params):
    reference_space =  mcc.get_reference_space()
    structure_values = []
    cmaps = {}
    for ix, group in enumerate(params['groups']):
        structure_vals, n = get_mean_value_per_structure(group, 
                                                     dat.structure_id.unique())
        rgb_vals, n = get_cmap(params['groups'][ix], params['scale'][ix])
        structure_values.append(structure_vals)
        cmaps[group] = rgb_vals
    imax = [index/8 for index in reference_space.annotation.shape]
    for ix in range(int(imax[0])):
        image = []
        for group in ['3mo', '8-9mo', '11-12mo']:
            new_im = reference_space.get_slice_image(0, ix*25*8, cmaps[group])
            image.append(new_im)
        composite_image = np.concatenate(np.array(image), axis = 1)
        plt.imshow(composite_image)
        plt.axis('off')
        plt.savefig(os.path.join(params['savepath'], 
                             '{0}.png'.format(ix)), 
            facecolor=rgb_vals[0],
            bbox_inches='tight', 
            pad_inches=0.3, 
            format='png', 
            dpi=500)
    structure_vals, n = get_mean_value_per_structure(params['group'], 
                                                     dat.structure_id.unique())
    rgb_vals, n = get_cmap(params['group'], params['scale'])
    image = [0,0,0]
    image[0] = reference_space.get_slice_image(0, params['xcoord'], rgb_vals) #posterior
    image[1] = np.flip(np.rot90(reference_space.get_slice_image(2, params['zcoord'], rgb_vals)), 0) #right
    image[2] = np.rot90(reference_space.get_slice_image(1, params['ycoord'], rgb_vals))   #inferior
    fig = plt.figure(figsize=(8, 3), facecolor=rgb_vals[0])
    columns = 3
    rows = 1
    for i in range(columns*rows):
        fig.add_subplot(rows, columns, i+1)
        f = plt.imshow(image[i], cmap = cm.hot)
        plt.axis('off')
    if not os.path.exists(params['savepath']):
        os.mkdir(params['savepath'])
    mouse_line_fname = params['group']
    plt.savefig(os.path.join(params['savepath'], 
                             '{0}_map_{1}-{2}-{3}.png'.format(mouse_line_fname,
                             params['xcoord'],
                             params['zcoord'],
                             params['ycoord'])), 
            facecolor=fig.get_facecolor(),
            bbox_inches='tight', 
            pad_inches=0.3, 
            format='png', 
            dpi=300)
    
if __name__ == '__main__':

    with open('movie_params.json', 'r') as data_file:    
        parameters = json.load(data_file)

    main(parameters)