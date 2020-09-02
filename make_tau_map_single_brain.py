# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:02 2018

@author: JENNIFERWH
Creates planar image sections from the 3D heatmap of plaque density per structure
as described in https://www.biorxiv.org/content/early/2018/08/18/395236. 
Data from the paper is included in the plaque_densities_per_structure.csv file.

Parameters
----------
image_params.json file that includes the mouse line to map (APP/PS1,
hAPP-J20, or Tg2576); age to map (5, 7, 9, 13, or 19 months), 
scale for the heatmap; x, y, and z coordinates for the
coronal, horizontal, and sagittal sections, respectively; and an output path
for saving images. 

Returns
-------
Saves a .png file containing the images at the requested x,y,z coordinates as
in Figures 9 and 10 in Whitesell et al., 2018 (accepted).
"""

import os
import numpy as np
import pandas as pd
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from matplotlib import cm
import operator
import json
import matplotlib.pyplot as plt

path = r'/Users/jenniferwh/Dropbox/Diamond_collaboration'
dat = pd.read_csv('/Users/jenniferwh/Dropbox (Allen Institute)/Diamond_collaboration/data_files/thresholded_data_LR.csv')
dat = pd.melt(dat, id_vars = ['SampleID', 'Condition'], value_name = 'density', var_name='acronym')
mcc = MouseConnectivityCache(manifest_file = '..connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
dat['acronym'] = [acronym[:-2] for acronym in dat['acronym'].values]
for structure in dat['acronym'].unique():
    dat.loc[dat['acronym'] == structure, 'structure_id'] = ia_map[structure]
    
dat.loc[dat['SampleID'] == 371, 'group'] = '371'
dat.loc[dat['group'].isnull(), 'group'] = 'other'

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
            structs.append(int(structure_id))
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
        rgb_vals[key] = tuple([255*i for i in cm.BuPu(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (0, 0, 0)
    print(rgb_vals)
    return rgb_vals, n

def main(params):
    reference_space =  mcc.get_reference_space()
    structure_vals, n = get_mean_value_per_structure(params['group'], 
                                                     dat.structure_id.unique())
    rgb_vals, n = get_cmap(params['group'], params['scale'])
    image = [0,0,0]
    print(params['xcoord'])
    image[0] = reference_space.get_slice_image(0, params['xcoord'], rgb_vals) #posterior
    image[1] = np.flip(np.rot90(reference_space.get_slice_image(2, params['zcoord'], rgb_vals)), 0) #right
    image[2] = np.rot90(reference_space.get_slice_image(1, params['ycoord'], rgb_vals))   #inferior
    fig = plt.figure(figsize=(8, 3), facecolor=rgb_vals[0])
    columns = 3
    rows = 1
    for i in range(columns*rows):
        fig.add_subplot(rows, columns, i+1)
        f = plt.imshow(image[i], cmap = cm.BuPu)
        plt.axis('off')
    cbar = fig.colorbar(f, fraction=0.046)
    cbar.set_label('pTau Probability (per mm$\mathregular{^{3}}$)', 
                   rotation=90, color = 'w', fontsize = 12)
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w', fontsize = 10)
    
    maxval = max(structure_vals.items(), key=operator.itemgetter(1))[1]
    cbar.ax.set_yticklabels([0, 
                             np.round(maxval*.2, 2), 
                             np.round(maxval*.4, 2),
                             np.round(maxval*.6, 2),
                             np.round(maxval*.8, 2),
                             np.round(maxval, 2)])
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

    with open('image_params_single_brain.json', 'r') as data_file:    
        parameters = json.load(data_file)

    main(parameters)