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

dat = pd.read_csv('plaque_densities_per_structure.csv')
dat = dat[dat['control'] == False]
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
dat['structure_id'] = [ia_map[structure] for structure in dat['structure_acronym']]

def get_mean_value_per_structure(group, age, structure_ids):
    means = []
    isids = dat[(dat['mouse_line'] == group) & (dat['age_group'] == age)]['image_series_id'].values
    for structure_id in structure_ids:
        str_mean = np.mean(dat[(dat['structure_id'] == structure_id) & 
                                 (dat['image_series_id'].isin(isids))]
                             ['plaque_density'])
        means.append(str_mean*100) #To convert from fraction to percent
    structuredat = dict(zip(structure_ids, means))
    return structuredat, len(isids)

def get_cmap(group, age, scale=1):
    structure_vals, n = get_mean_value_per_structure(group, age, dat['structure_id'].unique())
    rgb_vals = structure_vals.copy()
    for key in structure_vals:
        rgb_vals[key] = tuple([255*i for i in cm.hot(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (0, 0, 0)
    return rgb_vals, n

def main(params):
    reference_space =  mcc.get_reference_space()
    
    structure_vals, n = get_mean_value_per_structure(params['mouse_line'], 
                                                     params['age'], 
                                                     dat.structure_id.unique())
    rgb_vals, n = get_cmap(params['mouse_line'], params['age'], params['scale'])
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
    cbar = fig.colorbar(f, fraction=0.046)
    cbar.set_label('Plaque % Volume', rotation=90, color = 'w', fontsize = 12)
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w', fontsize = 10)
    
    maxval = max(structure_vals.iteritems(), key=operator.itemgetter(1))[1]
    cbar.ax.set_yticklabels([0, 
                             np.round(maxval*.125, 2), 
                             np.round(maxval*.25, 2), 
                             np.round(maxval*.375, 2),
                             np.round(maxval*.5, 2), 
                             np.round(maxval*.625, 2),
                             np.round(maxval*.75, 2),
                             np.round(maxval*.875, 2),
                             np.round(maxval, 2)])
    if not os.path.exists(params['savepath']):
        os.mkdir(params['savepath'])
    if params['mouse_line'] == 'APP/PS1':
        mouse_line_fname = 'APP_PS1'
    else:
        mouse_line_fname = params['mouse_line']
    plt.savefig(os.path.join(params['savepath'], 
                             '{0}_{1}_map_{2}-{3}-{4}.png'.format(mouse_line_fname,
                             params['age'],
                             params['xcoord'],
                             params['zcoord'],
                             params['ycoord'])), 
            facecolor=fig.get_facecolor(),
            bbox_inches='tight', 
            pad_inches=0.3, 
            format='png', 
            dpi=1000)
    
if __name__ == '__main__':

    with open('image_params.json', 'r') as data_file:    
        parameters = json.load(data_file)

    main(parameters)