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
    structure_values = []
    cmaps = {}
    for ix, mouse_line in enumerate(params['mouse_line']):
        structure_vals, n = get_mean_value_per_structure(mouse_line, 
                                                     params['age'][ix], 
                                                     dat.structure_id.unique())
        rgb_vals, n = get_cmap(params['mouse_line'][ix], params['age'][ix], params['scale'][ix])
        structure_values.append(structure_vals)
        cmaps[mouse_line] = rgb_vals
    imax = [index/8 for index in reference_space.annotation.shape]
    for ix in range(imax[0]):
        image = []
        for mouse_line in ['APP/PS1', 'Tg2576', 'hAPP-J20']:
            new_im = reference_space.get_slice_image(0, ix*25*8, cmaps[mouse_line])
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
    
if __name__ == '__main__':

    with open('plaque_movie_params.json', 'r') as data_file:    
        parameters = json.load(data_file)

    main(parameters)