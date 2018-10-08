# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:02 2018

@author: JENNIFERWH
Creates planar image sections from the 3D heatmap of false positive values
showing where segmentation artifacts were observed in the plaque map dataset.
Heatmap is overlaid on the template brain to assist in visualizing structures.
See Figure 2k,l in https://www.biorxiv.org/content/early/2018/08/18/395236.
Data from the paper is included in the plaque_densities_per_structure.csv file.

Parameters
----------
false_positive_params.json file that includes the mouse line to map (APP/PS1,
hAPP-J20, or both); scale for the heatmap; x, y, and z coordinates for the
coronal, horizontal, and sagittal sections, respectively; and an output path
for saving images. 

Returns
-------
Saves a .png file containing the images at the requested x,y,z coordinates as
in Figure 2 k,l in Whitesell et al., 2018 (accepted).
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
dat = dat[dat['control'] == True]
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
dat['structure_id'] = [ia_map[structure] for structure in dat['structure_acronym']]

def get_mean_value_per_structure(group, structure_ids):
    means = []
    if group == 'all':
        isids = dat['image_series_id'].values
    else:
        isids = dat[dat['mouse_line'] == group]['image_series_id'].values
    for structure_id in structure_ids:
        str_mean = np.mean(dat[(dat['structure_id'] == structure_id) & 
                                 (dat['image_series_id'].isin(isids))]['plaque_density'])
        means.append(str_mean*100) #To convert from fraction to percent
    structuredat = dict(zip(structure_ids, means))
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
    template = mcc.get_template_volume()[0]
    if params['mouse_line'] in ['APP/PS1', 'hAPP-J20']:
        mouse_line = params['mouse_line']
    else:
        mouse_line = 'all'
    structure_vals, n = get_mean_value_per_structure(mouse_line, 
                                                     dat.structure_id.unique())
    rgb_vals, n = get_cmap(mouse_line, params['scale'])
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
        if i == 0:
            template_index = int(params['xcoord']/25)
            g = plt.imshow(template[template_index], alpha = 0.5, cmap = 'gray')
        if i == 1:
            template_index = int(params['zcoord']/25)
            g = plt.imshow(np.flip(np.rot90(template[:,:,template_index]), 0), 
                           alpha = 0.5, cmap = 'gray')
        if i == 2:
            template_index = int(params['ycoord']/25)
            g = plt.imshow(np.rot90(template[:,template_index,:]), alpha = 0.5, cmap = 'gray')
        plt.axis('off')
    cbar = fig.colorbar(f, fraction=0.046)
    cbar.set_label('False Positive % Volume', rotation=90, color = 'w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    
    maxval = max(structure_vals.iteritems(), key=operator.itemgetter(1))[1]/params['scale']
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
    if mouse_line == 'APP/PS1':
        mouse_line_fname = 'APP_PS1'
    else:
        mouse_line_fname = mouse_line
    plt.savefig(os.path.join(params['savepath'], 
                             'False_positive_{0}_map_{1}-{2}-{3}.png'.format(
                                     mouse_line_fname,
                                     params['xcoord'],
                                     params['zcoord'],
                                     params['ycoord'])), 
            facecolor=fig.get_facecolor(),
            bbox_inches='tight', 
            pad_inches=0.3, 
            format='png', 
            dpi=1000)
    
if __name__ == '__main__':

    with open('false_positive_params.json', 'r') as data_file:    
        parameters = json.load(data_file)

    main(parameters)