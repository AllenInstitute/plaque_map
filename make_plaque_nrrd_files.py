# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:02 2018

@author: JENNIFERWH

Produce nrrd file containing Allen CCFv3 structures color-coded by plaque density
as described in https://www.biorxiv.org/content/early/2018/08/18/395236. 
Data from the paper is included in the plaque_densities_per_structure.csv file.
Separate nrrd files are created for each mouse line - age group combination.

Parameters
----------
nrrd_params.json file that includes the mouse lines to map (APP/PS1,
hAPP-J20, Tg2576), ages to map (5, 7, 9, 13, or 19 months), and an output path
for saving files. 

Returns
-------
Saves a nearly raw raster data (nrrd) file for each mouse line - age group 
combination.

"""
import os
import numpy as np
import pandas as pd
import json
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from matplotlib import cm
import nrrd

dat = pd.read_csv('plaque_densities_per_structure.csv')
mcc = MouseConnectivityCache(manifest_file = '../connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()
dat['structure_id'] = [ia_map[structure] for structure in dat['structure_acronym']]

c_dat = dat[dat['control'] == True]
dat = dat[dat['control'] == False]

def get_mean_value_per_structure(group, age, structure_ids):
    means = []
    if len(dat[(dat['mouse_line'] == group) & (dat['age_group'] == age)]) > 0:
        isids = dat[(dat['mouse_line'] == group) & (dat['age_group'] == age)]['image_series_id'].values
        for structure_id in structure_ids:
            density = np.mean(dat[(dat['structure_id'] == structure_id) & 
                                  (dat['image_series_id'].isin(isids))]['plaque_density'])
            means.append(density*100) #To convert from fraction to percent           
        structuredat = dict(zip(structure_ids, means))
    else:
        structuredat = dict(zip(structure_ids, np.zeros(len(structure_ids))))
        isids = []
    return structuredat, len(isids)

def get_cmap(group, age, colormap=cm.gray, scale=1):
    structure_vals, n = get_mean_value_per_structure(group, age, dat['structure_id'].unique())
    rgb_vals = structure_vals.copy()
    for key in structure_vals:
        rgb_vals[key] = tuple([255*i for i in colormap(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (0, 0, 0)
    return rgb_vals, n

def get_control_cmap(colormap=cm.gray, scale=1):
    means = []
    for structure_id in c_dat['structure_id'].unique():
        density = np.mean(c_dat[c_dat['structure_id'] == structure_id]['plaque_density'])
        means.append(density*100) #To convert from fraction to percent 
    structure_vals = dict(zip(c_dat['structure_id'].unique(), means))
    rgb_vals = structure_vals.copy()
    for key in structure_vals:
        rgb_vals[key] = tuple([255*i for i in colormap(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (0, 0, 0)
    return rgb_vals, len(c_dat['image_series_id'].unique())

def main(params):
    age_groups = []
    for age in params['ages']:
        age_groups.append(str(age) + ' mo')
        
    reference_space =  mcc.get_reference_space()
    
    for group in age_groups:
        for mouse_line in params['mouse_lines']:
            rgb_vals, n = get_cmap(mouse_line, group)
            image = []
            for ix in range(reference_space.annotation.shape[0]):
                new_im = reference_space.get_slice_image(0, ix*25, rgb_vals) #position in microns, 25 um resolution
                image.append(new_im)
            composite_image = np.array(image)
            if mouse_line == 'APP/PS1':
                mouse_line_fname = 'APP_PS1'
            else:
                mouse_line_fname = mouse_line
            if len(np.unique(composite_image)) > 1:
                fname = os.path.join(params['savepath'], 
                        '{0} {1}.nrrd'.format(mouse_line_fname, group))
                nrrd.write(fname, composite_image[:,:,:,0]) #three channels are equal for gray cmap 
    
    if len(c_dat) > 0:
        rgb_vals, n = get_control_cmap()
        image = []
        for ix in range(reference_space.annotation.shape[0]):
            new_im = reference_space.get_slice_image(0, ix*25, rgb_vals) #position in microns, 25 um resolution
            image.append(new_im)
        composite_image = np.array(image)
        nrrd.write(os.path.join(params['savepath'], 'control.nrrd'),
        composite_image[:,:,:,0]) #three channels are equal for gray cmap
    
if __name__ == '__main__':

    with open('nrrd_params.json', 'r') as data_file:    
        parameters = json.load(data_file)

    main(parameters)