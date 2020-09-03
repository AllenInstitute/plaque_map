# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:43:02 2018

@author: JENNIFERWH
Creates planar image sections from the 3D heatmap of plaque density per structure
as described in https://www.biorxiv.org/content/early/2018/08/18/395236. 
Ldata from the paper is included in the plaque_densities_per_structure.csv file.

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

mcc = MouseConnectivityCache(manifest_file = '..connectivity/mouse_connectivity_manifest.json')
st = mcc.get_structure_tree()
ia_map = st.get_id_acronym_map()

path = r'/Users/jenniferwh/Dropbox/Diamond_collaboration'
dat = pd.read_csv('/Users/jenniferwh/Dropbox (Allen Institute)/Diamond_collaboration/data_files/thresholded_data_LR.csv')

meta = dat[['SampleID', 'Condition']]
Ldat = dat[[acronym for acronym in dat.columns if '-L' in acronym]]
Ldat = pd.concat([meta, Ldat])
Rdat = dat[[acronym for acronym in dat.columns if '-R' in acronym]]
Rdat = pd.concat([meta, Rdat])

Ldat = pd.melt(Ldat, id_vars = ['SampleID', 'Condition'], value_name = 'density', var_name='acronym')
Ldat['acronym'] = [acronym[:-2] for acronym in Ldat['acronym'].values]
Rdat = pd.melt(Rdat, id_vars = ['SampleID', 'Condition'], value_name = 'density', var_name='acronym')
Rdat['acronym'] = [acronym[:-2] for acronym in Rdat['acronym'].values]
for structure in Ldat['acronym'].unique():
    Ldat.loc[Ldat['acronym'] == structure, 'structure_id'] = ia_map[structure]
    Rdat.loc[Rdat['acronym'] == structure, 'structure_id'] = ia_map[structure]

def get_mean_value_per_structure(group, structure_ids, hemisphere):
    means = []
    structs = []
    if hemisphere == 'L':
        df = Ldat
    elif hemisphere == 'R':
        df = Rdat
    isids = df[(df['SampleID'] == group)]['SampleID'].values
    print(isids)
    for structure_id in structure_ids:
        if len(df[(df['structure_id'] == structure_id) & 
                       (df['SampleID'].isin(isids))]) > 0:
            str_mean = df[(df['structure_id'] == structure_id) & 
                           (df['SampleID'].isin(isids))]['density'].mean()
            means.append(str_mean+1)
            structs.append(int(structure_id))
    
    # Not all structures are represented in this Ldataset. 
    # Color missing structures by their parent
    structures = pd.DataFrame(st.get_structures_by_set_id([184527634]))
    missing_strs = [structure for structure in structures['id'].values if 
                    structure not in structs]
    print(len(missing_strs))
    for structure in missing_strs:
        ancestors = structures[structures['id'] == structure]['structure_id_path'].values[0]
        for n in range(2, len(ancestors)+1):
            ancestor = ancestors[-n:][0]
            if ancestor in df['structure_id'].unique():
                str_mean = df[(df['structure_id'] == ancestor) & 
                           (df['SampleID'].isin(isids))]['density'].mean()
                means.append(str_mean)
                structs.append(structure)
                missing_strs.remove(structure)
    for structure in missing_strs: # fill in the rest with zeros
        means.append(0)
        structs.append(structure)

    structuredat = dict(zip(structs, means))
    # Color cortical layer structures by their parent
    layer_structs = st.get_structures_by_set_id([667481440, 667481441, 667481445, 667481446, 667481449])
    layer_ids = [structure['id'] for structure in layer_structs]
    layer_ids += [1121, 526, 20, 52, 543, 664, 92, 712, 139, 727, 28, 60, 743] #ENTl, ENTm
    for structure in layer_ids:
        parent = st.parent([structure])[0]['id']
        str_mean = df[(df['structure_id'] == parent) & 
                           (df['SampleID'].isin(isids))]['density'].mean()
        structuredat[structure] = str_mean
        structuredat[997] = 0
        structuredat[129] = 0
        structuredat[140] = 0
        structuredat[81] = 0
        structuredat[153] = 0
        structuredat[145] = 0
    return structuredat, len(isids)

def get_cmap(group, hemisphere, scale=1):
    structure_vals, n = get_mean_value_per_structure(group, Ldat['structure_id'].unique(), hemisphere)
    vals = np.array([value for key, value in structure_vals.items()])
    maxi = np.quantile(vals, 0.95)
    scaled_vals = [value/maxi for key, value in structure_vals.items()]
    rgb_vals = dict(zip(structure_vals.keys(), scaled_vals))
    for key in structure_vals:
        rgb_vals[key] = tuple([255*i for i in cm.BuPu(structure_vals[key]*scale)[:3]])
        rgb_vals[0] = (255, 255, 255)
    return rgb_vals, n

def main(params):
    reference_space =  mcc.get_reference_space()
    structure_vals_L, n = get_mean_value_per_structure(params['SampleID'], 
                                                     Ldat.structure_id.unique(), 'L')
    structure_vals_R, n = get_mean_value_per_structure(params['SampleID'], 
                                                     Rdat.structure_id.unique(), 'R')
    rgb_vals_L, n = get_cmap(params['SampleID'], 'L')
    rgb_vals_R, n = get_cmap(params['SampleID'], 'R')
    image = [0,0,0]
    image[0] = reference_space.get_slice_image(0, params['xcoord'], rgb_vals_L) #posterior
    print(image[0].shape)
    image[1] = np.flip(np.rot90(reference_space.get_slice_image(2, params['zcoord'], rgb_vals_L)), 0) #right
    image[2] = np.rot90(reference_space.get_slice_image(1, params['ycoord'], rgb_vals_L))   #inferior
    fig = plt.figure(figsize=(8, 3), facecolor='w')
    columns = 3
    rows = 1
    for i in range(columns*rows):
        fig.add_subplot(rows, columns, i+1)
        f = plt.imshow(image[i], cmap = cm.BuPu)
        plt.axis('off')
    cbar = fig.colorbar(f, fraction=0.046)
    cbar.set_label('pTau Probability (per mm$\mathregular{^{3}}$)', 
                   rotation=90, fontsize = 10)
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), fontsize = 10)
    print(structure_vals_L)
    maxval = max(structure_vals_L.items(), key=operator.itemgetter(1))[1]
    print(maxval)
    cbar_vals = [0, 
                             np.round(maxval*.2, -3), 
                             np.round(maxval*.4, -3),
                             np.round(maxval*.6, -3),
                             np.round(maxval*.8, -3),
                             np.round(maxval, -3)]
    cbar_vals = [int(value/1e3) for value in cbar_vals]
    cbar.ax.set_yticklabels(cbar_vals)
    if not os.path.exists(params['savepath']):
        os.mkdir(params['savepath'])
    mouse_line_fname = str(params['SampleID'])
    plt.savefig(os.path.join(params['savepath'], 
                             '{0}_map_{1}-{2}-{3}.png'.format(mouse_line_fname,
                             params['xcoord'],
                             params['zcoord'],
                             params['ycoord'])), 
            facecolor=fig.get_facecolor(),
            bbox_inches='tight', 
            format='png', 
            dpi=300)
    
if __name__ == '__main__':

    with open(os.path.join('params', 'image_params_single_brain.json'), 'r') as data_file:    
        parameters = json.load(data_file)
    parameters['SampleID'] = 297

    main(parameters)