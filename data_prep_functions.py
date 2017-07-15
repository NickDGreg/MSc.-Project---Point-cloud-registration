# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 16:26:55 2017

@author: User
"""
import numpy as np

def get_datasets_and_info(all_data):
    stat_names = list(filter(None, all_data[0]))
    stat_idxs = list(filter(None, all_data[1]))
    stat_idxs = [int(x) for x in stat_idxs]
    stat_names_idxs = dict(zip(stat_names,stat_idxs))
    stat_idxs = stat_idxs[:-1] #remove the final index which is not incl in X_train
    data = np.array(all_data[4:])
    data = data.astype(float)

    procedure_idxs = np.zeros((1,13))[0];
    current = 1;
    idx_found = 0
    for i in range(data.shape[0]):
        if data[i,0] != current:
            idx_found = idx_found + 1
            procedure_idxs[idx_found] = i
            current = data[i,0]
    procedure_idxs[-1] = data.shape[0]

    X_data = data[:,4:-3]
    y_data = data[:,-3:]
    return X_data, y_data, stat_idxs, stat_names_idxs, procedure_idxs
    
    
def get_trajectories(traj_file):
    traj_data = np.genfromtxt(traj_file, delimiter=',')
    traj_indices = traj_data[1,:]
    trajs = traj_data[2:,:]
    return traj_indices, trajs