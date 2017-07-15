# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 16:12:32 2017

@author: Nick
"""
import matlab
import numpy as np

def register_pointclouds(static_points, moving_points, eng, return_registered=False,noreg_error=False):
    """registers trajectory to aorta using method described by return_reg & CP_error
    """
    if isinstance(static_points, np.ndarray):
        static_points_list = static_points.tolist()
        static_points = matlab.double(static_points_list)
    if isinstance(moving_points, np.ndarray):
        moving_points_list = moving_points.tolist()
        moving_points = matlab.double(moving_points_list)
    if return_registered:
        registered_points, M, regError = eng.ICP_finite(static_points, moving_points, nargout=3)
        registered_points = np.array(registered_points._data).reshape(registered_points.size[::-1]).T
        return(regError,registered_points)
    elif noreg_error:
        return(eng.closestPoint_error(static_points, moving_points, nargout=1))
    else:
        _, _, regError = eng.ICP_finite(static_points, moving_points, nargout=3)
        return regError
        
def transform_points(M, points):
    transformed = np.zeros((points.shape))
    transformed[:,0] = points[:,0]*M[0,0] + points[:,1]*M[0,1] + points[:,2]*M[0,2] + M[0,3];
    transformed[:,1] = points[:,0]*M[1,0] + points[:,1]*M[1,1] + points[:,2]*M[1,2] + M[1,3];
    transformed[:,2] = points[:,0]*M[2,0] + points[:,1]*M[2,1] + points[:,2]*M[2,2] + M[2,3];
    return transformed
    
    
def shift_traj(centering_pos, traj):
    traj_center = np.mean(traj,axis=0)
    diff = centering_pos - traj_center
    return traj + diff 