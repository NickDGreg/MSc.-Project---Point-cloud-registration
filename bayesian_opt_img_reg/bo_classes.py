# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:09:04 2017

@author: Nick
"""
import numpy as np 
import matlab
from pointcloud_manipulation_functions import transform_points
from math import inf
class BoRegistration(object):
    """Stores variables to track bayesian optimisation process of trajectories and
       contains methods to perform the registration and error calculations
    
    Used for bayesian optimisation of sequential trajectories from t1,..,tN where tN 
    is full trajectory and t1,..,tN-1 are sequential subsets of the trajectory.
    The variables track the registration for each subset and allow previous trajectories
    and transformations to be used to alter how the current subset is registered.
    
    Attributes:
        aorta: the 3D point cloud of the aorta we are registering to 
        eng: matlab engine to call the point cloud registration functions
        perform_reg: matlab function e.g. eng.ICP_finite() that performs image registration
        constraints: an area with in which the centering position needs to be
        iteration: tracks number of bayesian optimisation iterations
        traj_list: stores all the trajectories seen to calculate centering 
                        position error with current trajectory and all past subsets
        traj: the most recent trajectory added
        traj_center: the center of the trajectory, used to calculate translation 
                     to centering position
        sparse: bool showing whether the trajectory is sparse or full
        sparse_idxs: bool array the same length as traj idicating which points are 
                     part of the sparse array
        sparse_error: used when the trajectory being registered is the full trajectory,
                      sparse version of the registered full trajectory is used to calculate
                      sparse_error. This way the error of the full trajectory registration
                      can be compared to the sparse trajectory registration. Otherwise even 
                      with same registration, full trajectory returns larger error as it 
                      has more points in center of aorta and therefore further from the walls
        GP: trained Gaussian Process used to predict sample locations  
        all_past_checks: boolean deciding if error measure is calculated for only 
                         the transformation and current trajectory or all past 
                         trajectories as well
        limited_past_checks: boolean deciding if error measure is for the transformation,
                             current trajectory and one previous trajectory
        reg_prev_best_transform: boolean deciding if when initialising the Bayesian
                                 optimisation the previous best transformation is 
                                 used as one if the initialisation points
        bo_opt_variables: dictionary storing the variables for the most recent Bayesian
                          optimisation round. 
                          opt_from_prev_wreg describes if the best found
                          came from the previous registration or a new one. 
                          opt_from_prev_noReg describes if best found came from 
                          the previous registration but without any registration applied
        global_opt_variables: same as above but for the entire BO registration process
                              not just the most recent
        
        Output:
            The important attributes are bo and global _opt_variables. These describe
            the quality of the registration, make comparisons to the ground truth and 
            allow different registration methods to be compared.
    """
    def __init__(self, aorta, eng, constraints=None):
        self.aorta = matlab.double(aorta.tolist())
        self.eng = eng
        self.perform_reg = None
        self.constraints = constraints 
        self.iteration = 0
        self.traj_list = []
        self.traj = None
        self.traj_center = None
        self.sparse = False
        self.sparse_idxs = None
        self.sparse_error = None
        self.GP = None    
        self.all_past_checks = False
        self.limited_past_checks = False
        self.reg_prev_best_transform = False
        self.bo_opt_variables = {'opt_error':inf, 'opt_center':np.zeros((1,3)), 
                                 'opt_transform':np.eye(4), 'opt_from_prev_wreg':False}
        self.global_opt_variables = {'opt_error':inf, 'opt_center':np.zeros((1,3)), 
                                     'opt_transform':np.eye(4), 'opt_from_prev_wreg':False, 
                                     'opt_from_prev_noReg':False}
        
    def add_traj(self, trajectory):
        self.traj_list.append(trajectory)
        self.traj = trajectory
        self.traj_center = np.mean(trajectory, axis=0)
        
    def set_prev_best_transform_check(self):
        self.reg_prev_best_transform = True
    
    def set_limited_past_checks(self):
        self.limited_past_checks = True
        self.all_past_checks = False
    
    def set_all_past_checks(self):
        self.limited_past_checks = False
        self.all_past_checks = True
        
    def set_GP(self, GP):
        self.GP = GP
        
    def set_global_opt_vars(self, opt_error, opt_center, opt_transform, opt_from_prev):
        self.global_opt_variables = {'opt_error':opt_error,'opt_center':opt_center, 
                                     'opt_transform':opt_transform, 'opt_from_prev_wreg':opt_from_prev}
        
    def set_bo_opt_vars(self, opt_error, opt_center, opt_transform, opt_from_prev):
        self.bo_opt_variables = {'opt_error':opt_error, 'opt_center':opt_center, 
                                 'opt_transform':opt_transform, 'opt_from_prev_wreg':opt_from_prev}
    
    def get_simple_reg_error_and_points(self, center_position):
        """register the current trajectory to the center position, return
        error and registered points. Only possible for one center position.
        """
        assert center_position.shape == (1,3)
        self.iteration += 1
        traj_matlab, difference = self.get_traj_matlabarray(center_position, self.traj, self.traj_center)
        errors = np.zeros((center_position.shape[0],1))
        idx = 0
        for trajectory in traj_matlab: #list comprehension wouldnt work on mlarray
            registered_points, _, regError = self.perform_reg(self.aorta, trajectory, nargout=3)
            errors[idx] = regError / trajectory.size[0]
            idx += 1
        return errors, registered_points
        
    def get_reg_error(self, center_positions):
        """registers trajectory with center_positions as init location. Contains 
        options to alter how the registration is performed. See class attributes.
        """
        print(self.iteration)
        self.iteration += 1
        num_trajs = len(self.traj_list)
        errors = np.zeros((1,center_positions.shape[0]))
    
        if self.reg_prev_best_transform:
            #saved transformation already moves points to CP but need to recalc CP for new larger trajectory 
            #as CP input is CP for the previous trajectory
            prev_opt_traj = transform_points(self.global_opt_variables['opt_transform'], self.traj)
            trajs_matlab, differences = self.get_traj_matlabarray(center_positions, prev_opt_traj, np.mean(prev_opt_traj, axis=0)) 
        else:
            trajs_matlab, differences = self.get_traj_matlabarray(center_positions, self.traj, self.traj_center) 
            
        diff = np.ones((1,4))
        idx = 0
        for tm in trajs_matlab: 
            errors_added = 0
            registered_points, transformation, reg_error = self.perform_reg(self.aorta, tm, nargout=3)
            errors[0,idx] += reg_error
            errors_added += 1
            transformation = np.array(transformation)
            
            if self.reg_prev_best_transform: 
                transformation = np.dot(transformation, self.global_opt_variables['opt_transform'])                          
            diff[0,:3] = differences[idx]
            translate_t = np.eye(4)
            translate_t[:,3] = diff
            transformation = np.dot(transformation, translate_t)
                
            if self.all_past_checks:
                for traj_idx in range(num_trajs - 1):
                    errors[0, idx] += self.get_past_traj_error(transformation, traj_idx)
                    errors_added += 1
            
            if self.limited_past_checks and num_trajs >= 3:
                if num_trajs % 2 == 0:
                    traj_idx = int(num_trajs / 2)
                else:
                    traj_idx = int((num_trajs + 1) / 2)
                errors[0, idx] += self.get_past_traj_error(transformation, traj_idx)
                errors_added += 1  
            
            errors[0,idx] = errors[0,idx]/ errors_added        
            if errors[0,idx] < self.bo_opt_variables['opt_error']:
                print('NEW OPTIMAL')
                if self.reg_prev_best_transform:
                    self.set_bo_opt_vars(errors[0,idx], center_positions[idx], transformation, True) 
                else:
                    self.set_bo_opt_vars(errors[0,idx], center_positions[idx], transformation, False) 
                if not self.sparse: 
                    self.sparse_error = self.get_sparse_error(transformation, self.traj)
                    
            if self.reg_prev_best_transform: self.reg_prev_best_transform = False #only use previous T once
            idx += 1
        return errors
    
    def align_trajs_with_CP(self, center_positions):
        print(self.iteration)
        self.iteration += 1
        trajs_matlab, differences = self.get_traj_matlabarary(center_positions, self.traj, self.traj_center)
        errors = np.zeros((center_positions.shape[0],1))
        idx = 0
        for tm in trajs_matlab: #list comprehension wouldnt work on mlarray
            regError = self.eng.closestPoint_error(self.aorta, tm, nargout=1)
            errors[idx] += regError / tm.size[0]
            if self.error_all_past_trajs:
                for traj_idx in range(len(self.traj_list) - 1):
                    past_traj = np.asarray([self.traj_list[idx] + differences[idx]])
                    past_traj_matlab = matlab.double(past_traj.tolist())
                    errors[idx] += self.eng.closestPoint_error(self.aorta, past_traj_matlab, nargout=1) / past_traj.shape[0]
                errors[idx] = errors[idx]/ len(self.traj_list) 
                idx += 1
        return(errors)
        
    def get_past_traj_error(self, transformation, traj_idx):
        transformed_traj = transform_points(transformation, self.traj_list[traj_idx])
        transformed_traj_matlab = matlab.double(transformed_traj.tolist())
        return self.eng.closestPoint_error(self.aorta, transformed_traj_matlab, nargout=1) / transformed_traj.shape[0] 
    
    def get_sparse_error(self, transformation, full_traj):
        sparse_traj = [points for points,select in zip(full_traj, self.sparse_idxs) if select]
        transformed_traj = transform_points(transformation, np.array(sparse_traj))
        transformed_traj_matlab = matlab.double(transformed_traj.tolist())
        return self.eng.closestPoint_error(self.aorta, transformed_traj_matlab, nargout=1) / transformed_traj.shape[0] 
    
    def update_global_opts(self):
        """checks if bayesian optimisation registration is better than simply taking the previous 
        best transformation found and applying that without any new registration. Updates 
        global optimal to be whichever of current transformation or previous is best
        """
        current_traj_globalT = transform_points(self.global_opt_variables['opt_transform'], self.traj)
        current_traj_globalT_list = current_traj_globalT.tolist()
        current_traj_globalT_matlab = matlab.double(current_traj_globalT_list)
        current_error = self.eng.closestPoint_error(self.aorta, current_traj_globalT_matlab, nargout=1) / self.traj.shape[0]
        #unnormalise previous error to calc new error using best found transformation so far
        updated_global_error = (((self.global_opt_variables['opt_error'] * (len(self.traj_list) - 1) ) + current_error))/ len(self.traj_list)
        if updated_global_error < self.bo_opt_variables['opt_error']:
            print('previous BO better than new BO')
            self.set_global_opt_vars(updated_global_error, 
                                     self.global_opt_variables['opt_center'],
                                     self.global_opt_variables['opt_transform'],
                                     False)
            self.global_opt_variables['opt_from_prev_noReg'] = True
        else:
            print('new BO better than previous')
            self.set_global_opt_vars(self.bo_opt_variables['opt_error'], 
                                     self.bo_opt_variables['opt_center'], 
                                     self.bo_opt_variables['opt_transform'],
                                     self.bo_opt_variables['opt_from_prev_wreg'])
            self.global_opt_variables['opt_from_prev_noReg'] = False
        
    def get_traj_matlabarray(self, center_positions, traj, traj_center):
        differences = center_positions - traj_center
        trajs = np.asarray([traj + diff for diff in differences])
        trajs_matlab = matlab.double(trajs.tolist())
        return trajs_matlab, differences