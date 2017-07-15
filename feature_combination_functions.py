# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 16:28:35 2017

@author: Nick
"""
import numpy as np
from time import time
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from itertools import combinations
from pointcloud_manipulation_functions import transform_points, register_pointclouds
import matplotlib.pyplot as plt
import GPyOpt
from math import ceil, inf

def bo_optimisation_errors(BO_reg, selected_error_function, traintest_data, 
                           all_trajs,sparse_idxs = None, sparse=False,
                           num_subsections=8, reg_prev_best_transform=False, 
                           plot_trajs=False, max_iter=5, aq_weight=40, 
                            rand_rotation = False, max_rotation=False):
    """Performs bayesian optimisation registration between the aorta and catheter points
    for n = num_subsections subsets of the catheter points to the full trajectory. Returns 
    the results of the registration.
    
    Attributes:
        BO_reg: Instance of the registration class containing the vars needed to perform 
                different versions of the registration
        selected_error_function: The method belonging to BO_reg chosen to perform the registration
        traintest_data: Dictionary containing X_train, X_test, y_train, y_test data
        all_trajs: Dictionary containing traj_indices (list indices in trajs & 
                    GT_trajs where each trajectory subset can be found), 
                    trajs (rotated trajectory and subsets),
                    GT_trajs (ground truth trajectory and subsets)
        sparse: bool showing whether the trajectory is sparse or full
        sparse_idxs: bool array the same length as traj idicating which points are 
                     part of the sparse array
        num_subsections: how many subsets of the trajectory to register and store the results of 
        max_iter: Number of iterations to run the BO for 
        aq_weight: Value to determine how much exploration is done by the BO
        rand_rotation: Bool indicating if trajectory has been randomly rotated 
        max_rotation: Bool indicating if trajectory was rotated by the maximum value (i.e. max of rand_rotation)
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
            bo_stats: Dictionary containing keys of the subset index 0-num_subsections and the registration 
                      results for that subset and a key of creation_params detailing the settings of the registration 
                      and the nature of the data used
    """
    num_trajs = int(all_trajs['trajs'].shape[1] / 3)
    steps = ceil(num_trajs/num_subsections)
    BO_reg.sparse = sparse
    bo_stats = {}
    traj_idx = -1
    for test_value in range(0, num_trajs, steps):
        traj_idx += 1
        traj, GT_traj = get_test_trajs(test_value, all_trajs)
        BO_reg.add_traj(traj)
        BO_reg.set_bo_opt_vars(inf, np.zeros((1,3)), np.eye(4), False)
        if reg_prev_best_transform: BO_reg.reg_prev_best_transform = True
        if not sparse:
            BO_reg.sparse_idxs = sparse_idxs[:traj.shape[0]]
        X_BO, y_BO = get_init_samples(BO_reg, selected_error_function, traintest_data)
        t1 = time()
        bounds = [{'name':'cathReg', 'type':'bandit', 'domain':traintest_data['y_train']}]
        Bopt = GPyOpt.methods.BayesianOptimization(f=selected_error_function,
                                                  domain=bounds,
                                                  X=X_BO,
                                                  Y=y_BO,
                                                  acquisition_type='LCB',
                                                  acquisition_weight=aq_weight)
        Bopt.run_optimization(max_iter)
        t2 = time()   
        if sparse:
            bo_opt_error = GPyOpt.util.general.best_value(Bopt.Y)[-1]
            bo_opt_full_error = 0
        else:
            bo_opt_error = BO_reg.sparse_error
            bo_opt_full_error = GPyOpt.util.general.best_value(Bopt.Y)[-1]
        print('Optimal X:' + str(Bopt.x_opt))
        BO_trans = transform_points(BO_reg.bo_opt_variables['opt_transform'], traj)
        diff_to_GT = register_pointclouds(GT_traj, BO_trans, BO_reg.eng, noreg_error=True) / BO_trans.shape[0]
        GT_error = register_pointclouds(BO_reg.aorta, GT_traj, BO_reg.eng, noreg_error=True) / GT_traj.shape[0]
        
        if plot_trajs:
            reg_error, BO_registered = BO_reg.get_simple_reg_error_and_points(np.array([Bopt.x_opt]))
            BO_registered = np.array(BO_registered)
            global_traj = transform_points(BO_reg.global_opt_variables['opt_transform'], traj)
            plot_registration(np.array(BO_reg.aorta), GT_traj, BO_registered, 
                              BO_trans, global_traj, point=Bopt.x_opt)
        
        BO_reg.update_global_opts()
        bo_stats[str(traj_idx)] = {'global_opt_error':BO_reg.global_opt_variables['opt_error'], 
                                    'bo_opt_error':bo_opt_error,
                                    'bo_opt_full_error': bo_opt_full_error,
                                    'diff_to_GT_error':diff_to_GT, 'GT_error':GT_error, 
                                    'comp_time':t2-t1, 'opt_from_prev_wreg':BO_reg.global_opt_variables['opt_from_prev_wreg'],
                                    'opt_from_prev_noReg':BO_reg.global_opt_variables['opt_from_prev_noReg'], 
                                    'opt_CP':np.ndarray.tolist(Bopt.x_opt)}

    bo_stats['creation_params'] = {'max_iters':max_iter,
                                   'reg_prev_best_transform':reg_prev_best_transform,
                                   'all_past_checks':BO_reg.all_past_checks,
                                   'limited_past_checks':BO_reg.limited_past_checks,
                                   'aquisition_weight':aq_weight,
                                   'sparse':sparse,
                                   'max_rotation':max_rotation,
                                   'rand_rotation':rand_rotation,
                                   'num_subsections':num_subsections}
    return bo_stats

def plot_registration(aorta, GT, BO_opt_traj=np.zeros((3,3)), 
                      BO_trans=np.zeros((3,3)), global_opt_traj=np.zeros((3,3)), 
                      extra=np.zeros((3,3)), point=np.array([0,0,0])):
                          
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(GT[:,0],GT[:,1],GT[:,2],
               marker="o",c='blue',s=10,label='Ground truth')    
    ax.scatter(BO_opt_traj[:,0],BO_opt_traj[:,1],BO_opt_traj[:,2],
               marker="o",c='g',s=20, label='BO registered') 
    ax.scatter(BO_trans[:,0],BO_trans[:,1],BO_trans[:,2],
               marker="o",c='r',s=20, label='BO transformed')           
    ax.scatter(global_opt_traj[:,0],global_opt_traj[:,1],global_opt_traj[:,2],
               marker="o",c='y',s=10, label='Previous global optimal') 
    ax.scatter(extra[:,0],extra[:,1],extra[:,2],
               marker="o",c='m',s=10, label='Extra') 
    ax.scatter(point[0],point[1],point[2],
               marker="o",c='r',s=30, label='Optimal centering location') 
    ax.scatter(aorta[:,0],aorta[:,1],aorta[:,2],
               marker="o",c='y',s=0.5,label='aorta')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1,
               ncol=2, mode="expand", borderaxespad=0.)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
      
def get_test_trajs(test_value, all_trajs):
    traj_col_idx = test_value * 3 
    traj_eidx = all_trajs['traj_indices'][traj_col_idx]
    traj = all_trajs['trajs'][:traj_eidx, traj_col_idx:traj_col_idx+3]
    GT_traj = all_trajs['GT_trajs'][:traj_eidx, traj_col_idx:traj_col_idx+3]
    return traj, GT_traj

def get_init_samples(BO_reg, selected_error, traintest_data=None, num_samples=5, GP=None):
    two_samples = False
    if GP is not None:
        print('to implement')
    else:
        if BO_reg.reg_prev_best_transform:
            prev_best = transform_points(BO_reg.global_opt_variables['opt_transform'], BO_reg.traj)
            prev_best_center = np.mean(np.array(prev_best), axis=0)
            y_BO_prev_best = selected_error(prev_best_center.reshape(1,3))
            two_samples = True
        X_BO = traintest_data['y_train'][np.random.choice(traintest_data['y_train'].shape[0], num_samples, replace=False), :]
        y_BO = selected_error(X_BO)
        if two_samples:
            X_BO = np.vstack((prev_best_center, X_BO))
            y_BO = np.concatenate((y_BO_prev_best, y_BO), axis=1)
        
    return X_BO, y_BO.reshape(-1,1)
        
def get_prediction(indices, stat_idxs, X_train, y_train, X_test):
    #alter to pass in a dict or object to reduce num of para,s
    #all data can definitely be grouped
    selected_data = np.array([])
    for i in range(len(indices)):
        selected_data = np.concatenate((selected_data, (np.arange(stat_idxs[indices[i]-1], stat_idxs[indices[i]]))))
    
    selected_data = [int(x) for x in selected_data]
    subset_Xtrain = X_train[:,selected_data]
        
    subset_Xtest = X_test[:,selected_data]

    lreg = LinearRegression()
    rreg = Ridge(alpha=0.5)
    lreg.fit(subset_Xtrain,y_train)
    rreg.fit(subset_Xtrain,y_train)
    lreg_predict = lreg.predict(subset_Xtest)
    rreg_predict = rreg.predict(subset_Xtest)
    
    return np.array([lreg_predict,rreg_predict])

def get_reg_error_fold_cross(predictions, y_test):
    errors = []
    for p in predictions:
        error = np.linalg.norm(np.mean(abs(p - y_test),axis=0))
        errors.append(error)    
    return np.asarray(errors)

def combi_errors_fold_cross(indices, stat_idxs,
                         X_data, y_data, cath_procedure_idx,
                         X_sim_data=None, y_sim_data=None, sim_procedure_idx=None,
                         traj_sim_train=False, sim_train=False):
    """Finds error for all possible feature combinations
    
    Performs K-fold cross validation in steps of 3 through the aortas. For each K,
    creates train and set sets of all possible feature combinations, trains 
    linear regression and ridge regression models and returns their errors.
    Final error is an average of the K test sets.
    Default assumption is that we train and test on only trajectory information.
    
    Args:
        indicies: array of 1,2,..,N with N being number of features
        trajSimTrain: True if want to train on trajectory info of other aortas
                      and simulated data of test aorta
        simTrain: True if only training on simulated data
    
    Returns:
        results: Dict containing feature set: error pairs. Error is an array of 
                 linear regression and ridge regression error
    """
    results = {}
    cross_val_loop = 0
    #perform aorta fold cross validation
    for i in range(1, len(cath_procedure_idx), 4):
        cross_val_loop += 1
        X_test = X_data[cath_procedure_idx[i-1]:cath_procedure_idx[i],:]
        y_test = y_data[cath_procedure_idx[i-1]:cath_procedure_idx[i],:]
        if sim_train:
            X_train = X_sim_data[sim_procedure_idx[i-1]:sim_procedure_idx[i],:]
            y_train = y_sim_data[sim_procedure_idx[i-1]:sim_procedure_idx[i],:]
        else:
            if i == 1:
                X_train = X_data[cath_procedure_idx[i]:,:]
                y_train = y_data[cath_procedure_idx[i]:,:]
            elif i == len(cath_procedure_idx):
                X_train = X_data[:cath_procedure_idx[i-1],:]
                y_train = y_data[:cath_procedure_idx[i-1],:]
            else:
                X_train = np.vstack((X_data[:cath_procedure_idx[i-1],:],
                                    X_data[cath_procedure_idx[i]:,:]))
                y_train = np.vstack((y_data[:cath_procedure_idx[i-1],:],
                                    y_data[cath_procedure_idx[i]:,:]))
            if traj_sim_train:
                X_train = np.vstack((X_train,
                                    X_sim_data[sim_procedure_idx[i-1]:sim_procedure_idx[i],:]))
                y_train = np.vstack((y_train,
                                    y_sim_data[sim_procedure_idx[i-1]:sim_procedure_idx[i]:,:]))

        for i in range(len(indices) + 1):
            for subset in combinations(indices, i):
                if len(subset) != 0:
                    predictions = get_prediction(subset, stat_idxs, X_train, y_train, X_test)
                    error = get_reg_error_fold_cross(predictions, y_test)
                    errorResult = results.setdefault(repr(subset), np.array([0.0, 0.0]))
                    errorResult += error
    
    results = {k: v / cross_val_loop for k, v in results.items()}
    return results
    
def get_data(aorta_idx, stat_idxs,
            X_data, y_data, cath_procedure_idx,
            X_sim_data=None, y_sim_data=None, sim_procedure_idx=None,
            traj_sim_train=False, sim_train=False):
    """Given number of aorta to use as test (1-13 excluding 7), 
    splits data into X,y train and test
    
    Not called from combiErrorsFoldCross as that must run through thousands of data rows and 
    the function call increases the time. This is used when visualising and testing different feature
    combinations and aortas as test.
    """
    X_test = X_data[cath_procedure_idx[aorta_idx-1]:cath_procedure_idx[aorta_idx],:]
    y_test = y_data[cath_procedure_idx[aorta_idx-1]:cath_procedure_idx[aorta_idx],:]
    if sim_train:
        X_train = X_sim_data[sim_procedure_idx[aorta_idx-1]:sim_procedure_idx[aorta_idx],:]
        y_train = y_sim_data[sim_procedure_idx[aorta_idx-1]:sim_procedure_idx[aorta_idx],:]
    else:
        if aorta_idx == 1:
            X_train = X_data[cath_procedure_idx[aorta_idx]:,:]
            y_train = y_data[cath_procedure_idx[aorta_idx]:,:]
        elif aorta_idx == len(cath_procedure_idx):
            X_train = X_data[:cath_procedure_idx[aorta_idx-1],:]
            y_train = y_data[:cath_procedure_idx[aorta_idx-1],:]
        else:
            X_train = np.vstack((X_data[:cath_procedure_idx[aorta_idx-1],:],
                                X_data[cath_procedure_idx[aorta_idx]:,:]))
            y_train = np.vstack((y_data[:cath_procedure_idx[aorta_idx-1],:],
                                y_data[cath_procedure_idx[aorta_idx]:,:]))
        if traj_sim_train:
            X_train = np.vstack((X_train,
                                X_sim_data[sim_procedure_idx[aorta_idx-1]:sim_procedure_idx[aorta_idx],:]))
            y_train = np.vstack((y_train,
                                y_sim_data[sim_procedure_idx[aorta_idx-1]:sim_procedure_idx[aorta_idx]:,:]))

    
    return X_train, y_train, X_test, y_test