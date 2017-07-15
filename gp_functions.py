# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 17:13:16 2017

@author: User
"""
import numpy as np

def posterior_samples_self_def(model,X,numSamples,full_cov=True):
    numSamplingPoints = X.shape[0]
    m,v = model._raw_predict(X, full_cov=full_cov) #can put a different kernel in here for the predictions? see cell above^
#     m,v = _raw_predict_TESTS(model,X, full_cov=full_cov)
    if model.normalizer is not None:
        m, v = model.normalizer.inverse_mean(m), model.normalizer.inverse_variance(v)
    
    def sim_one_dim(m, v):
        if not full_cov:
            return np.random.multivariate_normal(m.flatten(), np.diag(v.flatten()), numSamples).T
        else:
            return np.random.multivariate_normal(m.flatten(), v, numSamples).T
            
    fsim = np.empty((model.output_dim, numSamplingPoints, numSamples))
    for d in range(model.output_dim):
        if full_cov and v.ndim == 3:
            fsim[d] = sim_one_dim(m[:, d], v[:, :, d])
        elif (not full_cov) and v.ndim == 2:
            fsim[d] = sim_one_dim(m[:, d], v[:, d])
        else:
            fsim[d] =np.random.multivariate_normal(m[:, d].flatten(), v, numSamples).T

    return(fsim)