# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 14:15:20 2014

@author: ibackus
"""

import numpy as np
import scipy.optimize as opt

def fit(x, y, p0 = None, xtol=0.0001, ftol=0.0001,maxiter=None):
    
    if p0 is None:
        
        # Default initial guesses
        p0 = [1.0, 1.0, 0.0]
        
    def res(p):
        
        return ((y - p[0]*x**p[1] - p[2])**2).sum()
        
    p_f = opt.fmin(res, p0, xtol=xtol, ftol=ftol, maxiter=maxiter)
    
    return p_f