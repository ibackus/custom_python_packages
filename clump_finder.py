# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 12:10:08 2014

@author: ibackus
"""

import numpy as np
import pynbody
import isaac

import subprocess

import glob
import os

def find_clumps(f, n_smooth=32, param=None, arg_string=None, logfile_name='clump_finder.log'):
    """
    Uses skid (https://github.com/N-BodyShop/skid) to find clumps in a gaseous
    protoplanetary disk.  
    
    The linking length used is equal to the gravitational softening length of
    the gas particles.
    
    The density cut-off comes from the criterion that there are n_smooth particles
    within the Hill sphere of a particle.  This is formulated mathematically as:
    
        rho_min = 3*n_smooth*Mstar/R^3
        
    where R is the distance from the star.  The trick used here is to multiply
    all particle masses by R^3 before running skid so the density cut-off is:
    
        rho_min = 3*n_smooth*Mstar
        
    **ARGUMENTS**
    
    *f* : TipsySnap
    A tipsy snapshot loaded/created by pynbody
    
    *n_smooth* : int (optional)
    Number of particles used in SPH calculations.  Should be the same as used
    in the simulation.  Default = 32
    
    *param* : dict (optional)
    param dictionary (see isaac.configparser)
    
    *arg_string* : str (optional)
    Additional arguments to be passed to skid.  Cannot use -tau, -d, -m, -s, -o
    
    *logfile_name* : str (optional)
    Filename to log the stdout of skid
    
    **RETURNS**
    
    *clumps* : array, int-like
    Array containing the group number each particle belongs to, with star
    particles coming after gas particles.  A zero means the particle belongs
    to no groups
    """
        
    # Estimate the linking length as the gravitational softening length
    tau = f.g['eps'][0]
    
    # Calculate minimum density
    rho_min = 3*n_smooth*f.s['mass'][0]
    
    # Center on star.  This is done because R used in hill-sphere calculations
    # is relative to the star
    star_pos = f.s['pos'].copy()
    f['pos'] -= star_pos
    
    # Scale mass by R^3
    R = isaac.strip_units(f['rxy'])
    m0 = f['mass'].copy()
    f['mass'] *= (R+tau)**3
    
    # Save temporary snapshot
    f_prefix = str(np.random.randint(np.iinfo(int).max))
    f_name = f_prefix + '.std'
    
    if param is not None:
        
        param_name = f_prefix + '.param'
        isaac.configsave(param, param_name)
        
    f.write(filename=f_name, fmt=pynbody.tipsy.TipsySnap)
        
    f['pos'] += star_pos
    f['mass'] = m0
    
    # Open up a log file
    logfile = open(logfile_name, 'w')
    
    # Run skid
    command = 'totipnat < {} | skid -tau {:.2e} -d {:.2e} -m {:d} -s {:d} -o {}'\
    .format(f_name, tau, rho_min, n_smooth, n_smooth, f_prefix)
    print '\n', command
    p = subprocess.Popen(command, shell=True, stdout=logfile, stderr=logfile)
    p.wait()
    
    # Load clumps
    clumps = isaac.loadhalos(f_prefix + '.grp')
    
    # Cleanup
    logfile.close()
    for name in glob.glob(f_prefix + '*'):
        
        os.remove(name)
        
    return clumps