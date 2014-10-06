# -*- coding: utf-8 -*-
"""
Contains a routine to automatically estimate a reasonable time-step size for
ChaNGa.  The idea is to have about half the particles fall into the lowest
rung.  This is done by calculating the rung distribution for a large time step
by running ChaNGa

NOTE: this alpha version requires ChaNGa to kill itself after printing out
the rung distribution to the standard out

USAGE:

import find_dDelta

dDelta = find_dDelta.calc_dDelta('snapshot.param', dDelta0=10.0)

where dDelta0 is some large time step, big enough that all particles are above
the lowest rung

Created on Thu Sep 18 13:00:12 2014

@author: ibackus
"""

import numpy as np
#import isaac
import ICgen_utils
import re
#import os
#import signal

def calc_dDelta(param_name, dDelta0 = 10.0):
    """
    NOTE: Still in alpha and requires my version of ChaNGa which kills itself
    after printing the rung distribution
    
    Estimates a good time step (dDelta) for ChaNGa simulations
    """
    
    N = calc_rungdist(param_name, dDelta0)
    rungs = np.arange(len(N))
    med = histmedian(rungs, N)
    dDelta = dDelta0 * 2.0**(-med+1)
    
    return dDelta

def read_rungdist(p):
    """
    Parses the stdout froma ChaNGa subprocess to find the first instance of
    'rung dist'.  From this line it reads the rung distribution
    """
    
    rung_line = ''
    
    for line in iter(p.stdout.readline, ''):
        
        line = line.strip()
        print line
        
        if 'rung dist' in line.lower():
            
            rung_line = line
            #os.killpg(p.pid, signal.SIGTERM)
            break
        
    if rung_line == '':
        
        raise RuntimeError, 'Could not find rung distribution'
        
    rung_list = re.findall('\d+', rung_line)
    rung_array = np.array(rung_list).astype(int)
    
    return rung_array
    
def calc_rungdist(param_name, dDelta):
    """
    Calls ChaNGa to calculate the rung distribution
    """
        
    changa_args = '+n 1 -dt {}'.format(dDelta)
    command = ICgen_utils.changa_command(param_name, changa_args=changa_args, preset='isaac')
    p = ICgen_utils.changa_run(command, verbose=False)
    
    return read_rungdist(p)
        
def histmedian(binedges, N):
    """
    Calculate a median from a histogram
    
    if binedges is the same length as N, it's assumed only the left binedges are
    present and an extra binedge is added at the end.
    """
    
    # Make sure everything is floats to prevent stupid errors    
    if not issubclass(binedges.dtype.type, float):
        
        binedges = binedges.astype(float)
        
    if not issubclass(N.dtype.type, float):
        
        N = N.astype(float)
        
        
    # Cumulative sum of the histogram
    s = np.cumsum(N)
    Ntot = s[-1]
    
    # Bin widths
    w = binedges[1:] - binedges[0:-1]
    
    if len(binedges) == len(N):
        
        dummy = np.zeros(len(binedges) + 1)
        dummy[0:-1] = binedges
        dummy[-1] = binedges[-1] + w[-1]
        binedges = dummy
        w = binedges[1:] - binedges[0:-1]
        
    # Find first bin which gives us more than half the total number
    for i, n in enumerate(s):
        
        if n >= 0.5*Ntot:
            
            ind = i
            break
        
    med = binedges[ind] + w[ind]*(0.5*Ntot - s[ind-1])/N[ind]
    
    return med