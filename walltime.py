#!/usr/bin/env python
"""
Created on Fri Mar 14 15:25:36 2014

@author: ibackus
"""

import matplotlib.pyplot as plt
import numpy as np
import datetime
import sys

if len(sys.argv) < 2:
    
    print 'USAGE:     walltime filename'
    
else:
    
    fname = sys.argv[-1]
    log_file = np.genfromtxt(fname, comments='#', delimiter=' ')
    walltime_total = datetime.timedelta(seconds = log_file[:,-1].sum())
    walltime_avg = datetime.timedelta(seconds = log_file[:,-1].mean())

    print 'Total walltime: '
    print str(walltime_total)
    print 'Average walltime per step:'
    print str(walltime_avg)
    
    plt.plot(log_file[:,-1],'x')
    plt.show()