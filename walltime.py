#!/usr/bin/env python
"""
Created on Fri Mar 14 15:25:36 2014

@author: ibackus
"""

import time

t0 = time.time()
import matplotlib.pyplot as plt
import numpy as np
import datetime
import sys

t1 = time.time()
print 'Importing took {} s'.format(t1-t0)

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

    t2 = time.time()
    print 'Running took an extra {} s'.format(t2-t1)
    print 'For a total of {} s'.format(t2 - t0)
    
    plt.show()