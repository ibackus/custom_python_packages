# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 13:30:56 2014

@author: ibackus
"""

import os
import glob
import isaac
import re
import subprocess

def subber(directories, scriptname):
    
    # Make the directories list iterable
    if isinstance(directories, str):
        
        directories = [directories]
    
    # Get current working directory
    cwd = os.getcwd()
    
    # Change directories to full paths
    fullpaths = []    
    for directory in directories:
        
        fullpaths.append(os.path.abspath(directory))
     
    # Submission command
    command = 'qsub ' + scriptname
    
    # Submit all the scripts
    for fullpath in fullpaths:
        
        os.chdir(fullpath)
        
        p = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Print output
        for line in iter(p.stdout.readline, ''):
            
            print line,
            
        p.wait()
        
    # Return to current working directory
    os.chdir(cwd)

def make_continue_sub(simdirs='.', paramname='snapshot.param', simtime=1, \
oldsub='subber.sh', nSteps=None, newparam='continue.param', newsub='cont_subber.sh'):
    
    # Lazy man's way of dealing with files in another directory
    cwd = os.getcwd()
    
    # Make simdirs a list if needed
    if isinstance(simdirs, str):
        
        simdirs = [simdirs]
        
    # Iterate over all the simulation directories present
    for simdir in simdirs:
        
        os.chdir(simdir)
        
        # Load param file
        param = isaac.configparser(paramname, 'param')
        fprefix = param['achOutName']
        
        # Find all the outputs.  They should be of the format fprefix.000000
        search_exp = '^' + fprefix + '.(?:(?<!\d)\d{6}(?!\d))$'
        flist = []
        
        for fname in glob.glob(fprefix + '*'):
            
            if re.search(search_exp, fname) is not None:
                
                flist.append(fname)
        
        # Find the number of the last output (the last 6 chars should be an int)
        flist.sort()
        iStartStep = int(flist[-1][-6:])
        param['iStartStep'] = iStartStep
        param['achInFile'] = flist[-1]    
        
        # Set the number of steps to run
        if nSteps is None:
            
            nSteps = int(param['nSteps'] * simtime)
            
        param['nSteps'] = iStartStep + nSteps
        
        # Save new param file
        isaac.configsave(param, newparam, ftype='param')
        
        # Delete old checkpoints
        for checkpoint in glob.glob(fprefix + '.chk*'):
            
            print 'removing ' + checkpoint
            os.system('rm -rf ' + checkpoint)
            
        if os.path.exists('lastcheckpoint'):
            
            print 'removing lastcheckpoint'
            os.remove('lastcheckpoint')
        
        # Create a submission script for the simulation continuation
        oldsubfile = open(oldsub, 'r')
        newsubfile = open(newsub, 'w')
        
        for line in oldsubfile:
            
            newsubfile.write(line.replace(paramname, newparam))
            
        oldsubfile.close()
        newsubfile.close()
        
        # Make the submission script executable
        os.chmod(newsub, 0777)
        
        # Change back to original working directory
        os.chdir(cwd)
    
    return param