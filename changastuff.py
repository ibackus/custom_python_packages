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
    """
    Submits scriptname, contained in all the directories, to the submission
    queue using qsub.
    
    **ARGUMENTS**
    
    directories : str, list, ndarray, tuple
        The directory or directories containing the submission script
    scriptname : str
        Script to submit.  Should be present in all the directories
    """
    
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

def make_continue_sub(simdirs='.', paramname='snapshot.param', \
newparam='continue.param', simtime=1, nSteps=None, oldsub='subber.sh', \
newsub='cont_subber.sh'):
    """
    Makes a submission script for continuing a simulation from a previous output.
    Also makes a .param file for the continued simulation.  The simulation
    will be continued in the same directory, with snapshot numbering scheme
    for outputs being the same.
    
    Parameters for the original simulation cannot be altered (except the number
    of steps you want to continue the simulation by)
    
    Any checkpoints will be deleted.
    
    Requires a submission script be present for the original simulation
    
    
    **ARGUMENTS**
    
    simdirs : str or iterable containing strings
        The directory/directories containing the simulation(s) to continue
    paramname : str
        Filename of the .param file for the simulation(s).  Must be the same 
        for every simulation
    newparam : str
        filename for the .param file for the continued simulation
    simtime : float
        How many times the original nSteps to continue by.  IE, the new total
        number of steps (including the original simulation) is calculated as:
            nSteps_new = nSteps_old * (simtime + 1)
        if simtime=0, then the simulation will just continue and finish and
        the original number of steps.  (similar to checkpointing)
    nSteps : int
        OVERRIDES simtime.  Sets the total number of extra steps to simulate.    
    oldsub : str
        Filename for the original submission script (should be the same for
        every simdir)
    newsub : str
        Filename for the new submission script (should be the same for every 
        simdir)
    
    """
    
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
        
        # Delete old checkpointssubmitted 8 hours ago by Aikukun to /r/gifs

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