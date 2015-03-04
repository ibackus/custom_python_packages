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
import ICgen_utils
import ICglobal_settings
global_settings = ICglobal_settings.global_settings

def PBS_script(workdir, param='snapshot.param', nodes=1, ppn=12, walltime=48, \
jobname='PBS_job', backfill=True, email=None):
    """
    A not very robust function to generate PBS submission scripts for ChaNGa
    jobs on hyak.  Some of the requirements include:
    
    mpi version of ChaNGa
    gcc_4.4.7-ompi_1.6.5
    PBS (it's on hyak, don't worry)
    
    By default, any directory with 'lastcheckpoint' in it will be restarted!
    If you don't want to just restart a simulation, delete lastcheckpoint from
    the simulation directory!
    
    **ARGUMENTS**
    
    *Required*
    
    workdir : str
        Directory of the simulation
        
    *Optional*
    
    param : str
        Filename of the .param file (not full path)
    nodes : int
        Number of computation nodes to request
    ppn : int
        Number of cores per node
    walltime : float or int
        Walltime to request in hours
    jobname : str
    backfill : bool
        Boolean flag for whether to use the backfill (default is TRUE)
    email : str
        Email address for PBS to send updates to
        
    **RETURN**
    
    PBS_script : str
        A string for the PBS script.  Can be saved easily to file
    """
    
    # Setup filenames
    param_full = os.path.join(workdir, param)
    outfile = os.path.join(workdir, 'changa.out')
    fprefix = os.path.splitext(param)[0]
    
    # Get changa preset
    preset = global_settings['changa_presets']['mpi']
    
    # Set up the walltime for PBS
    hours = int(walltime)
    mins = int((walltime*60)%60)
    secs = int((walltime*3600)%60)
    walltime_str = '{0:d}:{1:02d}:{2:02d}'.format(hours,mins,secs)
    
    # Write the script!
    
    # Start with the shebang
    script = '#!/bin/bash\n'
    
    # Some PBS flags
    script += '#PBS -N {0}\n\
#PBS -j oe\n\
#PBS -l nodes={1}:ppn={2},feature={2}core\n\
#PBS -l walltime={3}\n\
#PBS -V\n'.format(jobname, nodes, ppn, walltime_str)
    
    # Email stuff
    if email is not None:
        
        script += '#PBS -M {0}\n'.format(email)
        script += '#PBS -m be\n'
        
    # Choose whether to use the backfill
    if backfill:
        
        script += '#PBS -q bf\n'
        
    # Runtime initialization
    script += 'module load gcc_4.4.7-ompi_1.6.5\n\
export MX_RCACHE=0\n\
workdir={0}\n\
cd $workdir\n\
changbin=$(which {1})\n'.format(workdir, preset[2])
    
    # Now assume that we want to restart if there is a checkpoint
    script += 'if [ -e "lastcheckpoint" ]\n\
then\n\
    echo "lastcheckpoint exists -- restarting simulation..."\n\
    last=`cat lastcheckpoint`\n\
    mpirun --mca mtl mx --mca pml cm $changbin +restart {0}.chk$last +balancer MultistepLB_notopo -wall {1} {2} >> {3} 2>&1\n'.format(fprefix,int(walltime*60), param_full, outfile)
    script += 'else\n\
    echo "lastcheckpoint doesnt exist -- starting new simulation..."\n\
    mpirun --mca mtl mx --mca pml cm $changbin -D 3 +consph +balancer MultistepLB_notopo -wall {0} {1} >& {2}\n\
fi\n'.format(int(walltime*60), param_full, outfile)
    
    return script

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