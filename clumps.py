# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 12:10:08 2014

@author: ibackus
"""

# Multiprocessing modules
from multiprocessing import Pool, cpu_count
import multiprocessing as mp

# Generic pacakges
import numpy as np
import pynbody
SimArray = pynbody.array.SimArray
import subprocess
import glob
import os
import matplotlib.pyplot as plt

# 'Internal' packages
import isaac

def multilink(link_list):
    
    n_links = len(link_list)
    
    clump_list = []
    
    iStep = 0
    
    while iStep < n_links:
        
        pairs0 = link_list[iStep]
        new_mask = pairs0[:,0] == -1
        new_iord = pairs0[new_mask,1]
        
        for iord0 in new_iord:
            
            t = iStep + 1
            clump = [ [iord0, t] ]
            
            while t < n_links:
                
                pairs1 = link_list[t]
                iord = clump[-1][0]
                
                new_ind = np.nonzero(pairs1[:,0] == iord)[0]
                if len(new_ind) > 0:
                    # The clump links to something in the next timestep
                    new_ind = int(new_ind)
                    # Increment the time step
                    t += 1
                    clump.append([pairs1[new_ind, 1], t])
                    
                else:
                    # The clump links to nothing.  It has died
                    clump_list.append(np.array(clump))
                    break
                
        iStep += 1
                
    return clump_list

def link_clumps(clump_pars1, clump_pars2, link_thresh = 0.3):
    
    n1 = len(clump_pars1['iord'])
    n2 = len(clump_pars2['iord'])
    
    iord_list1 = list(clump_pars1['iord'])
    iord_list2 = list(clump_pars2['iord'])
    
    # Used to store how many particles clumps have in common
    connections = np.zeros([n1, n2], dtype=int)
    
    # Calculate the number of connections common to the clumps in clump_pars1
    # and the clumps in clump_pars2
    # Loop over the first set of clumps
    for i, iord1 in enumerate(iord_list1):
        # Loop over the second set of clumps
        iord1.sort()
        npart1 = len(iord1)
        for j, iord2 in enumerate(iord_list2):
            
            # Find which particles are common to clump[i] in pars1 and clump[j]
            # in pars2
            intersect = np.intersect1d(iord1, iord2, assume_unique=True)
            # Save how many are shared in the 2 clumps
            connections[i,j] = len(intersect)
            # Now only retain particles that are not common to the clumps
            # IE, we know where these particles end up, we can stop checking
            # them
            iord1 = np.setdiff1d(iord1, intersect, assume_unique=True)
            iord2 = np.setdiff1d(iord2, intersect, assume_unique=True)
            iord_list2[j] = iord2
            
            if len(iord1) < 1:
                # There are no more particles to look at in the original clump                
                break
            
        # Now ignore any connections where number of particles shared between
        # clumps less than link_thresh * (num. part. in clump 1)
        thresh_mask = connections[i,:] < link_thresh * npart1
        connections[i, thresh_mask] = 0
    
    # Find the clump in clump_pars2 that shares the most number of members for
    # a clump in clump_pars1.  This gives us the children of the clumps in
    # clump_pars1
    col_ind = connections.argmax(1)
    # Set all others to 0
    mask = np.zeros([n1,n2], dtype=bool)
    row_ind = np.arange(n1)
    mask[row_ind, col_ind] = True
    connections[~mask] = 0
    
    # The clumps in clump_pars2 may have multiple parents.  Select the parent
    # which shares the most particles in common (note, if this isn't the
    # child of any clump, parent_index will always be 0)
    parent_index = connections.argmax(0)
    # and find the number of particles inherited from the parent
    n_inherit = connections.max(0)
    
    # Demand to inherit at least 1 particle from the parent
    parent_index[n_inherit < 1] = -1
    
    # Now create the clump pairs
    clump_pairs = np.array( [parent_index, np.arange(n2)] ).T

    return clump_pairs
        

def clump_im(f, clump_array, width, qty='rho', resolution=1200, clim=None, clump_min=None):
    """
    Plots an sph image from f with particles not in clumps colored red and 
    particles in clumps colored green.  Uses pynbody for a backend.
    
    **ARGUMENTS**
    
    f : TipsySnapshot (see pynbody) or str
        The snapshot to plot OR the filename of a snapshot to plot
    clump_array : numpy array
        A array (same length as f) such that 0 entries correspond to particles
        not in clumps and entries > 0 correspond to particles in clumps
    width : str, simarray
        See pynbody.plot.sph.image.  Width of the plot, ie '3 au'
    resolution : int
        Resolution in pixels of the plot.  The plot will be res by res pixels
    clim : tuple,list,array
        Density limits arranged as [low, high].  Any pixels below low are mapped
        to 0, any pixels above high are mapped to 1.
    clump_min : float
        Used to set a cutoff for masking the clumps.  Not needed
        
    **RETURNS**
    
    image : numpy nd-array
        Returns an NxNx3 numpy array of the color image plotted.
    """
    # Check if the input is a filename
    if isinstance(f, str):
        
        f = pynbody.load(f)
    
    # Get the current state for matplotlib (this still needs work, since an
    # extra window will in general be created)
    current_fig = plt.gcf()
    interactive_flag = plt.isinteractive()
    plt.ioff()
    
    # Intermediate figure, required for rendering the plots
    fig1 = plt.figure()
    
    # Render a grayscale image of all the particles
    im_all = pynbody.plot.sph.image(f.g, width=width,resolution=resolution, cmap='gray', qty=qty)
    xlim = plt.xlim()
    ylim = plt.ylim()
    extent = [xlim[0], xlim[1], ylim[0], ylim[1]]
    
    fig1.clf()
    
    # Initialize the color image
    im_color = np.zeros([resolution, resolution, 3])
    # Set the red channel of the color image to be the plot of all particles
    im_color[:,:,0] = np.log(im_all)
    
    # Check to see that there is at least one clump
    clump_flag = (clump_array.max() > 0)
    
    if clump_flag:
        
        # Get a sub-snap of just particles in a clump
        mask = clump_array > 0
        f2 = f[mask]
        # Render an image of just particles in a clump
        im_clump = pynbody.plot.sph.image(f2.g, width=width,resolution=resolution, cmap='gray',qty=qty)
        # Set the clump image as the green channel
        im_color[:,:,1] = np.log(im_clump)
        
        plt.clf()
        
        # Most of the clump image should be blank space: igore everything
        # below a threshold
        if clump_min is None:
        
            clump_min = im_clump.mean()
        
        mask2 = im_clump > clump_min
    
    # Set the color scaling
    if clim is None:
        
        clim = [im_all.min(), im_all.max()]
        
    log_clim = [np.log(clim[0]), np.log(clim[1])]
            
    im_color -= log_clim[0]
    im_color /= (log_clim[1] - log_clim[0])
    im_color[im_color < 0] = 0
    im_color[im_color > 1] = 1
    
    if clump_flag:
        
        # Set all pixels outside a clump (in the clump image) to 0
        im_color[~mask2,1] = 0
        # Set all pixels inside a clump (in the overall image) to 0
        im_color[mask2,0] = 0
    
    else:
            
        im_color[:,:,1] = 0
        
    im_color[:,:,2] = 0
    
    # Plot
    plt.figure(current_fig.number)
    
    if interactive_flag:
        plt.ion()
    
    plt.imshow(im_color, extent=extent, interpolation='none', aspect='equal')
    
    # Echo the color limits used
    print 'clims used: {}'.format(clim)
        
    return im_color
    
def _parallel_clump_pars(args):
    
    return calc_clump_pars(*args)
    
def batch_clump_pars(flist, clump_list):
    
    nproc = cpu_count()
    
    arg_list = zip(flist, clump_list)
    
    pool = Pool(nproc)
    
    clump_pars = pool.map(_parallel_clump_pars, arg_list)
    pool.close()
    pool.join()
    
    return clump_pars

def calc_clump_pars(f, clump_nums):
    
    if clump_nums.max() < 1:
        # Return none if there are no clumps
        return
    
    if isinstance(f, str):        
        f = pynbody.load(f)
        
    try:
        
        iorder = f['iord']
        
    except KeyError:
        
        print 'Warning.  iorder not found.  Assuming 0,1,2,3...'
        iorder = np.arange(len(f))
        
    particle_nums = np.arange(len(f))
        
    # Only include particles in a clump AND that are not star particles
    mask1 = clump_nums > 0
    n_star = len(f.s)
    mask1[-(n_star+1):-1] = False
    clump_nums1 = clump_nums[mask1]
    f1 = f[mask1]
    iorder1 = iorder[mask1]
    particle_nums1 = particle_nums[mask1]
    
    # Get units set up
    m_unit = f1['mass'].units
    l_unit = f1['pos'].units
    v_unit = f1['vel'].units
    rho_unit = f1['rho'].units
    
    # Get arrays of pointers to the required quantities
    f_mass = f1['mass']
    f_pos = f1['pos']
    f_v = f1['vel']
    f_T = f1['temp']
    f_rho = f1['rho']
    
    # Initialize arrays
    n_clumps = clump_nums1.max()
    
    m = SimArray(np.zeros(n_clumps), m_unit) # clump mass
    N = np.zeros(n_clumps, dtype=int) # Number of particles/clump
    pos = SimArray(np.zeros([n_clumps,3]), l_unit) # center of mass
    r = SimArray(np.zeros(n_clumps), l_unit) # center of mass radial position
    v = SimArray(np.zeros([n_clumps, 3]), v_unit) # center of mass velocity
    # Angular momentum around the center of mass rest frame
    L = SimArray(np.zeros([n_clumps, 3]), m_unit*l_unit*v_unit)
    T = SimArray(np.zeros(n_clumps), 'K') # mass averaged temperature
    rho = SimArray(np.zeros(n_clumps), rho_unit) # density
    r_clump = SimArray(np.zeros(n_clumps), l_unit) # clump radius (size)
    
    # index of each particle (in this file)
    particle_ids = []
    # universal identity of each particle
    particle_iord = []
    
    # loop over the clumps
    for i in range(n_clumps):
        
        mask2 = (clump_nums1 == i+1)
        
        # Mask the input arrays to look at only the current clump
        p_mass = f_mass[mask2]
        p_pos = f_pos[mask2]
        p_v = f_v[mask2]
        p_T = f_T[mask2]
        p_rho = f_rho[mask2]
        
        # Calculate properties of the clump
        N[i] = mask2.sum()
        m[i] = p_mass.sum()
        pos[i] = np.dot(p_pos.T, p_mass[:,None]).flatten()
        pos[i] /= float(m[i])
        r[i] = np.sqrt((pos[i]**2).sum())
        v[i] = np.dot(p_v.T, p_mass[:,None]).flatten()
        v[i] /= float(m[i])
        
        # position of all particles relative to center of mass
        cm_pos = p_pos - pos[i]
        # velocity of all particles relative to center of mass
        cm_v = p_v - v[i]
        # angular momentum of all particles relative to center of mass
        cm_momentum = (cm_v * p_mass[:,None])
        p_L = np.cross(cm_pos, cm_momentum)
        # Total angular momentum relative to center of mass
        L[i] = p_L.sum(0)
        
        T[i] = p_T.sum()/N[i]
        rho[i] = p_rho.sum()/N[i]
        
        particle_ids.append(particle_nums1[mask2])
        particle_iord.append(iorder1[mask2])
        
    out_dict = {'m':m, 'N':N, 'pos':pos, 'r':r, 'v':v, 'L':L, 'T':T, 'rho':rho,\
    'r_clump': r_clump, 'ids': particle_ids, 'iord': particle_iord}
    
    return out_dict

def _parallel_find_clumps(args):
    """
    A wrapper to parallelize find_clumps()
    """    
    return find_clumps(*args)
    
def batch_clumps2(f_list, n_smooth=32, param=None, arg_string=None):
    """
    A parallel implementation of find_clumps.  Since SKID is not parallelized
    this can be used to run find_clumps on a set of snapshots from one
    simulation.
    
    **ARGUMENTS**
    
    f_list : list
        A list containing the filenames of snapshots OR the tipsy snapshots
        
    **RETURNS**
    
    clumps : list
        A list containing the clumps for each snapshot in f_list
    """
    # Number of processes to create = number of cores
    n_proc = cpu_count()
    
    # Set up the arguments for calls to find_clumps
    arg_list = []
    
    for i, f_name in enumerate(f_list):
        
        arg_list.append([f_name, n_smooth, param, arg_string, i])
        
    print arg_list
    
    # Set up the pool
    pool = Pool(n_proc)
    
    # Run the job in parallel
    results = pool.map(_parallel_find_clumps, arg_list, chunksize=1)
    pool.close()
    pool.join()
    
    return results

def find_clumps(f, n_smooth=32, param=None, arg_string=None, seed=None):
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
    
    *f* : TipsySnap, or str
        A tipsy snapshot loaded/created by pynbody -OR- a filename pointing to a
        snapshot.
    
    *n_smooth* : int (optional)
        Number of particles used in SPH calculations.  Should be the same as used
        in the simulation.  Default = 32
    
    *param* : dict (optional)
        param dictionary (see isaac.configparser)
    
    *arg_string* : str (optional)
        Additional arguments to be passed to skid.  Cannot use -tau, -d, -m, -s, -o
    
    *seed* : int
        An integer used to seed the random filename generation for temporary
        files.  Necessary for multiprocessing and should be unique for each
        thread.
    
    **RETURNS**
    
    *clumps* : array, int-like
        Array containing the group number each particle belongs to, with star
        particles coming after gas particles.  A zero means the particle belongs
        to no groups
    """
    # Parse areguments
    if isinstance(f, str):
        
        f = pynbody.load(f)
        
    if seed is not None:
        
        np.random.seed(seed)
        
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
    
    command = 'totipnat < {} | skid -tau {:.2e} -d {:.2e} -m {:d} -s {:d} -o {}'\
    .format(f_name, tau, rho_min, n_smooth, n_smooth, f_prefix)
#    print '\n', command
    #p = subprocess.Popen(command, shell=True, stdout=logfile, stderr=logfile)
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    for line in iter(p.stdout.readline, ''):
        print line,
    p.wait()
    
    # Load clumps
    clumps = isaac.loadhalos(f_prefix + '.grp')
    
    # Cleanup
    for name in glob.glob(f_prefix + '*'):
        
        os.remove(name)
        
    return clumps