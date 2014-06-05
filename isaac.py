"""
 -----------------------------------------------
 Some simple python code to be easily imported from python
 -----------------------------------------------
 """
import pynbody
SimArray = pynbody.array.SimArray

import copy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import warnings
import glob
import os
import datetime

self_dir = os.path.dirname(os.path.realpath(__file__))
print os.path.realpath(__file__)

def walltime(filename):
    """
    Reads walltime information from a ChaNGa .log file.
    
    ** ARGUMENTS **
    
    filename : str
        Filename of the .log file to load
        
    ** RETURNS **
    
    wall_per_step : array
        Wall time per step in seconds
    """
    
    log_file = np.genfromtxt(filename, comments='#', delimiter=' ')
    wall_per_step = log_file[:,-1]
    walltime_total = datetime.timedelta(seconds = wall_per_step.sum())
    walltime_avg = datetime.timedelta(seconds = wall_per_step.mean())

    print 'Total walltime: '
    print str(walltime_total)
    print 'Average walltime per step:'
    print str(walltime_avg)
    
    return wall_per_step
    

def load_acc(filename, param_name = None):
    """
    Loads accelerations from a ChaNGa acceleration file (.acc2), ignoring the
    star particle.
    
    IF param_name is None, a .param file is searched for, otherwise param_name
    should be a string specifying a .param file name
    
    IF no param_file is found, the defaults are used:
        length unit: AU
        mass unit  : Msol
    """
    
    if param_name is None:
        
        prefix = filename.split('.')[0]
        
        param_list = glob.glob('*' + prefix +'*param')
        
        if len(param_list) > 0:
            
            param_name = param_list[0]
            
        elif len(glob.glob('*.param')) > 0:
            
            param_name = glob.glob('*.param')[0]
            
        else:
            
            warnings.warn('Could not find .param file.  Assuming default units')
            
    if param_name is not None:
        
        # If a param name is set or a param file has been found:
        print 'Loading param file: {}'.format(param_name)
        param = configparser(param_name, ftype='param')
        
    else:
        
        # Set the default parameters
        param = {}
        # Assume AU as length unit
        param['dKpcUnit'] = pynbody.units.au.ratio('kpc')
        # Assume mass units as Msol
        param['dMsolUnit'] = 1.0
    
    # Load acceleration file as numpy array
    acc_in = np.genfromtxt(filename)
    n_particles = int(acc_in[0])
    acc_in = acc_in[1:]
    
    if (len(acc_in) % n_particles) == 0:
        
        n_dim = len(acc_in)/n_particles
        
    else:
        
        raise IOError, 'Number of entries is not int multiple of num. particles'
        
    # Reshape the input accelerations, ignoring the star particle
    acc_2D = acc_in.reshape([n_particles, n_dim], order='F')
    acc_2D = acc_2D[0:-1]
    
    # Figure out units
    G = pynbody.units.G
    l_unit = param['dKpcUnit']*pynbody.units.kpc
    m_unit = param['dMsolUnit']*pynbody.units.Msol
    t_unit = ((l_unit**3) * G**-1 * m_unit**-1)**(1,2)
    a_unit = l_unit * t_unit**-2
    # Cast acc_2D as a SimArray with units a_unit
    acc_sim = match_units(acc_2D, a_unit)[0]
    
    return acc_sim
    
        
def sigma(snapshot, bins=100):
    """
    Calculates surface density vs r (relative to the star)
    """
    
    # Begin by subtracting off the star position
    star_pos = snapshot.star['pos'].copy()
    pos = snapshot.gas['pos'].copy()
    n_particles = pos.shape[0]
    # Make star_pos a matrix and subtract it off
    pos -= np.dot(np.ones([n_particles,1]), star_pos)
    r = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
    # particle mass
    m_gas = snapshot.gas['mass'][[0]]
    
    N, r_bins = np.histogram(r, bins=bins)
    r_bins = match_units(r_bins, r)[0]
    r_center = (r_bins[1:] + r_bins[0:-1])/2
    dr = r_bins[[1]] - r_bins[[0]]
    
    sig = N*m_gas/(2*np.pi*r_center*dr)
    
    return sig, r_center
    

def Q(snapshot, molecular_mass = 2.0, bins=100):
    """
    Calculates the Toomre Q, binned as a function of radius for a snapshot.
    snapshot should be either a filename or a pynbody simsnap.
    
    molecular_mass is the mean molecular mass divided by amu
    
    bins is either the number of radial bins to use or the radial bin edges
    
    RETURNS
    
    Returns a tuple (Q, r_bin_centers)
    where Q is evaluated at r_bin_centers
    """
    # If snapshot is a filename, load it.  Otherwise, assume it is a snapshot
    if isinstance(snapshot, str):
        
        snapshot = pynbody.load(snapshot)
        
    # Set up constants
    G = SimArray(1.0,'G')
    kB = SimArray(1.0,'k')
    # Begin by subtracting off the star position
    star_pos = snapshot.star['pos'].copy()
    pos = snapshot.gas['pos'].copy()
    n_particles = pos.shape[0]
    # Make star_pos a matrix and subtract it off
    pos -= np.dot(np.ones([n_particles,1]), star_pos)
    r = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
    # Load other quantities
    m_star = snapshot.star['mass']
    m_gas = snapshot.gas['mass'][[0]]
    T = snapshot.gas['temp']
    # Calculate sound speed
    m_mol = molecular_mass * SimArray(1.0,'m_p') # mean molecular mass in amu
    c_s = np.sqrt(kB*T/m_mol)
    # Find mass interior to every particle
    n_int = np.array(((r.argsort()).argsort()).tolist(), dtype='int')
    m_int = m_star + n_int*m_gas
    # Calculate sigma (surface density)
    n_per_bin, r_binedges = np.histogram(r, bins=bins)
    r_binedges = SimArray(r_binedges, r.units)
    dr = r_binedges[[1]] - r_binedges[[0]]
    r_bincenters = (r_binedges[0:-1] + r_binedges[1:])/2.0
    sigma = m_gas*n_per_bin/(2*np.pi*r_bincenters*dr)
    sigma_units = sigma.units
    sigma = interp.interp1d(r_bincenters, sigma, kind='linear', bounds_error=False)
    sigma = SimArray(sigma(r), sigma_units)
    # Calculate epicyclic frequency (kappa)
    k = np.sqrt(G * (m_int + 2*np.pi*sigma*r**2)/r**(3,1))
    # calculate Q and return as numpy array
    Q = c_s*k/(np.pi * G * sigma)
    Q.convert_units('1')
    Q = np.array(Q.tolist())
    # Bin Q and take the average
    Q_binned = binned_mean(r, Q, binedges=r_binedges)[1]
    
    return Q_binned, r_bincenters
    
def strip_units(x):
    """
    Removes the units from a SimArray and returns as a numpy array.  Note
    that x is copied so that it is not destroyed
    
    x can be a single SimArray or a tuple or list of SimArrays
    
    If any of the inputs are single number, they are returned as a number
    
    USAGE:
    
    array = strip_units(SimArray)
    
    array1, array2, ... = strip_units([SimArray1, SimArray2, ...])
    """
    if isinstance(x, (tuple,list)):
        
        # loop through and assign output
        x_out = []
        
        for x_i in x:
            
            if np.prod(x_i.shape) == 1:
                # There is only one element in x_i.  Make sure to return it as
                # a number  (not an array)
                if np.sum(x_i.shape) == 0:
                    # This is a zero dimensional SimArray
                    x_out.append(x_i.tolist())
                else:
                    # This is 1 dimensional SimArray
                    x_out.append(x_i[0])

            else:
                
                #This is a multi-element SimArray
                x_out.append(np.asarray(x_i.tolist()))
        
    else:
        
        if np.prod(x.shape) == 1:
            # There is only one element in x_i.  Return as a number
            if np.sum(x.shape) == 0:
                # This is a 0 dimensional SimArray
                x_out = x.tolist()
            else:
                # This a 1 dimensional SimArray
                x_out = x[0]
                
        else:
            
            x_out = np.asarray(x.tolist())
        
    return x_out
    
def set_units(x, units):
    """
    Sets the units of x to units.  If x has units, they are ignored.
    Does not destroy/alter x
    
    USAGE:
    
    SimArray = set_units(x, units)
    
    SimArray1, SimArray2, ... = set_units([x1, x2, ...], units)
    
    SimArray1, SimArray2, ... = set_units([x1, x2, ...], [units1, units2, ...])
    """
    if isinstance(x, (tuple,list)):
        
        x_out = []
        
        if not isinstance(units, (tuple, list)):
            
            units = [units]*len(x)
        
        for i in range(len(x)):
            
            x_i = x[i]
            
            if pynbody.units.has_units(x_i):
                
                x_i_array = strip_units(x_i)
                x_out.append(SimArray(x_i_array, units[i]))
                
            else:
                
                x_out.append(SimArray(x_i, units[i]))
            
    else:
        
        if pynbody.units.has_units(x):
            
            x_array = strip_units(x)
            x_out = SimArray(x_array, units)
            
        else:
            
            x_out = SimArray(x, units)
        
    return x_out
            
def make_param(snapshot, filename=None):
    """
    Generates a default param dictionary.  Can be saved using isaac.configsave
    
    EXAMPLE
    
    snapshot = pynbody.load('snapshot.std')  # Load snapshot
    param_dict = isaac.make_param(snapshot)  # Make default param dict
    isaac.configsave(param_dict, 'snapshot.param', ftype='param') # Save
    
    Optionally, the user can set the snapshot filename manually
    """
    fname_def = os.path.join(self_dir, 'default.param')
    param = configparser(fname_def, ftype='param')
    
    if filename is not None:
        
        param['achInFile'] = filename
        param['achOutName'] = os.path.splitext(filename)[0]
        
    elif snapshot.filename != '<created>':
        
        param['achInFile'] = snapshot.filename
        param['achOutName'] = os.path.splitext(snapshot.filename)[0]
        
    # Set up the length units
    param['dKpcUnit'] = snapshot['pos'].units.ratio('kpc')
    # Set up the mass units
    param['dMsolUnit'] = snapshot['mass'].units.ratio('Msol')
    # Set the mean molecular mass
    param['dMeanMolWeight'] = snapshot.gas['mu'][0]
    
    return param
    
def match_units(x, y):
    """
    Matches the units of x to y and returns x and y in the same units.
    
    IF x and y don't have units, they are unchanged
    
    IF one of x or y has units, the unitless quantity is returned as a 
    SimArray (see pynbody.array.SimArray) with the units of the other quantity.
    
    IF both have units, then an attempt is made to convert x into the units of
    y.  If this is not possible, an error is raised, for example if x is in
    units of 'au' and y is in units of 'Msol'
    
    x, y can be: scalar, array, SimArray, pynbody unit (eg pynbody.units.G),
        or a unit string (eg 'Msol a**-2')
    
    
    *** RETURNS ***
    
    x, y are returned as a tuple
    """
    # ----------------------------------------------
    # Check if either is a string
    # ----------------------------------------------
    if isinstance(x, str):
        
        x = SimArray(1.0, x)
        
    if isinstance(y,str):
        
        y = SimArray(1.0, y)
    
    # ----------------------------------------------
    # Check if one is a pynbody unit
    # ----------------------------------------------
    # If one is a named unit (eg pynbody.units.G), convert to SimArray
    if isinstance(x, pynbody.units.UnitBase):
        
        x = SimArray(1.0, x)
        
    if isinstance(y, pynbody.units.UnitBase):
        
        y = SimArray(1.0, y)
        
    # ----------------------------------------------
    # Check the units
    # ----------------------------------------------
    # If both have units, try to convert x to the units of y
    if (pynbody.units.has_units(x)) & (pynbody.units.has_units(y)):
        
        x_out = (x.in_units(y.units))
        y_out = y
    
    # If only x has units, make y a SimArray with the units of x
    elif (pynbody.units.has_units(x)):
        
        y_out = SimArray(y, x.units)
        x_out = x
        
    # If only y has units, make x a SimArray with the units of y
    elif (pynbody.units.has_units(y)):
        
        x_out = SimArray(x, y.units)
        y_out = y
    
    # Otherwise, neither has units
    else:
        
        x_out = x
        y_out = y
        
    # Try to copy so that changing x_out, y_out will not change x,y
    try: 
        
        x_out = x_out.copy()
        
    except AttributeError:
        
        pass
    
    try: 
        
        y_out = y_out.copy()
        
    except AttributeError:
        
        pass
    
    return x_out, y_out
    
def digitize_threshold(x, min_per_bin = 0, bins=10):
    
    """
    Digitizes x according to bins, similar to numpy.digitize, but requires
    that there are at least min_per_bin entries in each bin.  Bins that do not
    have enough entries are combined with adjacent bins until they meet the
    requirement.
    
    **ARGUMENTS**
    
    x : array_like
        Input array to be binned.  Must be 1-dimensional
    min_per_bin : int
        Minimum number of entries per bin.  Default = 0
    bins : int or sequence of scalars, optional
        [same as for np.histogram]
        If bins is an int, it defines the number of equal-width bins in the 
        given range (10, by default). If bins is a sequence, it defines the 
        bin edges, including the rightmost edge, allowing for non-uniform bin 
        widths.
        
    **RETURNS**
    
    A tuple containing:
    ind : array_like
        Indices of the bin each element of x falls into, such that:
        bin_edges[i] <= x[i] < bin_edges[i+1]
        (See np.digitize, this uses the same convention)
    bin_edges: array_like
        The edges of the bins
    """
    
    # Find number in each bin
    N, bin_edges = np.histogram(x, bins)
    
    if N.sum() < min_per_bin:
        
        raise RuntimeError,'Not enough particles within the bin range'
        
    n_bins = len(bin_edges) - 1
    
    # Find out which binedges to delete
    edge_mask = np.ones(len(bin_edges), dtype='bool')    
    
    for i in range(n_bins - 1):
        # Work forwards
        
        if N[i] < min_per_bin:
            
            # Set mask to not use the right bin edge
            edge_mask[i+1] = False
            # Combine the particles in current and next bin
            N[i] += N[i+1]
            N[i+1] = N[i]
            
    bin_mask = edge_mask[1:]
    N = N[bin_mask]
    bin_edges = bin_edges[edge_mask]
    edge_mask = np.ones(len(bin_edges), dtype='bool')
    n_bins = len(bin_edges) - 1
    
    for i in range(n_bins-1, 0, -1):
        # Work backwards
        
        if N[i] < min_per_bin:
            
            # Set mask to not use the left bin edge
            edge_mask[i] = False
            # Combine the particles in current and next bin
            N[i] += N[i-1]
            N[i-1] = N[i]
            
    bin_edges = bin_edges[edge_mask]
    ind = np.digitize(x, bin_edges)
    
    return ind, bin_edges

def binned_mean(x, y, nbins=10, binedges = None, weights=None,\
weighted_bins=False):
    """
    Bins y according to x and takes the average for each bin.  If binedges is
    specified, the x-bins are defined by binedges.  Otherwise the x-bins are
    determined by nbins
    
    If weights = None, equal weights are assumed for the average, otherwise
    weights for each data point should be specified
    
    y_err (error in y) is calculated as the standard deviation in y for each
    bin, divided by sqrt(N), where N is the number of counts in each bin
    
    IF weighted_bins is True, the bin centers are calculated as a center of
    mass
    
    NaNs are ignored for the input.  Empty bins are returned with nans
    
    RETURNS a tuple of (bin_centers, y_mean, y_err)
    """
        
    if binedges is not None:
        
        nbins = len(binedges) - 1
        
    else:
        
        binedges = np.linspace(x.min(), (1 + np.spacing(2))*x.max(), nbins + 1)
        
    if weights is None:
        
        weights = np.ones(x.shape)

    weights = strip_units(weights)
    
    # Pre-factor for weighted STD:
    A = 1/(1 - (weights**2).sum())
    
    
    # Initialize
    y_mean = np.zeros(nbins)
    y_std = np.zeros(nbins)
    # Find the index bins for each data point
    ind = np.digitize(x, binedges) - 1
    # Ignore nans
    nan_ind = np.isnan(y)
    N = np.histogram(x, binedges)[0]
    
    # Initialize bin_centers (try to retain units)
    bin_centers = 0.0*binedges[1:]
    
    for i in range(nbins):
        
        #Indices to use
        mask = (ind==i) & (~nan_ind)
        # Set up the weighting
        w = weights[mask].copy()
        w /= w.sum()
        A = 1/(1 - (w**2).sum())
        #y_mean[i] = np.nanmean(y[mask])
        y_mean[i] = (w * y[mask]).sum()
        var = A*(w*(y[mask] - y_mean[i])**2).sum()
        y_std[i] = np.sqrt(var)
        #y_std[i] = np.std(y[use_ind])
        
        if weighted_bins:
            # Center of mass of x positions
            bin_centers[i] = (w*x[mask]).sum()
        
    y_mean = match_units(y_mean, y)[0]
    y_err = y_std/np.sqrt(N)
    y_err = match_units(y_err, y)[0]

    y_mean[N==0] = np.nan
    y_err[N==0] = np.nan
    
    if not weighted_bins:
        
        bin_centers = (binedges[0:-1] + binedges[1:])/2.0
        bin_centers = match_units(bin_centers, x)[0]
        
    else:
        
        bin_centers[N==0] = np.nan
    
    return bin_centers, y_mean, y_err
    
def heatmap(x, y, z, bins=10, plot=True, output=False):
    """
    Creates a pcolor heatmap for z evaluated at (x,y).  z is binned and
    averaged according to x and y.  x, y, and z should be 1-D arrays with the
    same length.
    
    IF bins = N, a pcolor plot of shape (N,N) is returned
    IF bins = (M,N) [a tuple], a pcolor plot of shape (M,N) is returned
    
    IF plot = True (default) a plot is created.
    
    *** RETURNS ***
    IF output = False, nothing is returned (default)
    
    IF output = True:
    
    Returns x_mesh, y_mesh, z_binned
    
    x_mesh, y_mesh are the meshgrid x,y edges z is evaluted in.  z_binned is
    the average of z for each bin.
    """
    
    N, x_binedges, y_binedges = np.histogram2d(x, y, bins = bins)
    x_ind = np.digitize(x, x_binedges) - 1
    y_ind = np.digitize(y, y_binedges) - 1
    nx_bins = len(x_binedges) - 1
    ny_bins = len(y_binedges) - 1
    z_binned = np.zeros([nx_bins, ny_bins])
    
    for i in range(nx_bins):
        
        for j in range(ny_bins):
            
            z_binned[i,j] = z[(x_ind==i) & (y_ind==j)].mean()
            
    x_mesh, y_mesh = np.meshgrid(x_binedges, y_binedges, indexing = 'ij')
    
    if plot:
        
        cmap = copy.copy(matplotlib.cm.jet)
        cmap.set_bad('w',1.)
        masked_z = np.ma.array(z_binned, mask=np.isnan(z_binned))
        plt.pcolor(x_mesh, y_mesh, masked_z, cmap = cmap)
        plt.colorbar()

    if output:
        
        return x_mesh, y_mesh, z_binned     
    

def configparser(fname,ftype='auto'):
    """
     --------------------------------------------------
        parameters = configparser(fname,ftype='auto')
        
    Tries to parse ChaNGa configuration files
    
    ftype can be 'auto', 'param', or 'director'.  If auto, config parser will
    try to determine the filetype on its own.
    
    returns:
        dictionary 'parameters'.  The keys are the names of the parameters and 
        the values are the values defined in the file fname
     --------------------------------------------------
     """
    types = np.array(['param','director'])
    ftype = ftype.lower()
    param = {}
    if ftype == 'auto':
        # Try using extension to determine file type
        a = fname.split('.')
        ftype = a[-1].lower()
    if np.sum(types == ftype) == 0:
        # Could not find file type
        print ('Could not determine config filetype...exiting')
        return param
        # Try to determine filetype
    # Echo filetype
    print 'config filetype: {0}'.format(ftype)
    # --------------------------------------------------
    # Parse param file
    # --------------------------------------------------
    if ftype == 'param':
        farray = np.genfromtxt(fname,delimiter='=',dtype=None)
        for n in range(len(farray)):
            param[farray[n,0].strip()] = str2num(farray[n,1].strip())
    # --------------------------------------------------
    # Parse director file
    # --------------------------------------------------
    elif ftype == 'director':
        f = open(fname,'r')
        f.seek(0)
        dummy = 0
        for line in f:
            a = line.strip().split()
            if len(a) == 1:
                # we're dealing with a flag
                param[a[0]] = str2num(a[0])
            elif len(a) > 1:
                param[a[0]] = str2num(a[1:])
            else:
                # This is an empty line
                pass
        f.close()
    # --------------------------------------------------
    # Throw warning, return 'param' as empty
    # --------------------------------------------------
    else:
        warnings.warn('Still cannot determine filetype.')
    return param

def configsave(param,filename,ftype='auto'):
    """
     --------------------------------------------------    
    Saves parameters defined by param (see configparser) to filename.
    Possible ftypes are 'director' and 'param'.  If set to auto, configsave
    tries to guess file type from the extension.
     --------------------------------------------------
     """
    f = open(filename,'w')
    types = np.array(['param','director'])
    ftype = ftype.lower()
    if ftype == 'auto':
        a = filename.split('.')
        ftype = a[-1].lower()
    if ftype == 'param':
        pars = sorted(param.iteritems())
        for n in range(len(pars)):
            f.write('{0:25s}= {1}\n'.format(pars[n][0],pars[n][1]))
    elif ftype == 'director':
        values = param.values()
        keys = param.keys()
        for n in range(len(keys)):
            outstr = keys[n]
            if outstr == values[n]:
                # We just have a flag
                pass
            elif isinstance(values[n],(float,int,str)):
                outstr = outstr + ' {0}'.format(values[n])
            else:
                outstr = outstr + ' ' + ' '.join(map(str,values[n]))
            f.write('{0}\n'.format(outstr))
    else:
        #no file type
        warnings.warn('no such filetype {0}\nCould not save'.format(ftype))
    f.close()
    

# DEPRECATED (SLOW)
#def extrap1d(x,y):
#    """
#    Calculates a linear interpolation of x and y and does a linear
#    extrapolation for points outside of x and y.
#    Uses scipy.interpolate.interp1d
#    """
#    # Ignore nans
#    ind = (~np.isnan(x)) & (~np.isnan(y))
#    x = x[ind]
#    y = y[ind]
#    # calculate interpolation
#    yspline = interp.interp1d(x,y,kind='linear')
#    
#    def pointwise(x0):
#        if x0 < x.min():
#            return y[0] +  (x0 - x[0])*(y[1]-y[0])/(x[1]-x[0])
#        elif x0 > x.max():
#            return y[-1] + (x0 - x[-1])*(y[-1] - y[-2])/(x[-1] - x[-2])
#        else:
#            return yspline(x0)
#    
#    def ufunclike(x):
#        return np.array(map(pointwise,np.array(x)))
#    
#    return ufunclike
    
def extrap1d(x,y):
    """
    Calculates a linear interpolation of x and y and does a linear
    extrapolation for points outside of x and y.
    Uses scipy.interpolate.interp1d
    """
    # Ignore nans
    ind = (~np.isnan(x)) & (~np.isnan(y))
    x = x[ind]
    y = y[ind]
    # calculate interpolation
    yspline = interp.interp1d(x,y,kind='linear')
    
    def fcn(x0):
        
        mask1 = x0 < x.min()
        mask2 = x0 > x.max()
        out = np.zeros(len(x0))
        out[mask1] = y[0] +  (x0[mask1] - x[0])*(y[1]-y[0])/(x[1]-x[0])
        out[mask2] = y[-1] + (x0[mask2] - x[-1])*(y[-1] - y[-2])/(x[-1] - x[-2])
        mask3 = (~mask1) & (~mask2)
        out[mask3] = yspline(x0[mask3])
        
        return out        
    
    return fcn
            
def smoothstep(x,degree=5,rescale=False):
    """
    Calculates a smooth step function y(x) evaluated at the data points x.
    x should be a numpy array or float.  
    
    y(x) is a polynomial of order 'degree' (default is 5).  degree must be an
    odd number between 3 and 25 (inclusive).  The higher the order, the 
    sharper the step is.
    
    y(x) is defined by:
        y(0) = 0
        y(1) = 1
        The first (degree - 1)/2 derivatives are 0 at y = 0,1
        
    *** ARGUMENTS ***
    
    * x * Points at which to evaluate the smoothstep
    
    * degree * Degree of the smooth step.  Must be odd number between 3 and 25
        default = 5
        
    * rescale *  Rescale x to be between 0 and 1.  Default = False.  If True,
        x MUST be an array (greater than length 1)
    
    
    *** RETURNS ***
    
    """
    # -----------------------------------------------------------
    # Load up the hermite spline (polynomial) coefficients
    # -----------------------------------------------------------
    fname = os.path.join(self_dir,'hermite_spline_coeffs.dat')
    f =open(fname,'r')
    
    coeffs_list = []
    order_list = []
    
    for line in f:
        
        l = line.strip().split(',')
        order_list.append(int(l[0]))
        
        for n in range(len(l)):
            
            l[n] = float(l[n].strip())
            
        coeffs_list.append(np.array(l[1:],dtype='float'))
    
    order = np.array(order_list)
    coeffs = coeffs_list[(order==degree).argmax()]
    # -----------------------------------------------------------
    # Calculate the smooth step function y(x)
    # -----------------------------------------------------------
    n_coeffs = len(coeffs)
    
    if rescale:
        
        try:
            x = (x - x.min())/(x.max() - x.min())
        except:
            raise RuntimeError,'Could not rescale x.  Make sure x is an array'
    
    if isinstance(x, (int, long, float, complex)):
        
        # x is a number, handle accordingly
        y = 0.0
        
        if (x > 0) & (x < 1):
            # If 0<x<1, calculate the smooth step
            for n in range(n_coeffs):
                
                y += coeffs[n] * x**(degree - n)
                
        elif x <= 0:
            
            y = 0.0
            
        else:
            
            y = 1.0
        
    else:
        
        # Assume x is a numpy array
        y = np.zeros(x.shape)
        ind = (x > 0) & (x < 1)
        
        for n in range(n_coeffs):
            
            y[ind] += coeffs[n] * x[ind]**(degree-n)
            
        y[x >= 1] = 1
    
    return y
    

def str2num(string):
    """
     --------------------------------------------------
     Tries to see if 'string' is a number
     
     If 'string' is a string, returns:
       int(string) for integers
       float(string) for floats
       'string' otherwise
    
     If 'string' is a float or an integer, returns:
       string
    
     If none of the above, treats it like a list or tuple
     and returns for each entry of 'string' a float,int,
     or str as required.  Returns as a list
     --------------------------------------------------
     """
    if isinstance(string,int):
        output = string
    elif isinstance(string,float):
        output = string
    elif not isinstance(string,str):
        output = []
        for a in string:
            try:
                output.append(int(a))
            except:
                try:
                    output.append(float(a))
                except:
                    output.append(a)
        if len(output) == 1:
            output = output[0]
    else:
        output = string
        try:
            output = int(string)
        except:
            try:
                output = float(string)
            except:
                pass
    return output

def loadhalos(fname=''):
    """
     Load halo (.grp) file generated from fof
     Should be an ascii list of numbers, where the first row contains the
     total number of particles (gas+star+dark) and the remaining rows define
     which halo each particle belongs to
     """
    if fname == '':
        # Empty filename
        pass
    grp = np.loadtxt(fname,dtype=int)
    grp = grp[1:]   # (ignore the number of particles)
    
    return grp

def fof(fFilter,saveDir='',minMembers=8,linklen=0.01):
    """
     --------------------------------------------------
     A simple script that allows you to loop through calls to fof
     for many files in one directory
     --------------------------------------------------
     """
    flist = np.sort(glob.glob(fFilter))
    nfiles = len(flist)
    if (saveDir != '') and nfiles > 0:
        if ~os.path.isdir(saveDir):
            os.makedirs(saveDir)
    
    for n in range(nfiles):
        fname = flist[n]
        outname = os.path.join(saveDir,fname)
        os.system('totipnat < {0} | fof -g -m {1} -e {2} -o {3}'.format(fname,minMembers,linklen,outname))