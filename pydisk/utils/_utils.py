# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:19:03 2015

@author: ibackus
"""
import warnings
import logging
import os
import fnmatch
import numpy as np
import pynbody as pb
SimArray = pb.array.SimArray

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
    ftype = ftype.lower()
    if ftype == 'auto':
        # Try to figure out filetype
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
    
def units_from_param(param):
    """
    Figures out the simulation units from a .param file

    **ARGUMENTS**

    param : str or param dict (see configparser)
        Simulation .param file or param dict loaded by configparser
        Can also be a list or numpy array of these in which case a list
        of units dicts is returned

    **RETURNS**

    units : dict
        A dictionary of the units used in the simulation, returned as
        pynbody units
    """

    # Define function to load the units from a given param
    def _load_units(param):

        # Load param if necessary
        if isinstance(param, str):

            param = configparser(param, 'param')

        # Universal G
        G = pb.units.G

        # Load units
        dKpcUnit = param['dKpcUnit']
        dMsolUnit = param['dMsolUnit']

        # Set up pynbody units
        m_unit = pb.units.Unit('{0} Msol'.format(dMsolUnit))
        l_unit = pb.units.Unit('{0} kpc'.format(dKpcUnit))
        t_unit = (l_unit**3/(G*m_unit))**(1,2)

        # Convert the time unit to something sensible
        years = t_unit.in_units('yr')
        t_unit = pb.units.Unit('{0} yr'.format(years))

        # Return
        outdict = {'l_unit':l_unit, 'm_unit':m_unit, 't_unit':t_unit}
        return outdict

    # Iterate over param if necessary
    if isinstance(param, (list, np.ndarray)):

        outlist = []

        for par in param:

            outlist.append(_load_units(par))

        return outlist

    else:

        # Not iterable
        return _load_units(param)
        
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

            if pb.units.has_units(x_i):

                x_i_array = strip_units(x_i)
                x_out.append(SimArray(x_i_array, units[i]))

            else:

                x_out.append(SimArray(x_i, units[i]))

    else:

        if pb.units.has_units(x):

            x_array = strip_units(x)
            x_out = SimArray(x_array, units)

        else:

            x_out = SimArray(x, units)

    return x_out
    
def match_units(x, y):
    """
    Matches the units of x to y and returns x and y in the same units.

    IF x and y don't have units, they are unchanged

    IF one of x or y has units, the unitless quantity is returned as a
    SimArray (see pb.array.SimArray) with the units of the other quantity.

    IF both have units, then an attempt is made to convert x into the units of
    y.  If this is not possible, an error is raised, for example if x is in
    units of 'au' and y is in units of 'Msol'

    x, y can be: scalar, array, SimArray, pynbody unit (eg pb.units.G),
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
    # If one is a named unit (eg pb.units.G), convert to SimArray
    if isinstance(x, pb.units.UnitBase):

        x = SimArray(1.0, x)

    if isinstance(y, pb.units.UnitBase):

        y = SimArray(1.0, y)

    # ----------------------------------------------
    # Check the units
    # ----------------------------------------------
    # If both have units, try to convert x to the units of y
    if (pb.units.has_units(x)) & (pb.units.has_units(y)):

        x_out = (x.in_units(y.units))
        y_out = y

    # If only x has units, make y a SimArray with the units of x
    elif (pb.units.has_units(x)):

        y_out = SimArray(y, x.units)
        x_out = x

    # If only y has units, make x a SimArray with the units of y
    elif (pb.units.has_units(y)):

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
    
def findfiles(filefilter='*', basedir='.'):
    """
    Recursively find files according to filefilter

    ** ARGUMENTS **

    filefilter : str
        Filter for finding files.  ie, '*.jpg' or 'file.txt'

    basedir : str
        Base directory to search.  Default is the current directory

    ** RETURNS **

    files : list
        A list of the full path to all files matching filefilter

    """

    matches = []

    for root, dirnames, filenames in os.walk(basedir):

        for filename in fnmatch.filter(filenames, filefilter):
            fname = os.path.join(root, filename)
            fname = os.path.realpath(fname)

            matches.append(fname)

    return matches

def pbverbosity(cmd=None):
    """
    Changes and returns pynbody verbosity.  Works for different versions
    of pynbody.

    **ARGUMENTS**

    cmd
        -If None (default) current verbosity level is returned, nothing is done
        -If 'off', pynbody is silenced
        -If 'on', pynbody verbosity is set on
        -If something else, cmd is assumed to be a verbosity level

    **RETURNS**

    current_verbosity
        pynbody verbosity level before any changes were made

    **EXAMPLES**

    *Toggle pynbody verbosity*

        current_verbosity = pbverbosity('off')
        ...
        do stuff
        ...
        pbverbosity(current_verbosity)
    """

    # -----------------------------
    # Get current verbosity level
    # -----------------------------
    if hasattr(pb, 'logger'):
        # As of v0.30, pynbody uses python's logging to handle verbosity
        logger = True
        current_verbosity = pb.logger.getEffectiveLevel()
        pb.logger.setLevel(logging.ERROR)

    else:

        # For pynbody version < 0.3, verbosity is handled in the config
        logger = False
        current_verbosity = pb.config['verbose']

    # -----------------------------
    # Change verbosity
    # -----------------------------
    if cmd is None:
        # Don't change verbosity.  just return the current verbosity
        pass

    elif cmd == 'off':
        # Toggle verbosity off
        if logger:

            pb.logger.setLevel(logging.ERROR)

        else:

            pb.config['verbose'] = False

    elif cmd == 'on':
        # Toggle verbosity on
        if logger:

            pb.logger.setLevel(logging.DEBUG)

        else:

            pb.config['verbose'] = True

    else:
        # Set verbosity to the verbosity level specified by cmd
        if logger:

            pb.logger.setLevel(cmd)

        else:

            pb.config['verbose'] = cmd

    # Return the verbosity level before any changes were made
    return current_verbosity
    
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