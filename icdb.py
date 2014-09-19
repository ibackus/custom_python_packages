# -*- coding: utf-8 -*-
"""
An incredibly basic, not very flexible database for keeping track of all
the initial conditions that I create.  This is not user friendly, but mostly
for my personal use

EXAMPLE

import icdb

# Create a new database looking at 2 base directories (and all subdirectories)

# Finds all files matching '*_settings.p'

db = icdb.database('*_settings.p', '~/simulations', \
'/home/ibackus/hyak/gscratch/vsm/ibackus/simulations')

# Query the full paths to IC files:

db.data('filename')

# Query Qmin

db.data('sigma.Qmin')

# Access the settings/information for the first simulation

sim_settings = db.data[0]

# Refresh database

db.refresh()

Created on Mon Sep 15 16:06:30 2014

@author: ibackus
"""
import numpy as np
import ICgen_settings
import os
import fnmatch

class database():
    
    """
    IC database class.  Initialize with:
    
    db = icdb.database(filefilter, dir1, dir2, ...)
    
    filefilter is a string used to filter IC settings files, ie '*_settings.p'
    dir1, dir2 are base directories which are recursively searched to find
    settings files
    """
    
    def __init__(self, filefilter='IC_settings.p', *directories, **kwargs):
        
        print 'Building database...'
        if len(directories) == 0:
            
            directories = ['.']
        print 'Number of base directories: ', len(directories)
        self.dirs = directories
        self.filefilter = filefilter
        self._kwargs = kwargs
        
        matches = []
        
        for directory in directories:
            
            print 'Building for ', directory            
            new_matches = self._scan_dir(directory, filefilter)
            print '    Found {} files'.format(len(new_matches))
            matches += new_matches

        filelist = list(set(matches))
        datalist = []
        
        for i, fname in enumerate(filelist):
            
            a = ICgen_settings.settings()
            a.load(fname)
            a_dir = os.path.split(fname)[0]
            ICname = a.filenames.IC_file_name
            ICname = os.path.join(a_dir, ICname)
            ICname = os.path.realpath(ICname)
            
            if os.path.exists(ICname):
                
                a.dir = os.path.realpath(a_dir)
                a.ICname = ICname
                datalist.append(a)
            
        nfiles = len(datalist)
        print 'Found {} files in total'.format(nfiles)
        data = fancy_array(nfiles)
        data[:] = datalist
        self.data = data
        
    def nfiles(self):
        """
        Returns the current number of files in the database
        """
        
        return len(self.data)
    
    def _scan_dir(self, basedir='.', filefilter='*', kind='IC'):
        
        if kind is 'IC':
            
            matches = []
            
            for root, dirnames, filenames in os.walk(basedir):
                for filename in fnmatch.filter(filenames, filefilter):
                    fname = os.path.join(root, filename)
                    fname = os.path.realpath(fname)
                    
                    matches.append(fname)
                    
            return matches
                    
    def refresh(self):
        """
        Refreshes the database by re-initializing it
        """
        
        self.__init__(self.filefilter, *self.dirs, **self._kwargs)

class fancy_array(np.ndarray):
    
    def __new__(subtype, shape, buffer=None, offset=0, \
    strides=None, order=None):
        
        if not isinstance(shape, int):
            
            raise ValueError, 'Shape must be an integer.  1D arrays only'
            
        dtype = object
        
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, \
        strides, order)
        
        return obj
        
    def __array_finalize__(self, obj):
        
        if obj is None: 
            
            return
            
        #self.info = getattr(obj, 'info', None)
        #self.info = obj.info
        
    def __call__(self, attr):
        
        outarray = np.ndarray(self.shape, dtype=object)
        
        for i, entry in enumerate(self):
            
            if hasattr_nested(entry, attr):
                
                outarray[i] = getattr_nested(entry, attr)
            
            else:
                
                outarray[i] = None
                
        return outarray
            
        
def hasattr_nested(obj, attr):
    """
    Check whether an object contains a specified (possibly nested) attribute
    
    ie:
    
    hasattr_nested(obj, 'x.y.z') will check if obj.x.y.z exists
    """
    
    parts = attr.split('.')
    
    for part in parts:
        
        if hasattr(obj, part):
            obj = getattr(obj, part)
        else:
            return False
    else:
        return True
        
def getattr_nested(obj, attr):
    """
    Get attribute (possibly nested) from object
    
    ie:
    
    a = getattr_nested(obj, 'x.y.z')
    
    is equivalent to a = obj.x.y.z
    """
    
    parts = attr.split('.')
    
    for part in parts:
        
        if hasattr(obj, part):
            obj = getattr(obj, part)
        else:
            return False
    else:
        return obj         
            