##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import os, string, gzip
import numpy
import time

# gdal_available = True #GDAL disabled for now, but should work- BHM 01/04/12
# try:
    # from osgeo import gdal_array, gdal
    # from osgeo.gdalconst import *
    # #print 'GDAL AVAILABLE'
# except ImportError:
   # gdal_available = False

# Disable GDAL as it is error-prone for some cases for now. VS - 4/5/09
gdal_available = False
    
from string import split
from numpy import loadtxt, where
    
from numpy import *
from scipy import sparse

# TODO: make IO object
def read_header(filename):
    """Reads header for ASCII grids (standard input) or numpy arrays (used

    for faster read/write when calling Circuitscape from ArcGIS python code).
    
    """    
    if os.path.isfile(filename)==False:
        raise RuntimeError('File "'  + filename + '" does not exist')
    fileBase, fileExtension = os.path.splitext(filename) 
    if fileExtension == '.npy': #numpy array will have an associated header file
        filename = fileBase + '.hdr'
    
    f = open(filename, 'r')
    try:
        [ign, ncols] = string.split(f.readline())
    except ValueError:
            raise  RuntimeError('Unable to read ASCII grid: "'  + filename + '". If file is a text list, please use .txt extension.')
    ncols = int(ncols)
    [ign, nrows] = string.split(f.readline())
    nrows = int(nrows)
    [ign, xllcorner] = string.split(f.readline())
    xllcorner = float(xllcorner)
    [ign, yllcorner] = string.split(f.readline())
    yllcorner = float(yllcorner)
    [ign, cellsize] = string.split(f.readline())
    cellsize = float(cellsize)
   
    try:
        [ign, nodata] = string.split(f.readline())
        try:
            nodata= int(nodata)
        except ValueError:
            nodata= float(nodata)
    except ValueError:
        nodata=False
  
    f.close()
 
    # print 'header',ncols, nrows, xllcorner, yllcorner, cellsize, nodata 
    return ncols, nrows, xllcorner, yllcorner, cellsize, nodata 

def reader(filename, type):
    """Reads rasters saved as ASCII grids or numpy arrays into Circuitscape."""    
    if os.path.isfile(filename)==False:      
        raise RuntimeError('File "'  + filename + '" does not exist')
    (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)

    fileBase, fileExtension = os.path.splitext(filename)     
    if fileExtension == '.npy': 
        map = numpy.load(filename, mmap_mode=None)
        map = map.astype('float64')
        
    elif gdal_available == True:
        map = numpy.float64(gdal_array.LoadFile(filename))  
        if nodata!=False:    
            map = where(map==nodata, -9999, map)
    else:
        if nodata==False:
            map = loadtxt(filename, skiprows=5, dtype=type)
        else:
            map = loadtxt(filename, skiprows=6, dtype=type)
            map = where(map==nodata, -9999, map)
    if nrows==1:
        temp=numpy.zeros((1,map.size))
        temp[0,:]=map
        map=temp
    if ncols==1:
        temp=numpy.zeros((map.size,1))
        temp[:,0]=map
        map=temp       
  
    return map

    
def writer(file, data, state, compress):  
    """Writes rasters to ASCII grid or numpy formats."""     
    outputBase, outputExtension = os.path.splitext(file) 
    
    if outputExtension == '.npy': # Data came in as numpy array, so write same.
        numpy.save(file, data)
        return
        
    if gdal_available == True:
        format = "MEM"
        driver = gdal.GetDriverByName( format )      
        dst_ds = driver.Create( file, len(data[0]), len(data),1,gdal.GDT_Float32)

        ull=state['yllcorner']+(state['cellsize'])*(len(data))
        dst_ds.SetGeoTransform([state['xllcorner'],  # left x
                             state['cellsize'],   # w-e pixel resolution
                             0,                   # rotation
                             ull,                 # upper left corner y
                             0,                   # rotation
                             state['cellsize']])   # n-s pixel resolution
                             
   
        dst_ds.GetRasterBand(1).WriteArray(data)
        format = 'AAIGrid'
        driver = gdal.GetDriverByName(format)
        dst_ds_new = driver.CreateCopy(file, dst_ds) #STILL GETTING LEADING SPACES.
        dst_ds = None
        
        
    else:
        f = False
        if compress == True:
            file = file + '.gz'
            f = gzip.open(file, 'w')
        else:
            f = open(file, 'w')

        f.write('ncols         ' + str(state['ncols']) + '\n')
        f.write('nrows         ' + str(state['nrows']) + '\n')
        f.write('xllcorner     ' + str(state['xllcorner']) + '\n')
        f.write('yllcorner     ' + str(state['yllcorner']) + '\n')
        f.write('cellsize      ' + str(state['cellsize']) + '\n')
        f.write('NODATA_value  ' + str(state['nodata']) + '\n')
        
        delimiter = ''
        fmt = ['%.6f ']*state['ncols'] 
        format = delimiter.join(fmt)
        for row in data:
            f.write(format % tuple(row) + '\n')

        f.close()

