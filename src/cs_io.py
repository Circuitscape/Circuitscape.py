##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import os, string, gzip, time
import numpy as np
from cs_util import *
#from scipy import sparse

# gdal_available = True #GDAL disabled for now, but should work- BHM 01/04/12
# try:
    # from osgeo import gdal_array, gdal
    # from osgeo.gdalconst import *
    # #print 'GDAL AVAILABLE'
# except ImportError:
   # gdal_available = False

# Disable GDAL as it is error-prone for some cases for now. VS - 4/5/09
gdal_available = False

# TODO: make IO object
def read_header(filename):
    """Reads header for ASCII grids (standard input) or numpy arrays (used
    for faster read/write when calling Circuitscape from ArcGIS python code).
    """    
    if not os.path.isfile(filename):
        raise RuntimeError('File "'  + filename + '" does not exist')
    file_base, file_extension = os.path.splitext(filename)
    if file_extension == '.npy': #numpy array will have an associated header file
        filename = file_base + '.hdr'
    
    with open(filename, 'r') as f:
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
                nodata = int(nodata)
            except ValueError:
                nodata = float(nodata)
        except ValueError:
            nodata=False

    # print 'header',ncols, nrows, xllcorner, yllcorner, cellsize, nodata 
    return ncols, nrows, xllcorner, yllcorner, cellsize, nodata 


def reader(filename, type):
    """Reads rasters saved as ASCII grids or numpy arrays into Circuitscape."""
    if not os.path.isfile(filename):      
        raise RuntimeError('File "'  + filename + '" does not exist')
    (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)

    file_base, file_extension = os.path.splitext(filename)     
    if file_extension == '.npy': 
        map = np.load(filename, mmap_mode=None)
        map = map.astype('float64')
        
    elif gdal_available == True:
        map = np.float64(gdal_array.LoadFile(filename))  
        if nodata != False:    
            #map = np.where(map==nodata, -9999, map)
            map[np.where(map==nodata)] = -9999
    else:
        if nodata == False:
            map = np.loadtxt(filename, skiprows=5, dtype=type)
        else:
            map = np.loadtxt(filename, skiprows=6, dtype=type)
            map[np.where(map==nodata)] = -9999
            #map = np.where(map==nodata, -9999, map)
    
    if nrows == 1:
        temp = np.zeros((1,map.size))
        temp[0,:] = map
        map = temp
        
    if ncols == 1:
        temp = np.zeros((map.size,1))
        temp[:,0] = map
        map = temp       
  
    return map

    
def writer(file, data, state, compress):  
    """Writes rasters to ASCII grid or numpy formats."""     
    outputBase, outputExtension = os.path.splitext(file) 
    
    if outputExtension == '.npy': # Data came in as numpy array, so write same.
        np.save(file, data)
        return
        
    if gdal_available == True:
        format = "MEM"
        driver = gdal.GetDriverByName( format )      
        dst_ds = driver.Create( file, len(data[0]), len(data),1,gdal.GDT_Float32)

        ull=state.yllcorner+(state.cellsize)*(len(data))
        dst_ds.SetGeoTransform([state.xllcorner,  # left x
                             state.cellsize,   # w-e pixel resolution
                             0,                   # rotation
                             ull,                 # upper left corner y
                             0,                   # rotation
                             state.cellsize])   # n-s pixel resolution
                             
   
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

        f.write('ncols         ' + str(state.ncols) + '\n')
        f.write('nrows         ' + str(state.nrows) + '\n')
        f.write('xllcorner     ' + str(state.xllcorner) + '\n')
        f.write('yllcorner     ' + str(state.yllcorner) + '\n')
        f.write('cellsize      ' + str(state.cellsize) + '\n')
        f.write('NODATA_value  ' + str(state.nodata) + '\n')
        
        delimiter = ''
        fmt = ['%.6f ']*state.ncols 
        format = delimiter.join(fmt)
        for row in data:
            f.write(format % tuple(row) + '\n')

        f.close()

def load_graph(filename):
    """Returns data for arbitrary graph or focal node list from file."""
    if not os.path.isfile(filename):
        raise RuntimeError('File "'  + filename + '" does not exist')

    try:    
        graph_object = np.loadtxt(filename, dtype = 'Float64', comments='#') 
    except:
        try:
            graph_object = np.loadtxt(filename, dtype = 'Float64', comments='#', delimiter=',')
        except:
            raise RuntimeError('Error reading file "' + filename + '". Please check file format.')
    return graph_object


def read_point_strengths(filename):
    """
    Reads list of variable source strengths from disk.
    This code also used for reading file for reclassifying input data.
    """
    if not os.path.isfile(filename):
        raise RuntimeError('File "'  + filename + '" does not exist')
    
    try:
        point_strengths = np.loadtxt(filename)
    except ValueError:
        raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
       
    i = argsort(point_strengths[:,0])
    return point_strengths[i]


#     def writeGraph(self,filename,graph,nodeNames):
#         """Save graph to disk in 3-column format."""  
#         graphNcol = self.convert_graph_to_3_col(graph,nodeNames)
#         savetxt(filename,graphNcol)
#         return
            
def resample_map(filename, reading_mask, resample_to, header, map):
    """Code to crudely resample input raster if raster headers don't match (i.e. different extents or cell sizes used)."""  
    try:
        print("resample_map called")
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = header
        
        if reading_mask==True:
            map = np.where(map>0, 1, 0)      
            map = 1-map #now zeros are areas to keep, ones are masked out
            
        (rows, cols) = np.where(map > 0)
        
        xcoords = (cols + 0.5) * cellsize + xllcorner
        ycoords = (nrows - rows - 0.5) * cellsize + yllcorner
        
        values = np.zeros(rows.shape, dtype='int32')
        for i in range(0, rows.size):
            values[i] = map[rows[i], cols[i]]
        map_coords = c_[values, xcoords, ycoords]

        i = argsort(map_coords[:,0])
        map_coords = map_coords[i]
        
        #From map_coords to map_rc:
        map_rc = np.zeros(map_coords.shape, dtype='int32')
        map_rc[:,0] = map_coords[:,0]
        map_rc[:,1] = ceil((resample_to.nrows - (map_coords[:,2] - resample_to.yllcorner) / resample_to.cellsize)) - 1
        map_rc[:,2] = ceil(((map_coords[:,1] - resample_to.xllcorner) / resample_to.cellsize)) - 1
        i = argsort(map_rc[:,0])
        map_rc = map_rc[i]

        rows = map_rc[:,1]
        delrows = asarray(np.where(rows < 0))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
            
        rows = map_rc[:,1] 
        delrows = asarray(np.where(rows > resample_to.nrows - 1))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
            
        cols = map_rc[:,2]
        delrows = asarray(np.where(cols < 0))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
            
        cols = map_rc[:,2]
        delrows = asarray(np.where(cols > resample_to.ncols - 1))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
        del delrows
        del delrows2

        #From map_rc to map:
        map = np.zeros((resample_to.nrows, resample_to.ncols), int)
        map[map_rc[:,1], map_rc[:,2]] = map_rc[:,0] 
        if reading_mask == True:
            map = 1-map #now zeros are areas to mask out, ones are kept
        
    except: 
        raise RuntimeError('Error resampling focal node, mask, or short-circuit region locations to match habitat map cell size and extent.  We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.')   

    return map

def read_poly_map(filename, reading_mask, resample_to, file_type):
    """Reads short-circuit region map (aka polygon map) from disk."""  
    (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = header = read_header(filename)
    
    msg = '\n********\nWarning: %s raster has different \n%s than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'

    map = reader(filename, 'int32')
    map[np.where(map==nodata)] = 0        
    #map = np.where(map == nodata, 0, map)

    if cellsize!= resample_to.cellsize:
        print(msg % (file_type, "cell size",))
        map = resample_map(filename, reading_mask, resample_to, header, map)
    elif ncols!= resample_to.ncols:
        print(msg % (file_type, "number of columns",))
        map = resample_map(filename, reading_mask, resample_to, header, map)
    elif nrows!= resample_to.nrows:
        print(msg % (file_type, "number of rows",))
        map = resample_map(filename, reading_mask, resample_to, header, map)
    elif xllcorner!= resample_to.xllcorner:
        print(msg % (file_type, "xllcorner",))
        map = resample_map(filename, reading_mask, resample_to, header, map)
    elif yllcorner!= resample_to.yllcorner:
        print(msg % (file_type, "yllcorner",))
        map = resample_map(filename, reading_mask, resample_to, header, map)

    if reading_mask==True:
        map[np.where(map < 0)] = 0
        #map = np.where(map < 0, 0, map)        

    return map
