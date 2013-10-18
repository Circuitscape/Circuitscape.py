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

def check_file_exixts(filename):
    if not os.path.isfile(filename):
        raise RuntimeError('File "'  + filename + '" does not exist')

# TODO: make IO object
def read_header(filename):
    """Reads header for ASCII grids (standard input) or numpy arrays (used for faster read/write when calling Circuitscape from ArcGIS python code)."""
    check_file_exixts(filename)
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
    check_file_exixts(filename)
    (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = read_header(filename)

    file_base, file_extension = os.path.splitext(filename)     
    if file_extension == '.npy': 
        map = np.load(filename, mmap_mode=None)
        map = map.astype('float64')
        
    elif gdal_available == True:
        map = np.float64(gdal_array.LoadFile(filename))  
        if nodata != False:    
            map = np.where(map==nodata, -9999, map)
    else:
        if nodata == False:
            map = np.loadtxt(filename, skiprows=5, dtype=type)
        else:
            map = np.loadtxt(filename, skiprows=6, dtype=type)
            map = np.where(map==nodata, -9999, map)
    
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
    check_file_exixts(filename)
    try:    
        graph_object = np.loadtxt(filename, dtype = 'Float64', comments='#') 
    except:
        try:
            graph_object = np.loadtxt(filename, dtype = 'Float64', comments='#', delimiter=',')
        except:
            raise RuntimeError('Error reading file "' + filename + '". Please check file format.')
    return graph_object


def read_point_strengths(filename):
    """Reads list of variable source strengths from disk.
    
    This code also used for reading file for reclassifying input data.
    """
    check_file_exixts(filename)
    try:
        point_strengths = np.loadtxt(filename)
    except ValueError:
        raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
       
    i = np.argsort(point_strengths[:,0])
    return point_strengths[i]

def read_txt_list(filename, type, habitat_size):
    check_file_exixts(filename)
    try:
        points = np.loadtxt(filename)
    except ValueError:
        raise RuntimeError('File "'  + filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                
    
    points_rc = np.zeros(points.shape, dtype=type)
    try:
        points_rc[:,0] = points[:,0]
        points_rc[:,1] = np.ceil(habitat_size.nrows - (points[:,2] - habitat_size.yllcorner) / habitat_size.cellsize) - 1
        points_rc[:,2] = np.ceil((points[:,1] - habitat_size.xllcorner) / habitat_size.cellsize) - 1
    except IndexError:
        raise RuntimeError('Error extracting locations from .txt file. Please check file format.')
    return points_rc
    

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
        map_coords = np.c_[values, xcoords, ycoords]

        i = np.argsort(map_coords[:,0])
        map_coords = map_coords[i]
        
        #From map_coords to map_rc:
        map_rc = np.zeros(map_coords.shape, dtype='int32')
        map_rc[:,0] = map_coords[:,0]
        map_rc[:,1] = np.ceil((resample_to.nrows - (map_coords[:,2] - resample_to.yllcorner) / resample_to.cellsize)) - 1
        map_rc[:,2] = np.ceil(((map_coords[:,1] - resample_to.xllcorner) / resample_to.cellsize)) - 1
        i = np.argsort(map_rc[:,0])
        map_rc = map_rc[i]

        rows = map_rc[:,1]
        delrows = np.asarray(np.where(rows < 0))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
            
        rows = map_rc[:,1] 
        delrows = np.asarray(np.where(rows > resample_to.nrows - 1))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
            
        cols = map_rc[:,2]
        delrows = np.asarray(np.where(cols < 0))
        delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
        delrows2[:] = delrows[:]
        if delrows2 != []:
            map_rc = deleterow(map_rc, delrows2)
            
        cols = map_rc[:,2]
        delrows = np.asarray(np.where(cols > resample_to.ncols - 1))
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

def read_poly_map(filename, reading_mask, nodata_as, habitat_size, resample, file_type, data_type):
    """Reads short-circuit region map (aka polygon map) from disk."""  
    (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = header = read_header(filename)
    
    msg_resample = '\n********\nWarning: %s raster has different \n%s than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
    msg_no_resample = '%s raster must have same %s as habitat raster'

    map = reader(filename, data_type)
    if nodata_as != None:
        map = np.where(map==nodata, nodata_as, map)

    if cellsize!= habitat_size.cellsize:
        if resample:
            print(msg_resample % (file_type, "cell size",))
            map = resample_map(filename, reading_mask, habitat_size, header, map)
        else:
            raise RuntimeError(msg_no_resample%(file_type, "cell_size"))            
    elif ncols!= habitat_size.ncols:
        if resample:
            print(msg_resample % (file_type, "number of columns",))
            map = resample_map(filename, reading_mask, habitat_size, header, map)
        else:
            raise RuntimeError(msg_no_resample%(file_type, "number of columns"))            
    elif nrows!= habitat_size.nrows:
        if resample:
            print(msg_resample % (file_type, "number of rows",))
            map = resample_map(filename, reading_mask, habitat_size, header, map)
        else:
            raise RuntimeError(msg_no_resample%(file_type, "number of rows"))            
    elif xllcorner!= habitat_size.xllcorner:
        if resample:
            print(msg_resample % (file_type, "xllcorner",))
            map = resample_map(filename, reading_mask, habitat_size, header, map)
        else:
            raise RuntimeError(msg_no_resample%(file_type, "xllcorner"))            
    elif yllcorner!= habitat_size.yllcorner:
        if resample:
            print(msg_resample % (file_type, "yllcorner",))
            map = resample_map(filename, reading_mask, habitat_size, header, map)
        else:
            raise RuntimeError(msg_no_resample%(file_type, "yllcorner"))            

    if reading_mask==True:
        map = np.where(map < 0, 0, map)        

    return map

def read_point_map(filename, file_type, habitat_size):
    """Reads map or text list of focal nodes from disk.
    
    File extension is used to determine whether format is ascii grid, numpy array, or text list.
    """
    base, extension = os.path.splitext(filename)
    
    if extension not in [".txt", ".asc", ".npy"]:
        raise RuntimeError('%s file must have a .txt, .asc or .npy extension'%(file_type,))
    
    if extension == ".txt":
        points_rc = read_txt_list(filename, 'int32', habitat_size)
    elif extension == ".asc" or extension == ".npy": # We use Numpy format for quickly passing grids between ArcGIS and Circuitscape.
        point_map = read_poly_map(filename, False, 0, habitat_size, True, file_type, 'int32')
        (rows, cols) = np.where(point_map > 0)

        values = np.zeros(rows.shape, dtype='int32') 
        for i in range(0, rows.size):
            values[i] = point_map[rows[i], cols[i]]
        points_rc = np.c_[values, rows, cols]
    else:
        raise RuntimeError('Focal node file must have a .txt or .asc extension')

    try:            
        i = np.argsort(points_rc[:,0])
        points_rc = points_rc[i]
    except IndexError:
        raise RuntimeError('Error extracting focal node locations. Please check file format.')                
    
    # Check to make sure points fall within cellmap
    if (np.min(points_rc[:,1]) < 0) or (np.min(points_rc[:,2]) < 0) \
        or (np.max(points_rc[:,1]) > (habitat_size.nrows - 1)) \
        or (np.max(points_rc[:,2]) > (habitat_size.ncols - 1)):
        raise RuntimeError('At least one focal node location falls outside of habitat map')
    
    if (np.unique(np.asarray(points_rc[:,0]))).shape[0] < 2:
        raise RuntimeError('Less than two valid focal nodes found. Please check focal node location file.')                    
    
    return points_rc
    
    
def read_source_and_ground_maps(source_filename, ground_filename, habitat_size, options): 
    """Reads srouce and ground raster maps from disk."""
    #print("reading source and ground maps:\n\t[%s]\n\t[%s]"%(source_filename,ground_filename))  
    #FIXME: reader does not currently handle infinite inputs for ground conductances.
    base, extension = os.path.splitext(source_filename)
    if extension not in [".txt", ".asc"]:
        raise RuntimeError('Current source files must have a .txt or .asc extension')
    
    if extension == ".txt":  
        #FIXME: probably want to roll code used for reading source, ground and point text files into single utility
        sources_rc = read_txt_list(source_filename, 'int32', habitat_size)
        source_map = np.zeros((habitat_size.nrows, habitat_size.ncols), dtype='float64')
        source_map[sources_rc[:,1], sources_rc[:,2]] = sources_rc[:,0]
    elif extension=='.asc':
        source_map = read_poly_map(source_filename, False, 0, habitat_size, False, "Current source", 'float64')
        source_map = np.where(source_map==-9999, 0, source_map)

    if options.use_unit_currents == True:
        source_map = np.where(source_map, 1, source_map)


    base, extension = os.path.splitext(ground_filename)
    if extension not in [".txt", ".asc"]:
        raise RuntimeError('Ground files must have a .txt or .asc extension')
    
    if extension == ".txt":
        grounds_rc = read_txt_list(ground_filename, 'int32', habitat_size)
        ground_map_raw = -9999 * np.ones((habitat_size.nrows, habitat_size.ncols), dtype = 'float64')
        ground_map_raw[grounds_rc[:,1], grounds_rc[:,2]] = grounds_rc[:,0]
    elif extension=='.asc':
        ground_map_raw = read_poly_map(ground_filename, False, None, habitat_size, False, "Ground", 'float64')

    if options.ground_file_is_resistances==True:
        ground_map = 1 / ground_map_raw
        ground_map = np.where(ground_map_raw == -9999, 0, ground_map)
    else:
        ground_map = np.where(ground_map_raw == -9999, 0, ground_map_raw)
    
    if options.use_direct_grounds == True:
        ground_map = np.where(ground_map, Inf, 0)

    conflicts = np.logical_and(source_map, ground_map)
    if options.remove_src_or_gnd in ['rmvsrc', 'rmvall']:
        source_map = np.where(conflicts, 0, source_map)
    if options.remove_src_or_gnd in ['rmvgnd', 'rmvall']:
        ground_map = np.where(conflicts, 0, ground_map)

    if np.size(np.where(source_map)) == 0:
        raise RuntimeError('No valid sources detected. Please check source file') 
    if np.size(np.where(ground_map)) == 0:
        raise RuntimeError('No valid grounds detected. Please check ground file') 
    return source_map, ground_map

def read_included_pairs(filename):
    """Reads matrix denoting node pairs to include/exclude from calculations.
    
    FIXME: matrices are an inconvenient way for users to specify pairs.  Using a 
    2- or 3-column format would be easier.
    
    """
    check_file_exixts(filename)
    
    try:
        with open(filename, 'r') as f:
            [ign, minval] = string.split(f.readline())
            [ign, maxval] = string.split(f.readline())
            minval = float(minval)
            maxval = float(maxval)
            
        included_pairs = np.loadtxt(filename, skiprows=2, dtype='Float64')
        point_ids = included_pairs[:,0]
        included_pairs = np.where(included_pairs>maxval, 0, included_pairs)
        included_pairs = np.where(included_pairs<minval, 0, 1)             
        included_pairs[:,0] = point_ids
        included_pairs[0,:] = point_ids
        included_pairs[0,0] = -1
        i = np.argsort(included_pairs[:,0])
        included_pairs = included_pairs[i]
        included_pairs = included_pairs.T            
        i = np.argsort(included_pairs[:,0])
        included_pairs = included_pairs[i] 
    except:
        raise RuntimeError('Error reading focal node include/exclude matrix. Please check file format.')                

    return included_pairs
