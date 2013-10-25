##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import os, string, gzip, logging
import numpy as np
from cs_util import deleterow

# gdal_available = True #GDAL disabled for now, but should work- BHM 01/04/12
# try:
    # from osgeo import gdal_array, gdal
    # from osgeo.gdalconst import *
    # #print 'GDAL AVAILABLE'
# except ImportError:
# gdal_available = False

# Disable GDAL as it is error-prone for some cases for now. VS - 4/5/09
# below defines null packages for PyDev IDE to ignore missing imports
gdal = None
gdal_array = None

class CSIO:
    gdal_available = False
    
    @staticmethod
    def _check_file_exixts(filename):
        if not os.path.isfile(filename):
            raise RuntimeError('File "'  + filename + '" does not exist')

    @staticmethod
    def _read_header(filename):
        """Reads header for ASCII grids (standard input) or numpy arrays (used for faster read/write when calling Circuitscape from ArcGIS python code)."""
        CSIO._check_file_exixts(filename)
        file_base, file_extension = os.path.splitext(filename)
        if file_extension == '.npy': #numpy array will have an associated header file
            filename = file_base + '.hdr'
        
        with open(filename, 'r') as f:
            try:
                ncols = string.split(f.readline())[1]
            except ValueError:
                raise  RuntimeError('Unable to read ASCII grid: "'  + filename + '". If file is a text list, please use .txt extension.')
            ncols = int(ncols)
            nrows = string.split(f.readline())[1]
            nrows = int(nrows)
            xllcorner = string.split(f.readline())[1]
            xllcorner = float(xllcorner)
            yllcorner = string.split(f.readline())[1]
            yllcorner = float(yllcorner)
            cellsize = string.split(f.readline())[1]
            cellsize = float(cellsize)
           
            try:
                [_ign, nodata] = string.split(f.readline())
                try:
                    nodata = int(nodata)
                except ValueError:
                    nodata = float(nodata)
            except ValueError:
                nodata=False
    
        # print 'header',ncols, nrows, xllcorner, yllcorner, cellsize, nodata 
        return ncols, nrows, xllcorner, yllcorner, cellsize, nodata 


    @staticmethod
    def _reader(filename, data_type):
        """Reads rasters saved as ASCII grids or numpy arrays into Circuitscape."""
        CSIO._check_file_exixts(filename)
        (ncols, nrows, _xllcorner, _yllcorner, _cellsize, nodata) = CSIO._read_header(filename)
    
        _file_base, file_extension = os.path.splitext(filename)     
        if file_extension == '.npy': 
            pmap = np.load(filename, mmap_mode=None)
            pmap = pmap.astype('float64')
            
        elif CSIO.gdal_available == True:
            pmap = np.float64(gdal_array.LoadFile(filename))  
            if nodata != False:    
                pmap = np.where(pmap==nodata, -9999, pmap)
        else:
            if nodata == False:
                pmap = np.loadtxt(filename, skiprows=5, dtype=data_type)
            else:
                pmap = np.loadtxt(filename, skiprows=6, dtype=data_type)
                pmap = np.where(pmap==nodata, -9999, pmap)
        
        if nrows == 1:
            temp = np.zeros((1,pmap.size))
            temp[0,:] = pmap
            pmap = temp
            
        if ncols == 1:
            temp = np.zeros((pmap.size,1))
            temp[:,0] = pmap
            pmap = temp       
      
        return pmap


    @staticmethod    
    def _writer(file_name, data, state, compress):  
        """Writes rasters to ASCII grid or numpy formats."""     
        _out_base, out_extn = os.path.splitext(file_name) 
        
        if out_extn == '.npy': # Data came in as numpy array, so write same.
            np.save(file_name, data)
            return
            
        if CSIO.gdal_available == True:
            fmt = "MEM"
            driver = gdal.GetDriverByName(fmt)
            dst_ds = driver.Create(file_name, len(data[0]), len(data), 1, gdal.GDT_Float32)
    
            ull = state.yllcorner +  state.cellsize * len(data)
            dst_ds.SetGeoTransform([state.xllcorner,  # left x
                                 state.cellsize,   # w-e pixel resolution
                                 0,                   # rotation
                                 ull,                 # upper left corner y
                                 0,                   # rotation
                                 state.cellsize])   # n-s pixel resolution
                                 
       
            dst_ds.GetRasterBand(1).WriteArray(data)
            fmt = 'AAIGrid'
            driver = gdal.GetDriverByName(fmt)
#             dst_ds_new = driver.CreateCopy(file_name, dst_ds) #STILL GETTING LEADING SPACES.
#             dst_ds = None
            
        else:
            f = gzip.open(file_name+'.gz', 'w') if compress else open(file_name, 'w')    
            f.write('ncols         ' + str(state.ncols) + '\n')
            f.write('nrows         ' + str(state.nrows) + '\n')
            f.write('xllcorner     ' + str(state.xllcorner) + '\n')
            f.write('yllcorner     ' + str(state.yllcorner) + '\n')
            f.write('cellsize      ' + str(state.cellsize) + '\n')
            f.write('NODATA_value  ' + str(state.nodata) + '\n')
            
            delimiter = ''
            fmt = ['%.6f ']*state.ncols 
            fmt = delimiter.join(fmt)
            for row in data:
                f.write(fmt % tuple(row) + '\n')
    
            f.close()

    @staticmethod
    def load_graph(filename):
        """Returns data for arbitrary graph or focal node list from file."""
        CSIO._check_file_exixts(filename)
        try:    
            graph_object = np.loadtxt(filename, dtype='Float64', comments='#') 
        except:
            try:
                graph_object = np.loadtxt(filename, dtype='Float64', comments='#', delimiter=',')
            except:
                raise RuntimeError('Error reading file "' + filename + '". Please check file format.')
        return graph_object


    @staticmethod
    def read_point_strengths(filename):
        """Reads list of variable source strengths from disk.
        
        This code also used for reading file for reclassifying input data.
        """
        CSIO._check_file_exixts(filename)
        try:
            point_strengths = np.loadtxt(filename)
        except ValueError:
            raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
           
        i = np.argsort(point_strengths[:,0])
        return point_strengths[i]

    @staticmethod
    def _read_txt_list(filename, data_type, habitat_size):
        CSIO._check_file_exixts(filename)
        try:
            points = np.loadtxt(filename)
        except ValueError:
            raise RuntimeError('File "'  + filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                
        
        points_rc = np.zeros(points.shape, dtype=data_type)
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

    @staticmethod            
    def _resample_map(filename, reading_mask, resample_to, header, pmap):
        """Code to crudely resample input raster if raster headers don't match (i.e. different extents or cell sizes used)."""  
        try:
            (_ncols, nrows, xllcorner, yllcorner, cellsize, _nodata) = header
            
            if reading_mask==True:
                pmap = np.where(pmap>0, 1, 0)      
                pmap = 1-pmap #now zeros are areas to keep, ones are masked out
                
            (rows, cols) = np.where(pmap > 0)
            
            xcoords = (cols + 0.5) * cellsize + xllcorner
            ycoords = (nrows - rows - 0.5) * cellsize + yllcorner
            
            values = np.zeros(rows.shape, dtype='int32')
            for i in range(0, rows.size):
                values[i] = pmap[rows[i], cols[i]]
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
    
            #From map_rc to pmap:
            pmap = np.zeros((resample_to.nrows, resample_to.ncols), int)
            pmap[map_rc[:,1], map_rc[:,2]] = map_rc[:,0] 
            if reading_mask == True:
                pmap = 1-pmap #now zeros are areas to mask out, ones are kept
            
        except: 
            raise RuntimeError('Error resampling focal node, mask, or short-circuit region locations to match habitat pmap cell size and extent.  We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.')   
    
        return pmap

    @staticmethod
    def read_poly_map(filename, reading_mask, nodata_as, habitat_size, resample, file_type, data_type):
        """Reads raster maps for short-circuit regions (aka polygon poly_map), focal nodes, masks, current sources or grounds from disk."""  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = header = CSIO._read_header(filename)
        
        msg_resample = '%s raster has different %s than habitat raster. Circuitscape will try to crudely resample the raster. We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
        msg_no_resample = '%s raster must have same %s as habitat raster'
    
        poly_map = CSIO._reader(filename, data_type)
        if nodata_as != None:
            poly_map = np.where(poly_map==nodata, nodata_as, poly_map)
    
        if cellsize != habitat_size.cellsize:
            if resample:
                logging.warning(msg_resample % (file_type, "cell size",))
                poly_map = CSIO._resample_map(filename, reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(msg_no_resample%(file_type, "cell_size"))            
        elif ncols != habitat_size.ncols:
            if resample:
                logging.warning(msg_resample % (file_type, "number of columns",))
                poly_map = CSIO._resample_map(filename, reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(msg_no_resample%(file_type, "number of columns"))            
        elif nrows != habitat_size.nrows:
            if resample:
                logging.warning(msg_resample % (file_type, "number of rows",))
                poly_map = CSIO._resample_map(filename, reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(msg_no_resample%(file_type, "number of rows"))            
        elif xllcorner != habitat_size.xllcorner:
            if resample:
                logging.warning(msg_resample % (file_type, "xllcorner",))
                poly_map = CSIO._resample_map(filename, reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(msg_no_resample%(file_type, "xllcorner"))            
        elif yllcorner != habitat_size.yllcorner:
            if resample:
                logging.warning(msg_resample % (file_type, "yllcorner",))
                poly_map = CSIO._resample_map(filename, reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(msg_no_resample%(file_type, "yllcorner"))            
    
        if reading_mask==True:
            poly_map = np.where(poly_map < 0, 0, poly_map)        
    
        return poly_map

    @staticmethod
    def read_point_map(filename, file_type, habitat_size):
        """Reads map or text list of focal nodes from disk.
        
        File extension is used to determine whether format is ascii grid, numpy array, or text list.
        """
        _base, extension = os.path.splitext(filename)
        
        if extension not in [".txt", ".asc", ".npy"]:
            raise RuntimeError('%s file must have a .txt, .asc or .npy extension'%(file_type,))
        
        if extension == ".txt":
            points_rc = CSIO._read_txt_list(filename, 'int32', habitat_size)
        elif extension == ".asc" or extension == ".npy": # We use Numpy format for quickly passing grids between ArcGIS and Circuitscape.
            point_map = CSIO.read_poly_map(filename, False, 0, habitat_size, True, file_type, 'int32')
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
    

    @staticmethod    
    def read_source_and_ground_maps(source_filename, ground_filename, habitat_size, options): 
        """Reads srouce and ground raster maps from disk."""
        #print("reading source and ground maps:\n\t[%s]\n\t[%s]"%(source_filename,ground_filename))  
        #FIXME: reader does not currently handle infinite inputs for ground conductances.
        extension = os.path.splitext(source_filename)[1]
        if extension not in [".txt", ".asc"]:
            raise RuntimeError('Current source files must have a .txt or .asc extension')
        
        if extension == ".txt":  
            #FIXME: probably want to roll code used for reading source, ground and point text files into single utility
            sources_rc = CSIO._read_txt_list(source_filename, 'int32', habitat_size)
            source_map = np.zeros((habitat_size.nrows, habitat_size.ncols), dtype='float64')
            source_map[sources_rc[:,1], sources_rc[:,2]] = sources_rc[:,0]
        elif extension=='.asc':
            source_map = CSIO.read_poly_map(source_filename, False, 0, habitat_size, False, "Current source", 'float64')
            source_map = np.where(source_map==-9999, 0, source_map)
    
        if options.use_unit_currents == True:
            source_map = np.where(source_map, 1, source_map)
    
    
        _base, extension = os.path.splitext(ground_filename)
        if extension not in [".txt", ".asc"]:
            raise RuntimeError('Ground files must have a .txt or .asc extension')
        
        if extension == ".txt":
            grounds_rc = CSIO._read_txt_list(ground_filename, 'int32', habitat_size)
            ground_map_raw = -9999 * np.ones((habitat_size.nrows, habitat_size.ncols), dtype = 'float64')
            ground_map_raw[grounds_rc[:,1], grounds_rc[:,2]] = grounds_rc[:,0]
        elif extension=='.asc':
            ground_map_raw = CSIO.read_poly_map(ground_filename, False, None, habitat_size, False, "Ground", 'float64')
    
        if options.ground_file_is_resistances==True:
            ground_map = 1 / ground_map_raw
            ground_map = np.where(ground_map_raw == -9999, 0, ground_map)
        else:
            ground_map = np.where(ground_map_raw == -9999, 0, ground_map_raw)
        
        if options.use_direct_grounds == True:
            ground_map = np.where(ground_map, np.Inf, 0)
    
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

    @staticmethod
    def read_included_pairs(filename):
        """Reads matrix denoting node pairs to include/exclude from calculations.
        
        FIXME: matrices are an inconvenient way for users to specify pairs.  Using a 
        2- or 3-column format would be easier.
        """
        CSIO._check_file_exixts(filename)
        
        try:
            with open(filename, 'r') as f:
                minval = string.split(f.readline())[1]
                maxval = string.split(f.readline())[1]
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

    @staticmethod
    def write_aaigrid(grid_type, fileadd, data, options, state):
        """Writes ASCII grid.  This is main raster output format for Circuitscape."""  
        if grid_type == 'voltmap':
            if options.write_volt_maps == False: 
                return
            if options.set_null_voltages_to_nodata == True:
                data = np.where(state.g_map == 0, state.nodata, data)
        elif grid_type in ['curmap', 'cum_curmap', 'max_curmap']:
            if options.write_cur_maps == False: 
                return
            if options.set_null_currents_to_nodata == True:
                data = np.where(state.g_map == 0, state.nodata, data)
        else:
            return
    
        out_base = os.path.splitext(options.output_file)[0]
        inp_extn = os.path.splitext(options.habitat_file)[1]
        out_file = out_base + '_' + grid_type + fileadd + ('.npy' if (inp_extn == '.npy') else '.asc')
        CSIO._writer(out_file, data, state, options.compress_grids)

    @staticmethod
    def read_cell_map(habitat_map, is_resistances, reclass_file, state):
        """Reads resistance or conductance raster into memory, converts former to conductance format."""  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = CSIO._read_header(habitat_map)
        state.ncols = ncols
        state.nrows = nrows
        state.xllcorner = xllcorner
        state.yllcorner = yllcorner
        state.cellsize = cellsize
        state.nodata = -9999 if (nodata == False) else nodata 
    
        cell_map = CSIO._reader(habitat_map, 'float64')
    
        # Reclassification code
        if reclass_file != None:
            try:
                reclass_table = CSIO.read_point_strengths(reclass_file)
            except:
                raise RuntimeError('Error reading reclass table')
            for i in range (0, reclass_table.shape[0]):
                cell_map = np.where(cell_map==reclass_table[i,0], reclass_table[i,1], cell_map)
            logging.info('Reclassified habitat map using %s'%(reclass_file,))
        
        if is_resistances == True:
            zeros_in_resistance_map = (np.where(cell_map==0, 1, 0)).sum() > 0
            if zeros_in_resistance_map == True: 
                #FIXME: Should be easy to accomodate zeros in resistance map, just treat them like polygons.
                raise RuntimeError('Error: zero resistance values are not currently supported for habitat maps.  Use a short-circuit region file instead.')
            g_map = 1 / cell_map  
            g_map = np.where(cell_map == -9999, 0, g_map)
        else:
            g_map = np.where(cell_map == -9999, 0, cell_map)    
        state.g_map = np.where(g_map < 0, 0, g_map)    

    @staticmethod
    def save_incomplete_resistances(outfile_template, resistances):
        """Saves resistances from ongoing calculations.  
        
        Helpful for debugging or recovering partial results after crash or user abort.
        """  
        out_base, out_ext = os.path.splitext(outfile_template)
        out_file = out_base + '_resistances_incomplete' + out_ext
        np.savetxt(out_file, resistances)
        return

    @staticmethod
    def write_resistances_3columns(outfile_template, resistances_3columns):    
        """Writes effective resistances to disk in 3 column format."""  
        out_base, out_ext = os.path.splitext(outfile_template)
        out_file = out_base + '_resistances_3columns' + out_ext
        np.savetxt(out_file, resistances_3columns)
        return         

    @staticmethod
    def write_resistances_one_to_all(outfile_template, resistances, string, scenario):
        """Saves effective resistances from one-to-all calculations to disk."""
        out_base, out_extn = os.path.splitext(outfile_template)
        out_file = out_base + ('_resistances' if (scenario == 'one-to-all') else '_results') + string + out_extn
        np.savetxt(out_file, resistances) 
        
        #remove partial result file        
        if string=='':
            old_file = out_base + ('_resistances_incomplete' if (scenario == 'one-to-all') else '_results_incomplete') + string + out_extn
            try:
                os.remove(old_file)
            except:
                pass 
        return

    @staticmethod
    def write_resistances(outfile_template, resistances, resistances_3columns):
        """Writes resistance file to disk."""
        out_base, out_extn = os.path.splitext(outfile_template)
        out_file = out_base + '_resistances' + out_extn
        np.savetxt(out_file, resistances)
        
        CSIO.write_resistances_3columns(outfile_template, resistances_3columns)
        
        #remove partial result file        
        old_file = out_base + '_resistances_incomplete' + out_extn
        try:
            os.remove(old_file)
        except:
            pass 
    
    @staticmethod
    def write_currents(outfile_template, branch_currents, node_currents, fileadd):
        """Writes currents from network operations.
        
           Inputs are arrays with node names.
        """
        out_base, _out_ext = os.path.splitext(outfile_template)
        
        if branch_currents!=None:
            filename = out_base + '_branch_currents_' + fileadd + '.txt'
            np.savetxt(filename, branch_currents)
        filename = out_base + '_node_currents_' + fileadd + '.txt'
        np.savetxt(filename, node_currents)      

    @staticmethod
    def write_voltages(outfile_template, voltages, node_names, fileadd):
        """Saves voltage values from solving arbitrary graphs to disk."""  
        output_voltages = np.zeros((len(voltages),2), dtype='float64')
        if node_names != None:
            output_voltages[:,0] = node_names[:]
        output_voltages[:,1] = voltages[:]

        out_base, _out_extn = os.path.splitext(outfile_template)
        filename = out_base + '_voltages_' + fileadd + '.txt'
        np.savetxt(filename, output_voltages)      

