import os, string, gzip
import numpy as np
from scipy import sparse
from profiler import print_rusage

class CSIO:
    FILE_TYPE_NPY = 1
    FILE_TYPE_AAGRID = 2
    FILE_TYPE_TXTLIST = 3
    FILE_TYPE_INCL_PAIRS_AAGRID = 4
    FILE_TYPE_INCL_PAIRS = 5
    
    FILE_HDR_GZIP = '\x1f\x8b\x08'
    FILE_HDR_NPY = '\x93NUMPY'
    FILE_HDR_AAGRID = 'ncols'
    FILE_HDR_INCL_PAIRS_AAGRID = 'min'
    FILE_HDR_INCL_PAIRS = 'mode'

    MSG_RESAMPLE = '%s raster has different %s than habitat raster. Circuitscape will try to crudely resample the raster. We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
    MSG_NO_RESAMPLE = '%s raster must have same %s as habitat raster'
    
    logger = None
        
    @staticmethod
    def _check_file_exists(filename):
        if not os.path.isfile(filename):
            raise RuntimeError('File "'  + filename + '" does not exist')

    @staticmethod
    def _open_auto_uncompress(filename):
        f = open(filename, 'r')
        #is_compressed = False
        try:
            hdr = f.read(3)
            if hdr.startswith(CSIO.FILE_HDR_GZIP):  # gzip header
                f.close()
                f = gzip.open(filename, 'r')
                #is_compressed = True
            else:
                f.seek(0)
        except:
            f.close()
            f = None
        return f

    @staticmethod
    def _guess_file_type(f):
        is_filename = isinstance(f, str) or isinstance(f, unicode)
        if is_filename:
            f = CSIO._open_auto_uncompress(f)
        
        hdr = f.read(10)

        if is_filename:
            f.close()
        else:
            f.seek(0)
            
        if hdr.startswith(CSIO.FILE_HDR_NPY):
            filetype = CSIO.FILE_TYPE_NPY
        elif hdr.startswith(CSIO.FILE_HDR_AAGRID):
            filetype = CSIO.FILE_TYPE_AAGRID
        elif hdr.startswith(CSIO.FILE_HDR_INCL_PAIRS_AAGRID):
            filetype = CSIO.FILE_TYPE_INCL_PAIRS_AAGRID
        elif hdr.startswith(CSIO.FILE_HDR_INCL_PAIRS):
            filetype = CSIO.FILE_TYPE_INCL_PAIRS
        else:
            filetype = CSIO.FILE_TYPE_TXTLIST
            
        return filetype

    @staticmethod
    def _ascii_grid_read_header(filename):
        """Reads header for ASCII grids (standard input) or numpy arrays (used for faster read/write when calling Circuitscape from ArcGIS python code)."""
        CSIO._check_file_exists(filename)
        
        with CSIO._open_auto_uncompress(filename) as f:
            file_type = CSIO._guess_file_type(f)
            if file_type == CSIO.FILE_TYPE_NPY:
                file_base, _file_extension = os.path.splitext(filename)
                filename = file_base + '.hdr' # numpy array will have an associated header file
                ncols, nrows, xllcorner, yllcorner, cellsize, nodata, _file_type = CSIO._ascii_grid_read_header(filename)
            else:
                try:
                    ncols = int(string.split(f.readline())[1])
                    nrows = int(string.split(f.readline())[1])
                    xllcorner = float(string.split(f.readline())[1])
                    yllcorner = float(string.split(f.readline())[1])
                    cellsize = float(string.split(f.readline())[1])
                except ValueError:
                    raise  RuntimeError('Unable to read ASCII grid: "'  + filename + '".')
               
                try:
                    [_ign, nodata] = string.split(f.readline())
                    try:
                        nodata = int(nodata)
                    except ValueError:
                        nodata = float(nodata)
                except ValueError:
                    nodata = False
    
        return ncols, nrows, xllcorner, yllcorner, cellsize, nodata, file_type 


    @staticmethod
    def _ascii_grid_reader(filename, data_type):
        """Reads rasters saved as ASCII grids or numpy arrays into Circuitscape."""
        (ncols, nrows, _xllcorner, _yllcorner, _cellsize, nodata, filetype) = CSIO._ascii_grid_read_header(filename)
    
        if filetype == CSIO.FILE_TYPE_NPY: 
            pmap = np.load(filename, mmap_mode=None)
            pmap = pmap.astype('float64')
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
    def _ascii_grid_writer(file_name, file_type, data, state, compress):
        """Writes rasters to ASCII grid or numpy formats."""     
        if file_type == CSIO.FILE_TYPE_NPY:
            np.save(file_name, data)
            return            
        else:
            f = gzip.open(file_name+'.gz', 'w') if compress else open(file_name, 'w')    
            f.write('ncols         ' + str(state.ncols) + '\n')
            f.write('nrows         ' + str(state.nrows) + '\n')
            f.write('xllcorner     ' + str(state.xllcorner) + '\n')
            f.write('yllcorner     ' + str(state.yllcorner) + '\n')
            f.write('cellsize      ' + str(state.cellsize) + '\n')
            f.write('NODATA_value  ' + str(state.nodata) + '\n')
             
            delimiter = ''
            fmt = ['%.10g ']*state.ncols 
            fmt = delimiter.join(fmt)
            fmt += '\n'
            for row in data:
                f.write(fmt % tuple(row))
     
            f.close()


    @staticmethod
    def _txt_list_reader(filename, data_type, habitat_size):
        CSIO._check_file_exists(filename)
        try:
            points = np.loadtxt(filename)
        except ValueError:
            raise RuntimeError('File "'  + filename + '" appears to be a text list file, but is not in the correct format.')                
        
        pts_remapped = np.zeros(points.shape, dtype=data_type)
        try:
            pts_remapped[:,0] = points[:,0]
            pts_remapped[:,1] = np.ceil(habitat_size.nrows - (points[:,2] - habitat_size.yllcorner) / habitat_size.cellsize) - 1
            pts_remapped[:,2] = np.ceil((points[:,1] - habitat_size.xllcorner) / habitat_size.cellsize) - 1
        except IndexError:
            raise RuntimeError('Error extracting locations from text list file. Please check file format.')
        return pts_remapped


    @staticmethod
    def load_graph(filename):
        """Returns data for arbitrary graph or focal node list from file."""
        CSIO._check_file_exists(filename)
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
        CSIO._check_file_exists(filename)
        try:
            point_strengths = np.loadtxt(filename)
            if len(point_strengths.shape) == 1:
                point_strengths = np.array([point_strengths])
        except ValueError:
            raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
        
        i = np.argsort(point_strengths[:,0])
        return point_strengths[i]
    
    @staticmethod
    def deleterow(A, delrow):
        m = A.shape[0]
        n = A.shape[1]
        keeprows = np.delete(np.arange(0, m), delrow)
        keepcols = np.arange(0, n)
        return A[keeprows][:,keepcols]

#     def writeGraph(self,filename,graph,nodeNames):
#         """Save graph to disk in 3-column format."""  
#         graphNcol = self.convert_graph_to_3_col(graph,nodeNames)
#         savetxt(filename,graphNcol)
#         return

    @staticmethod            
    def _resample_map(reading_mask, resample_to, header, pmap):
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
                map_rc = CSIO.deleterow(map_rc, delrows2)
                
            rows = map_rc[:,1] 
            delrows = np.asarray(np.where(rows > resample_to.nrows - 1))
            delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2 != []:
                map_rc = CSIO.deleterow(map_rc, delrows2)
                
            cols = map_rc[:,2]
            delrows = np.asarray(np.where(cols < 0))
            delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2 != []:
                map_rc = CSIO.deleterow(map_rc, delrows2)
                
            cols = map_rc[:,2]
            delrows = np.asarray(np.where(cols > resample_to.ncols - 1))
            delrows2 = np.zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2 != []:
                map_rc = CSIO.deleterow(map_rc, delrows2)
            del delrows
            del delrows2
    
            #From map_rc to pmap:
            pmap = np.zeros((resample_to.nrows, resample_to.ncols), int)
            pmap[map_rc[:,1], map_rc[:,2]] = map_rc[:,0] 
            if reading_mask == True:
                pmap = 1-pmap #now zeros are areas to mask out, ones are kept
            
        except: 
            raise RuntimeError('Error resampling focal node, mask, or short-circuit region locations to match resistance map cell size and extent.  We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.')   
    
        return pmap

    @staticmethod
    def read_poly_map(filename, reading_mask, nodata_as, habitat_size, resample, file_type, data_type):
        """Reads raster maps for short-circuit regions (aka polygon poly_map), focal nodes, masks, current sources or grounds from disk."""  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata, _filetype) = header = CSIO._ascii_grid_read_header(filename)
        
        poly_map = CSIO._ascii_grid_reader(filename, data_type)
        if nodata_as != None:
            poly_map = np.where(poly_map==nodata, nodata_as, poly_map)
    
        if cellsize != habitat_size.cellsize:
            if resample:
                CSIO.logger.warning(CSIO.MSG_RESAMPLE % (file_type, "cell size",))
                poly_map = CSIO._resample_map(reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(CSIO.MSG_NO_RESAMPLE % (file_type, "cell_size"))            
        elif ncols != habitat_size.ncols:
            if resample:
                CSIO.logger.warning(CSIO.MSG_RESAMPLE % (file_type, "number of columns",))
                poly_map = CSIO._resample_map(reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(CSIO.MSG_NO_RESAMPLE % (file_type, "number of columns"))            
        elif nrows != habitat_size.nrows:
            if resample:
                CSIO.logger.warning(CSIO.MSG_RESAMPLE % (file_type, "number of rows",))
                poly_map = CSIO._resample_map(reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(CSIO.MSG_NO_RESAMPLE % (file_type, "number of rows"))            
        elif xllcorner != habitat_size.xllcorner:
            if resample:
                CSIO.logger.warning(CSIO.MSG_RESAMPLE % (file_type, "xllcorner",))
                poly_map = CSIO._resample_map(reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(CSIO.MSG_NO_RESAMPLE % (file_type, "xllcorner"))            
        elif yllcorner != habitat_size.yllcorner:
            if resample:
                CSIO.logger.warning(CSIO.MSG_RESAMPLE % (file_type, "yllcorner",))
                poly_map = CSIO._resample_map(reading_mask, habitat_size, header, poly_map)
            else:
                raise RuntimeError(CSIO.MSG_NO_RESAMPLE % (file_type, "yllcorner"))            
    
        if reading_mask==True:
            poly_map = np.where(poly_map < 0, 0, poly_map)        
    
        return poly_map

    @staticmethod
    def read_point_map(filename, file_type, habitat_size):
        """Reads map or text list of focal nodes from disk.
        
        File extension is used to determine whether format is ascii grid, numpy array, or text list.
        """
        filetype = CSIO._guess_file_type(filename)
        
        if filetype == CSIO.FILE_TYPE_TXTLIST:
            points_rc = CSIO._txt_list_reader(filename, 'int32', habitat_size)
        elif (filetype == CSIO.FILE_TYPE_AAGRID) or (filetype == CSIO.FILE_TYPE_NPY): # We use Numpy format for quickly passing grids between ArcGIS and Circuitscape.
            point_map = CSIO.read_poly_map(filename, False, 0, habitat_size, True, file_type, 'int32')
            (rows, cols) = np.where(point_map > 0)
    
            values = np.zeros(rows.shape, dtype='int32') 
            for i in range(0, rows.size):
                values[i] = point_map[rows[i], cols[i]]
            points_rc = np.c_[values, rows, cols]
        else:
            raise RuntimeError('Focal node file must be in one of text list, ascii grid or numpy array format.')
    
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
        #FIXME: reader does not currently handle infinite inputs for ground conductances.
        
        filetype = CSIO._guess_file_type(source_filename)

        if filetype == CSIO.FILE_TYPE_TXTLIST:  
            sources_rc = CSIO._txt_list_reader(source_filename, 'float64', habitat_size)
            source_map = np.zeros((habitat_size.nrows, habitat_size.ncols), dtype='float64')
            sources_rc_int = sources_rc.astype('int32')
            source_map[sources_rc_int[:,1], sources_rc_int[:,2]] = sources_rc[:,0]
        elif filetype == CSIO.FILE_TYPE_AAGRID:
            source_map = CSIO.read_poly_map(source_filename, False, 0, habitat_size, False, "Current source", 'float64')
            source_map = np.where(source_map==-9999, 0, source_map)
        else:
            raise RuntimeError('Current source files must either be in text list or ascii grid format.')
        
        if options.use_unit_currents == True:
            source_map = np.where(source_map, 1, source_map)
    
    
        filetype = CSIO._guess_file_type(ground_filename)
        
        if filetype == CSIO.FILE_TYPE_TXTLIST:
            grounds_rc = CSIO._txt_list_reader(ground_filename, 'float64', habitat_size) 
            ground_map_raw = -9999 * np.ones((habitat_size.nrows, habitat_size.ncols), dtype = 'float64')
            grounds_rc_int = grounds_rc.astype('int32')
            ground_map_raw[grounds_rc_int[:,1], grounds_rc_int[:,2]] = grounds_rc[:,0]
        elif filetype == CSIO.FILE_TYPE_AAGRID:
            ground_map_raw = CSIO.read_poly_map(ground_filename, False, None, habitat_size, False, "Ground", 'float64')
        else:
            raise RuntimeError('Ground files must either be in text list or ascii grid format.')
    
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
        
        """
        CSIO._check_file_exists(filename)
        filetype = CSIO._guess_file_type(filename)
        
        try:
            if filetype == CSIO.FILE_TYPE_INCL_PAIRS_AAGRID:
                with open(filename, 'r') as f:
                    minval = float(string.split(f.readline())[1])
                    maxval = float(string.split(f.readline())[1])
                    
                included_pairs = np.loadtxt(filename, skiprows=2, dtype='float')
                point_ids = included_pairs[1:,0]
                included_pairs = included_pairs[1:, 1:]
                included_pairs = np.where(included_pairs>maxval, 0, included_pairs)
                pair_list = np.where(included_pairs >= minval)
                pair_list = (point_ids[pair_list[0]], point_ids[pair_list[1]])
                pair_list = np.vstack(pair_list).T
                
                mode = "include"
            elif filetype == CSIO.FILE_TYPE_INCL_PAIRS:
                with open(filename, 'r') as f:
                    mode = string.split(f.readline())[1]
                pair_list = np.loadtxt(filename, skiprows=1, dtype='int32')
                point_ids = np.unique(pair_list)
                
            I = pair_list[:,0]
            J = pair_list[:,1]
            V = np.ones(pair_list.shape[0])
            max_node = np.max(pair_list)+1
            included_pairs = sparse.csr_matrix((V, (I, J)), shape=(max_node, max_node))
        except Exception as ex:
            CSIO.logger.exception(ex)
            raise RuntimeError('Error reading focal node include/exclude matrix. Please check file format.')
    
        return (mode, point_ids, included_pairs)


    @staticmethod
    @print_rusage
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
    
        file_type = CSIO._guess_file_type(options.habitat_file)        
        out_base = os.path.splitext(options.output_file)[0]
        out_file = out_base + '_' + grid_type + fileadd + ('.npy' if (file_type == CSIO.FILE_TYPE_NPY) else '.asc')
        CSIO._ascii_grid_writer(out_file, file_type, data, state, options.compress_grids)


    @staticmethod
    def read_cell_map(habitat_map, is_resistances, reclass_file, state):
        """Reads resistance or conductance raster into memory, converts former to conductance format."""  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata, _filetype) = CSIO._ascii_grid_read_header(habitat_map)
        state.ncols = ncols
        state.nrows = nrows
        state.xllcorner = xllcorner
        state.yllcorner = yllcorner
        state.cellsize = cellsize
        state.nodata = -9999 if (nodata == False) else nodata 
    
        cell_map = CSIO._ascii_grid_reader(habitat_map, 'float64')
    
        # Reclassification code
        if reclass_file != None:
            try:
                reclass_table = CSIO.read_point_strengths(reclass_file)
            except:
                raise RuntimeError('Error reading reclass table')
            for i in range (0, reclass_table.shape[0]):
                cell_map = np.where(cell_map==reclass_table[i,0], reclass_table[i,1], cell_map)
            CSIO.logger.info('Reclassified habitat map using %s'%(reclass_file,))
        
        if is_resistances == True:
            zeros_in_resistance_map = (np.where(cell_map==0, 1, 0)).sum() > 0
            if zeros_in_resistance_map == True: 
                raise RuntimeError('Error: zero resistance values are not currently supported for habitat maps.  Use a short-circuit region file instead.')
            g_map = 1 / cell_map  
            g_map = np.where(cell_map == -9999, 0, g_map)
        else:
            g_map = np.where(cell_map == -9999, 0, cell_map)    
        state.g_map = np.where(g_map < 0, 0, g_map)    


    @staticmethod
    def write_resistances_3columns(outfile_template, resistances_3columns):    
        """Writes effective resistances to disk in 3 column format."""  
        out_base, out_ext = os.path.splitext(outfile_template)
        out_file = out_base + '_resistances_3columns' + out_ext
        np.savetxt(out_file, resistances_3columns, fmt='%.10g')
        return         

    @staticmethod
    def write_resistances(outfile_template, resistances, resistances_3columns=None, incomplete=False):
        """Writes resistance file to disk."""
        out_base, out_extn = os.path.splitext(outfile_template)
        out_file = out_base + '_resistances' + ('_incomplete' if incomplete else '') + out_extn
        np.savetxt(out_file, resistances, fmt='%.10g')
        
        if resistances_3columns != None:
            CSIO.write_resistances_3columns(outfile_template, resistances_3columns)
        
        #remove partial result file        
        old_file = out_base + '_resistances_incomplete' + out_extn
        try:
            os.remove(old_file)
        except:
            pass 
    
    @staticmethod
    @print_rusage
    def write_currents(outfile_template, branch_currents, node_currents, fileadd, options):
        """Writes currents from network operations.
        
           Inputs are arrays with node names.
        """
        out_base, _out_ext = os.path.splitext(outfile_template)
        
        if fileadd != '':
            fileadd = ('_' + fileadd)
        elif options.scenario != 'advanced':
            fileadd = '_cum' #For backward compatibility
        if branch_currents != None:
            filename = out_base + '_branch_currents' + fileadd + '.txt'
            np.savetxt(filename, branch_currents, fmt='%.10g')
        filename = out_base + '_node_currents' + fileadd + '.txt'
        np.savetxt(filename, node_currents, fmt='%.10g') 

    @staticmethod
    @print_rusage
    def write_voltages(outfile_template, voltages, node_names, fileadd):
        """Saves voltage values from solving arbitrary graphs to disk."""
        output_voltages = np.zeros((len(voltages),2), dtype='float64')
        if node_names != None:
            output_voltages[:,0] = node_names[:]
        output_voltages[:,1] = voltages[:]

        out_base, _out_extn = os.path.splitext(outfile_template)
        if fileadd != '':
            fileadd = ('_' + fileadd)
        filename = out_base + '_voltages' + fileadd + '.txt'
        np.savetxt(filename, output_voltages, fmt='%.10g')      

    @staticmethod
    def match_headers(t_file, match_files):
        (ncols, nrows, xllcorner, yllcorner, cellsize, _nodata, _filetype) = CSIO._ascii_grid_read_header(t_file)
        
        for m_file in match_files:
            (ncols1, nrows1, xllcorner1, yllcorner1, cellsize1, _nodata1, _filetype) = CSIO._ascii_grid_read_header(m_file)
            if (ncols1 != ncols) or (nrows1 != nrows) or (abs(xllcorner1 - xllcorner) > cellsize/3) or (abs(yllcorner1 - yllcorner) > cellsize/3) or (cellsize1 != cellsize):
                return False
        return True

    @staticmethod
    def problem_size(data_type, habitat_file):
        # TODO: enhance to consider other parameters and be more accurate
        if data_type == 'network':
            graph_list = CSIO.load_graph(habitat_file)
            max_node = int(max(max(graph_list[:,0]), max(graph_list[:,2])))
            return max_node * max_node
        else:
            (ncols, nrows, _xllcorner, _yllcorner, _cellsize, _nodata, _filetype) = CSIO._ascii_grid_read_header(habitat_file)
            return ncols * nrows
        
            
        

