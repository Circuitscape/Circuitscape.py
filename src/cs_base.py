##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import sys, time, gc, traceback, logging, inspect
import numpy as np
from scipy.sparse.linalg import cg
from scipy import sparse
from pyamg import smoothed_aggregation_solver

from cs_cfg import CSConfig
from cs_state import CSState

#from numpy import *

print_timings_spaces = 0
print_timings = False

def print_timing_enabled(is_enabled):
    """Enables or disables the print_timings decorator."""
    global print_timings
    print_timings = is_enabled

def print_timing(func):
    """Prints time elapsed for functions with print_timings decorator."""  
    def wrapper(*arg):
        if not print_timings:
            return func(*arg)
        global print_timings_spaces
        print_timings_spaces +=  2
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print_timings_spaces -=  2
        print("%10d sec: %s%s"%((t2-t1), " "*print_timings_spaces, func.func_name))
        sys.stdout.flush()
        return res
    return wrapper

wx_available = True
try:
    import wx
except ImportError:
    wx_available = False
    wx = None

class CSBase(object):    
    """Circuitscape base class, common across all circuitscape modules"""
    def __init__(self, configFile, logger_func):
        gc.enable()
        np.seterr(invalid='ignore')
        np.seterr(divide='ignore')
        
        self.state = CSState()
        self.state.amg_hierarchy = None
        self.options = CSConfig(configFile)
        
        self.options.use_reclass_table = False
        self.options.reclass_file = './reclass.txt'        

        print_timing_enabled(self.options.print_timings)
        #print_timing_enabled(True)
        
        if logger_func == 'Screen':
            self.options.screenprint_log = True
            logger_func = None
        else:
            self.options.screenprint_log = False
        #self.options.screenprint_log = True
        
        self.logger = logger_func
        logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                            datefmt='%m/%d/%Y %I.%M.%S.%p',
                            level=logging.DEBUG)


    # TODO: should ultimately be replaced with the logging module.
    # logging to UI should be handled with another logger or handler
    def log(self, text, col):
        """Prints updates to GUI or python window."""
        if None == self.state.last_gui_yield_time:
            self.state.last_gui_yield_time = time.time()
        else: # Force update every 10 secs
            (hours,mins,secs) = self.elapsed_time(self.state.last_gui_yield_time)
            if (secs > 10) or (mins > 0) or (hours > 0):
                self.state.last_gui_yield_time = time.time()
                try:
                    if wx_available: 
                        wx.SafeYield(None, True)  
                        wx.GetApp().Yield(True)
                except:
                    pass
            
        if self.logger or (self.options.screenprint_log == True and len(text) > 1):
            text = '%s%s'%(' '*len(inspect.stack()), str(text))
        
            if self.logger:
                self.logger(text, col)
                
            if self.options.screenprint_log == True and len(text) > 1:
                if col == 1:
                    logging.info(text)
                else:
                    logging.debug(text)
                sys.stdout.flush()
    
        
    def logCompleteJob(self):
        """Writes total time elapsed at end of run."""
        (hours,mins,secs) = self.elapsed_time(self.state.start_time)
        if hours>0:
            self.log('Job took ' + str(hours) +' hours ' + str(mins) + ' minutes to complete.',2)
        else:
            self.log('Job took ' + str(mins) +' minutes ' + str(secs) + ' seconds to complete.',2)

        
    def enable_low_memory(self, restart):
        """Runs circuitscape in low memory mode.  Not incredibly helpful it seems."""  
        self.state.amg_hierarchy = None
        gc.collect()
        if self.options.low_memory_mode==True:
            if restart==False: #If this module has already been called
                raise MemoryError
        self.options.low_memory_mode = True
        print'\n**************\nMemory error reported.'        

        ex_type, value, tb = sys.exc_info()
        info = traceback.extract_tb(tb)
        print'Full traceback:'
        print info
        print'***************'
        filename, lineno, function, _text = info[-1] # last line only
        print"\n %s:%d: %s: %s (in %s)" %\
              (filename, lineno, ex_type.__name__, str(value), function)

        ex_type = value = tb = None # clean up
        print'\nWARNING: CIRCUITSCAPE IS RUNNING LOW ON MEMORY.'
        if restart==True:
            print'Restarting in low memory mode, which will take somewhat longer to complete.'
        else:
            print'Switching to low memory mode, which will take somewhat longer to complete.'            
        print'CLOSING OTHER PROGRAMS CAN HELP FREE MEMORY RESOURCES.'
        print'Please see the user guide for more information on memory requirements.\n'               
        if restart==True:
            print'***Restarting in low memory mode***\n'
        else:
            print'***Continuing in low memory mode***\n'
        return

        
    @print_timing
    def solve_linear_system(self, G, rhs): 
        """Solves system of equations."""  
        gc.collect()
        # Solve G*x = rhs
        x = []
        if self.options.solver == 'cg+amg':
            ml = self.state.amg_hierarchy
            G.psolve = ml.psolve
            (x, flag) = cg(G, rhs, tol = 1e-6, maxiter = 100000)
            if flag !=  0 or np.linalg.norm(G*x-rhs) > 1e-3:
                raise RuntimeError('CG did not converge. May need more iterations.') 

        if self.options.solver == 'amg':
            ml = self.state.amg_hierarchy
            x = ml.solve(rhs, tol = 1e-6);

        return x 
 
         
    @staticmethod
    @print_timing
    def laplacian(G): 
        """Returns Laplacian of graph."""  
        n = G.shape[0]

        # FIXME: Potential for memory savings, if assignment is used
        G = G - sparse.spdiags(G.diagonal(), 0, n, n)
        G = -G + sparse.spdiags(G.sum(0), 0, n, n)

        return G

        
    @print_timing
    def create_amg_hierarchy(self, G): 
        """Creates AMG hierarchy."""  
        if self.options.solver == 'amg' or self.options.solver == 'cg+amg':
            self.state.amg_hierarchy = None
            # construct the MG hierarchy
            ml = []
            #  scipy.io.savemat('c:\\temp\\graph.mat',mdict={'d':G})
            ml = smoothed_aggregation_solver(G)
            self.state.amg_hierarchy = ml


    @staticmethod
    def grid_to_graph (x, y, node_map):
        """Returns node corresponding to x-y coordinates in input grid."""  
        return node_map[x, y] - 1
        
    
    @staticmethod
    def elapsed_time(startTime): 
        """Returns elapsed time given a start time."""    
        now = time.time()
        elapsed = now-startTime
        secs = int(elapsed)
        mins = int(elapsed/60)
        hours = int(mins/60)
        mins = mins - hours*60
        secs = secs - mins*60 - hours*3600
        return hours, mins, secs
    
    @staticmethod
    def deleterow(A, delrow):
        m = A.shape[0]
        n = A.shape[1]
        keeprows = np.delete(np.arange(0, m), delrow)
        keepcols = np.arange(0, n)
        return A[keeprows][:,keepcols]
            
    @staticmethod
    def deletecol(A, delcol):
        m = A.shape[0]
        n = A.shape[1]
        keeprows = np.arange(0, m)
        keepcols = np.delete(np.arange(0, n), delcol)
        return A[keeprows][:,keepcols]
    
    @staticmethod
    def deleterowcol(A, delrow, delcol):
        m = A.shape[0]
        n = A.shape[1]
    
        keeprows = np.delete(np.arange(0, m), delrow)
        keepcols = np.delete(np.arange(0, n), delcol)
    
        return A[keeprows][:,keepcols]
    
    @staticmethod
    def relabel(oldlabel, offset=0):
        newlabel = np.zeros(np.size(oldlabel), dtype='int32')
        s = np.sort(oldlabel)
        perm = np.argsort(oldlabel)
        f = np.where(np.diff(np.concatenate(([s[0]-1], s))))
        newlabel[f] = 1
        newlabel = np.cumsum(newlabel)
        newlabel[perm] = np.copy(newlabel)
        return newlabel-1+offset





class CSFocalPoints:
    """Represents a set of focal points and associate logic to work with them"""
    MODE_RASTER = 1
    MODE_NETWORK = 2
    def __init__(self, points, included_pairs):
        self.included_pairs = included_pairs
        pdims = len(points.shape)
        if pdims == 2:
            # raster mode
            self.mode = CSFocalPoints.MODE_RASTER
            self.points_rc = points
            self._prune_included_pairs()
            self.point_ids = self._get_point_ids()
        elif pdims == 1:
            # network mode
            self.mode = CSFocalPoints.MODE_NETWORK
            self.points_rc = None
            self.point_ids = points
            self._prune_included_pairs()
        else:
            raise RuntimeError('Focal points should either be coordinates or node names. Got an array with %d dimensions'%(pdims,))

    
    def _prune_included_pairs(self):
        """Remove excluded points from focal node list when using extra file that lists pairs to include/exclude."""
        
        if self.included_pairs == None:
            return
        
        include_list = list(self.included_pairs[0,:])
        point = 0
        _drop_flag = False
        
        if self.mode == CSFocalPoints.MODE_NETWORK:
            while point < self.point_ids.size: #Prune out any points not in includeList
                if self.point_ids[point] in include_list: #match
                    point = point+1
                else:
                    _drop_flag = True   
                    self.point_ids = np.delete(self.point_ids, point)
             
            include_list = list(self.point_ids[:])            
        else:
            while point < self.points_rc.shape[0]: #Prune out any points not in include_list
                if self.points_rc[point,0] in include_list: #match
                    point = point+1
                else:
                    _drop_flag = True   
                    self.points_rc = self.deleterow(self.points_rc, point)  
             
            include_list = list(self.points_rc[:,0])
        
        num_connection_rows = self.included_pairs.shape[0]
        row = 1
        while row < num_connection_rows: #Prune out any entries in include_list that are not in points_rc
            if self.included_pairs[row,0] in include_list: #match
                row = row+1
            else:
                self.included_pairs = CSBase.deleterowcol(self.included_pairs, delrow=row, delcol=row)   
                _drop_flag = True
                num_connection_rows = num_connection_rows-1

        if _drop_flag==True:
            logging.info('\nNOTE: Code to exclude pairwise calculations is activated and some entries did not match with focal node file. Some focal nodes may have been dropped.')      


    def _get_point_ids(self):
        """Return a list of unique focal node IDs"""
        return np.unique(np.asarray(self.points_rc[:,0]))

    def num_points(self):
        if self.mode == CSFocalPoints.MODE_NETWORK:
            return self.point_ids.size
        else:
            return self.points_rc.shape[0]
    
    def get_unique_coordinates(self):
        """Return a list of unique focal node IDs and x-y coordinates."""
        if self.mode == CSFocalPoints.MODE_NETWORK:
            raise RuntimeError('Not available in network mode')
          
        points_rc_unique = np.zeros((self.point_ids.size, 3), int)
        for i in range(0, self.point_ids.size):
            for j in range(0, self.points_rc.shape[0]):
                if self.points_rc[j,0] == self.point_ids[i]:
                    points_rc_unique[i,:] = self.points_rc[j,:] 
                    break                    
        return points_rc_unique          


    def get_coordinates(self, pt_idx=None):
        """Returns a list of focal node IDs and x-y coordinates or only x-y coordinates if an index or ID is specified"""
        if self.mode == CSFocalPoints.MODE_NETWORK:
            raise RuntimeError('Not available in network mode')
          
        if pt_idx != None:
            return (self.points_rc[pt_idx,1], self.points_rc[pt_idx,2])
            
        return self.points_rc

    def get_subset(self, idx_list):
        """Returns a subset of focal point coordinates for supplied indexes"""
        if self.mode == CSFocalPoints.MODE_NETWORK:
            raise RuntimeError('Not available in network mode')
        
        ncoords = len(idx_list)
        sub_coords = np.zeros((ncoords,3), int)
        for idx in range(0, ncoords):
            sub_coords[idx,:] = self.points_rc[idx_list[idx], :]
        
        return CSFocalPoints(sub_coords, self.included_pairs)

    @staticmethod
    def grid_to_graph(x, y, node_map):
        """Returns node corresponding to x-y coordinates in input grid."""  
        return node_map[x, y] - 1


    def exists_points_in_component(self, comp, components, node_map):
        """Checks to see if there are focal points in a given component.
        
        In network mode, components and node_map both are vectors.
        In raster mode, components and node_map are matrices.
        """
        if self.mode == CSFocalPoints.MODE_NETWORK:
            indices = np.where(components == comp)
            nodes_in_component = node_map[indices]            
            include_list = nodes_in_component.tolist()
            for pt_id in self.point_ids:
                if pt_id in include_list:
                    return True
        else:
            numpoints = self.points_rc.shape[0]
            for pt1 in range(0, numpoints): 
                src = self.grid_to_graph(self.points_rc[pt1,1], self.points_rc[pt1,2], node_map)
                for pt2 in range(pt1+1, numpoints):
                    dst = self.grid_to_graph(self.points_rc[pt2,1], self.points_rc[pt2,2], node_map)
                    if (src >=  0 and components[src] == comp) and (dst >=  0 and components[dst] == comp):
                        return True
        return False
    
    
    def get_graph_node_idx(self, focal_point_idx, node_map):
        """Returns the index of the focal point node in the graph.
        
        In network mode, node_map is a vector. It can be a pruned node map representing only one component of the map.
        In raster mode, node_map is a matrix.
        """
        if self.mode == CSFocalPoints.MODE_NETWORK:
            return node_map.tolist().index(self.point_ids[focal_point_idx])
        else:
            return self.grid_to_graph(self.points_rc[focal_point_idx,1], self.points_rc[focal_point_idx,2], node_map)
    
    
    def point_id(self, idx):
        if self.mode == CSFocalPoints.MODE_NETWORK:
            return self.point_ids[idx]
        else:
            return self.points_rc[idx,0]    
    
    def point_pair_idxs(self):
        """Returns pairs of point indices across all components
        
        Returns (x, -1) to denote end of each first node number.
        """
        if self.mode == CSFocalPoints.MODE_NETWORK:
            numpoints = self.point_ids.size
            for n1_idx in range(0, numpoints-1):
                for n2_idx in range(n1_idx+1, numpoints):
                    yield (n1_idx, n2_idx)
                yield(n1_idx, -1)
        else:        
            numpoints = self.points_rc.shape[0]
            for pt1_idx in range(0, numpoints): 
                for pt2_idx in range(pt1_idx+1, numpoints):
                    # tan: is this check required? 
                    # after pruning focal points for included pairs, this should always point to a valid pair
                    if (None != self.included_pairs) and (self.included_pairs[pt1_idx+1, pt2_idx+1] != 1):
                        continue
                    yield (pt1_idx, pt2_idx)
                yield(pt1_idx, -1)
        
    def point_pair_idxs_in_component(self, comp, components, node_map):
        """Returns pairs of point indices that belong to the same component.
        
        Returns (x, -1) to denote end of each first node number.
        """
        if self.mode == CSFocalPoints.MODE_NETWORK:
            indices = np.where(components == comp)
            nodes_in_component = node_map[indices]            
            include_list = list(nodes_in_component[:])
            
            numpoints = self.point_ids.size
            for n1_idx in range(0, numpoints-1):
                if self.point_ids[n1_idx] not in include_list:
                    continue
                for n2_idx in range(n1_idx+1, numpoints):
                    if self.point_ids[n2_idx] in include_list:
                        yield (n1_idx, n2_idx)
                yield(n1_idx, -1)
        else:        
            numpoints = self.points_rc.shape[0]
            for pt1_idx in range(0, numpoints): 
                dst = self.get_graph_node_idx(pt1_idx, node_map)
                if (dst <  0 or components[dst] != comp):
                    continue
                for pt2_idx in range(pt1_idx+1, numpoints):
                    # tan: is this check required? 
                    # after pruning focal points for included pairs, this should always point to a valid pair
                    if (None != self.included_pairs) and (self.included_pairs[pt1_idx+1, pt2_idx+1] != 1):
                        continue
                    src = self.get_graph_node_idx(pt2_idx, node_map)
                    if (src >=  0 and components[src] == comp):
                        yield (pt1_idx, pt2_idx)
                yield(pt1_idx, -1)
                