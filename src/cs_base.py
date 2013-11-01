##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import sys, time, gc, traceback, logging, inspect
import numpy as np
from scipy.sparse.linalg import cg
from scipy import sparse
from pyamg import smoothed_aggregation_solver
from scipy.sparse.csgraph import connected_components

from cs_cfg import CSConfig
from cs_state import CSState
from cs_io import CSIO

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
    
        
    def log_complete_job(self):
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
                    self.points_rc = CSBase.deleterow(self.points_rc, point)  
             
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


    def exists_points_in_component(self, comp, habitat):
        """Checks to see if there are focal points in a given component.
        
        Last parameter can either be a habitat object or a tuple of components and node_map.
        In network mode, components and node_map both are vectors.
        In raster mode, components is a vector and node_map is a matrix.
        """
        if type(habitat) == tuple:
            components, node_map = habitat
        else:
            components = habitat.components
            node_map = habitat.node_map
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
        
    def point_pair_idxs_in_component(self, comp, habitat):
        """Returns pairs of point indices that belong to the same component.
        
        Returns (x, -1) to denote end of each first node number.
        """
        if type(habitat) == tuple:
            components, node_map = habitat
        else:
            components = habitat.components
            node_map = habitat.node_map
        
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




class CSHabitatGraph:
    MODE_RASTER = 1
    MODE_NETWORK = 2
    
    def __init__(self, g_map=None, poly_map=None, connect_using_avg_resistances=False, connect_four_neighbors_only=False, g_graph=None, node_names=None):
        if None != g_map:
            self.mode = CSHabitatGraph.MODE_RASTER
            self.g_map = g_map
            self.poly_map = poly_map
            self.connect_using_avg_resistances = connect_using_avg_resistances
            self.connect_four_neighbors_only = connect_four_neighbors_only
            
            self.node_map = CSHabitatGraph._construct_node_map(g_map, poly_map)
            (component_map, components) = CSHabitatGraph._construct_component_map(g_map, self.node_map, connect_using_avg_resistances, connect_four_neighbors_only)
            self.component_map = component_map
            self.components = components
            
            self.num_components = components.max()
            self.num_nodes = self.node_map.max()
        else:
            self.mode = CSHabitatGraph.MODE_NETWORK
            self.g_graph = g_graph          # is the sparse CSR matrix 
            self.node_map = node_names    # list of node names
            
            (_num_components, C) = connected_components(g_graph)
            C += 1
            self.components = C
            
            self.num_components = C.max()
            self.num_nodes = self.node_map.size
    
    def prune_nodes_for_component(self, keep_component):
        """Removes nodes outside of component being operated on.
        
        Returns node map and adjacency matrix that only include nodes in keep_component.
        """
        if self.mode == CSHabitatGraph.MODE_RASTER:
            selector = self.component_map == keep_component
            
            g_map_pruned = selector * self.g_map
            poly_map_pruned = []
            if self.poly_map !=  []:
                poly_map_pruned = selector * self.poly_map
    
            node_map_pruned = CSHabitatGraph._construct_node_map(g_map_pruned, poly_map_pruned)
            G_pruned = CSHabitatGraph._construct_g_graph(g_map_pruned, node_map_pruned, self.connect_using_avg_resistances, self.connect_four_neighbors_only)
             
            return (G_pruned, node_map_pruned)
        else:
            del_indices = np.where(self.components != keep_component)
            pruned_graph = CSBase.deleterowcol(self.g_graph, delrow=del_indices, delcol=del_indices)
            indices = np.where(self.components == keep_component)
            nodes_in_component = self.node_map[indices]
            return (pruned_graph, nodes_in_component)                
        
    @staticmethod
    @print_timing
    def _construct_node_map(g_map, poly_map):
        """Creates a grid of node numbers corresponding to raster pixels with non-zero conductances."""  
        node_map = np.zeros(g_map.shape, dtype='int32')
        node_map[g_map.nonzero()] = np.arange(1, np.sum(g_map>0)+1, dtype='int32')

        if poly_map == []:
            return node_map

        # Remove extra points from poly_map that are not in g_map
        poly_map_pruned = np.zeros(g_map.shape, dtype='int32')
        poly_map_pruned[np.where(g_map)] = poly_map[np.where(g_map)]
        
        polynums = np.unique(poly_map)
   
        for i in range(0, polynums.size):
            polynum = polynums[i]
            if polynum !=  0:
                (pi, pj) = np.where(poly_map_pruned == polynum) #
                (pk, pl) = np.where(poly_map == polynum) #Added 040309 BHM                
                if len(pi) > 0:  
                    node_map[pk, pl] = node_map[pi[0], pj[0]] #Modified 040309 BHM  
        node_map[np.where(node_map)] = CSBase.relabel(node_map[np.where(node_map)], 1) #BHM 072409

        return node_map

    @staticmethod
    @print_timing
    def _construct_component_map(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only):
        """Assigns component numbers to grid corresponding to pixels with non-zero conductances.
        
        Nodes with the same component number are in single, connected components.
        
        """  
        G = CSHabitatGraph._construct_g_graph(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only) 
        (_num_components, C) = connected_components(G)
        C += 1 # Number components from 1

        (I, J) = np.where(node_map)
        nodes = node_map[I, J].flatten()

        component_map = np.zeros(node_map.shape, dtype = 'int32')
        component_map[I, J] = C[nodes-1]

        return (component_map, C)
    
    @staticmethod
    @print_timing
    def _construct_g_graph(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only):
        """Construct sparse adjacency matrix given raster maps of conductances and nodes."""
        numnodes = node_map.max()
        (node1, node2, conductances) = CSHabitatGraph._get_conductances(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only)
        return CSHabitatGraph._make_sparse_csr(node1, node2, conductances, numnodes)
        
    @staticmethod
    def _make_sparse_csr(node1, node2, conductances, numnodes):
        G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes)) # Memory hogging operation?
        g_graph = G + G.T
        
        return g_graph


    @staticmethod
    def _neighbors_horiz(g_map):
        """Returns values of horizontal neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        g_map_l = g_map[:, 0:(n-1)]
        g_map_r = g_map[:, 1:n]
        g_map_lr = np.double(np.logical_and(g_map_l, g_map_r))
        s_horiz = np.where(np.c_[g_map_lr, np.zeros((m,1), dtype='int32')].flatten())
        t_horiz = np.where(np.c_[np.zeros((m,1), dtype='int32'), g_map_lr].flatten())

        return (s_horiz, t_horiz)

        
    @staticmethod
    def _neighbors_vert(g_map):
        """Returns values of vertical neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        g_map_u = g_map[0:(m-1), :]
        g_map_d = g_map[1:m    , :]
        g_map_ud = np.double(np.logical_and(g_map_u, g_map_d))
        s_vert = np.where(np.r_[g_map_ud, np.zeros((1,n), dtype='int32')].flatten())
        t_vert = np.where(np.r_[np.zeros((1,n), dtype='int32'), g_map_ud].flatten())
        
        return (s_vert, t_vert)

        
    @staticmethod
    def _neighbors_diag1(g_map):
        """Returns values of 1st diagonal neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        z1 = np.zeros((m-1, 1), dtype='int32')
        z2 = np.zeros((1  , n), dtype='int32')
        
        g_map_ul  = g_map[0:m-1, 0:n-1]
        g_map_dr  = g_map[1:m  , 1:n  ]
        g_map_udr = np.double(np.logical_and(g_map_ul, g_map_dr)) 
        s_dr      = np.where(np.r_[np.c_[g_map_udr, z1], z2].flatten())
        t_dr      = np.where(np.r_[z2, np.c_[z1, g_map_udr]].flatten())
        
        return (s_dr, t_dr)

        
    @staticmethod
    def _neighbors_diag2(g_map):
        """Returns values of 2nd diagonal neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        z1 = np.zeros((m-1, 1), dtype='int32')
        z2 = np.zeros((1  , n), dtype='int32')

        g_map_ur  = g_map[0:m-1, 1:n  ]
        g_map_dl  = g_map[1:m  , 0:n-1]
        g_map_udl = np.double(np.logical_and(g_map_ur, g_map_dl)) 
        s_dl      = np.where(np.r_[np.c_[z1, g_map_udl], z2].flatten())
        t_dl      = np.where(np.r_[z2, np.c_[g_map_udl, z1]].flatten())
                        
        return (s_dl, t_dl)
        
    @staticmethod
    def _get_conductances(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only):
        """Calculates conductances between adjacent nodes given a raster conductance map.
        
        Returns an adjacency matrix with values representing node-to-node conductance values.
        
        """  
        (s_horiz, t_horiz) = CSHabitatGraph._neighbors_horiz(g_map)
        (s_vert,  t_vert)  = CSHabitatGraph._neighbors_vert(g_map)

        s = np.c_[s_horiz, s_vert].flatten()
        t = np.c_[t_horiz, t_vert].flatten()

        # Conductances
        g1 = g_map.flatten()[s]
        g2 = g_map.flatten()[t]

        if connect_using_avg_resistances == False:
            conductances = (g1+g2)/2
        else:
            conductances = 1 /((1/g1+1/g2)/2)

        if connect_four_neighbors_only == False:
            (s_dr, t_dr) = CSHabitatGraph._neighbors_diag1(g_map)
            (s_dl, t_dl) = CSHabitatGraph._neighbors_diag2(g_map)

            sd = np.c_[s_dr, s_dl].flatten()
            td = np.c_[t_dr, t_dl].flatten()

            # Conductances
            g1 = g_map.flatten()[sd]
            g2 = g_map.flatten()[td]

            if connect_using_avg_resistances == False:
                conductances_d = (g1+g2) / (2*np.sqrt(2))
            else:
                conductances_d =  1 / (np.sqrt(2)*(1/g1 + 1/g2) / 2)

            conductances = np.r_[conductances, conductances_d]

            s = np.r_[s, sd].flatten()
            t = np.r_[t, td].flatten()

        # Nodes in the g_graph. 
        # Subtract 1 for Python's 0-based indexing. Node numbers start from 1
        node1 = node_map.flatten()[s]-1
        node2 = node_map.flatten()[t]-1
        
        return (node1, node2, conductances)


class CSOutput:
    """Handles output of current and voltage maps"""
    def __init__(self, options, state, report_status, g_shape=None):
        self.options = options
        self.state = state
        self.report_status = report_status
        self.g_shape = g_shape
        
        self.is_network = (self.options.data_type == 'network')
        self.scenario = self.options.scenario
        
        if self.scenario == 'pairwise':
            if self.is_network:
                self._network_module_alloc_current_maps()
            else:
                self._pairwise_module_alloc_current_maps()

    
    def write_cum_current_map(self, node_map):
        if not self.options.write_cur_maps:
            return
        
        if self.is_network:
            self.full_branch_currents = CSOutput._convert_graph_to_3_col(self.full_branch_currents, node_map)
            self.full_node_currents =  CSOutput._append_names_to_node_currents(self.full_node_currents, node_map)
            
            ind = np.lexsort((self.full_branch_currents[:, 1], self.full_branch_currents[:, 0]))
            self.full_branch_currents = self.full_branch_currents[ind]                        
                
            ind = np.lexsort((self.full_node_currents[:, 1], self.full_node_currents[:, 0]))
            self.full_node_currents = self.full_node_currents[ind]                        

            CSIO.write_currents(self.options.output_file, self.full_branch_currents, self.full_node_currents, 'cum')
        else:
            CSIO.write_aaigrid('cum_curmap', '', self._log_transform(self.cum_current_map), self.options, self.state)
            
            if self.options.write_max_cur_maps:      
                CSIO.write_aaigrid('max_curmap', '', self._log_transform(self.max_current_map), self.options, self.state)


    def write_current_map(self, G, local_src, local_dst, voltages, node_map, frompoint, topoint, point_number, num_points):
        if not self.options.write_cur_maps:
            return

        if self.report_status==True:
            #self.log ('writing current map ' + str(point_number) + ' of ' + str(num_points) + '.',1)
            logging.info('writing current map ' + str(point_number) + ' of ' + str(num_points) + '.')
            
        finitegrounds = [-9999] #create dummy value for pairwise case
        
        if self.is_network:
            (node_currents, branch_currents) = self._create_current_maps(voltages, G, finitegrounds)
            
            if not self.options.write_cum_cur_map_only:                
                # Append node names and convert to array format
                branch_currents_array = CSOutput._convert_graph_to_3_col(branch_currents, node_map)
                node_currents_array = CSOutput._append_names_to_node_currents(node_currents, node_map)                
                CSIO.write_currents(self.options.output_file, branch_currents_array, node_currents_array, str(frompoint) + '_' + str(topoint))

            self.full_node_currents[node_map] += node_currents
            
            branch_currents_array = CSOutput._convert_graph_to_3_col(branch_currents, node_map)
            
            self.full_branch_currents = self.full_branch_currents + sparse.csr_matrix((branch_currents_array[:,2], (branch_currents_array[:,0], branch_currents_array[:,1])), shape=self.full_branch_currents.shape) 
        else:
            try:
                current_map = self._create_current_maps(voltages, G, finitegrounds, node_map)   
            except MemoryError:
                CSBase.enable_low_memory(False)
                current_map = self._create_current_maps(voltages, G, finitegrounds, node_map)                                                                                    
    
            if self.options.set_focal_node_currents_to_zero==True:
                # set source and target node currents to zero
                focal_node_pair_map = np.where(node_map == local_src+1, 0, 1)
                focal_node_pair_map = np.where(node_map == local_dst+1, 0, focal_node_pair_map)                                                
                current_map = np.multiply(focal_node_pair_map, current_map)
                del focal_node_pair_map
                
            self.cum_current_map = self.cum_current_map + current_map
             
            if self.options.write_max_cur_maps:
                self.max_current_map = np.maximum(self.max_current_map, current_map)
                 
            if not self.options.write_cum_cur_map_only:
                CSIO.write_aaigrid('curmap', '_' + str(frompoint) + '_' + str(topoint), self._log_transform(current_map), self.options, self.state)



    def write_voltage_map(self, voltages, node_map, frompoint, topoint, point_number, num_points):
        if not self.options.write_volt_maps:
            return
        
        if self.report_status==True:
            logging.info('writing voltage map ' + str(point_number) + ' of ' + str(num_points) + '.')

        if self.is_network:
            CSIO.write_voltages(self.options.output_file, voltages, node_map, str(frompoint) + '_' + str(topoint))
        else:
            voltage_map = self._create_voltage_map(node_map, voltages) 
            CSIO.write_aaigrid('voltmap', '_' + str(frompoint) + '_' + str(topoint), voltage_map, self.options, self.state)

    def _log_transform(self, map_to_transform):
        if self.options.log_transform_maps:
            map_to_transform = np.where(map_to_transform > 0, np.log10(map_to_transform), self.state.nodata)
        return map_to_transform

    def _network_module_alloc_current_maps(self):
        self.full_branch_currents = sparse.csr_matrix(self.g_shape)
        self.full_node_currents = np.zeros((self.g_shape[0], 1), dtype='float64')

        
    def _pairwise_module_alloc_current_maps(self):
        self.cum_current_map = self.max_current_map = []
        if self.options.write_cur_maps == True:
            self.cum_current_map = np.zeros((self.state.nrows, self.state.ncols), dtype='float64') 
            if self.options.write_max_cur_maps == True:
                self.max_current_map = self.cum_current_map

    @print_timing
    def _create_current_maps(self, voltages, G, finitegrounds, node_map=None):
        """In raster mode, returns raster current map given node voltage vector, adjacency matrix, etc.
        In network mode returns node and branch currents given voltages in arbitrary graphs.
        """  
        gc.collect()
        G =  G.tocoo()
        node_currents = CSOutput._get_node_currents(voltages, G, finitegrounds)
        
        if self.is_network:
            node_currents_col = np.zeros((node_currents.shape[0],1), dtype='float64')
            node_currents_col[:,0] = node_currents[:]
            branch_currents = CSOutput._get_branch_currents(G, voltages, True) 
            branch_currents = np.absolute(branch_currents) 
            return node_currents_col, branch_currents
        else:
            (rows, cols) = np.where(node_map)
            vals = node_map[rows, cols]-1
            current_map = np.zeros((self.state.nrows, self.state.ncols), dtype='float64')
            current_map[rows,cols] = node_currents[vals]    
            return current_map


    @staticmethod
    def _get_node_currents(voltages, G, finitegrounds):
        """Calculates currents at nodes."""  
        node_currents_pos = CSOutput._get_node_currents_posneg(G, voltages, finitegrounds, True) 
        node_currents_neg = CSOutput._get_node_currents_posneg(G, voltages, finitegrounds, False)
        node_currents = np.where(node_currents_neg > node_currents_pos, node_currents_neg, node_currents_pos)

        return np.asarray(node_currents)[0]

    
    @staticmethod
    def _get_node_currents_posneg(G, voltages, finitegrounds, pos):
        """Calculates positive or negative node currents based on pos flag."""  
        branch_currents = CSOutput._get_branch_currents(G, voltages, pos)
        branch_currents = branch_currents - branch_currents.T #Can cause memory error
        
        branch_currents = branch_currents.tocoo() #Can cause memory error, but this and code below more memory efficient than previous version.
        mask = branch_currents.data > 0
        row  = branch_currents.row[mask]
        col  = branch_currents.col[mask]
        data = branch_currents.data[mask]
        del mask
        n = G.shape[0]
        branch_currents = sparse.csr_matrix((data, (row, col)), shape = (n,n))
           
        if finitegrounds[0]!= -9999:  
            finiteground_currents = np.multiply(finitegrounds, voltages)
            if pos:
                finiteground_currents = np.where(finiteground_currents < 0, -finiteground_currents, 0)
            else:
                finiteground_currents = np.where(finiteground_currents > 0, finiteground_currents, 0)  
            n = G.shape[0]
            branch_currents = branch_currents + sparse.spdiags(finiteground_currents.T, 0, n, n)        

        return branch_currents.sum(0)
    
    
    @staticmethod
    def _get_branch_currents(G, voltages, pos):    
        """Calculates branch currents."""  
        branch_currents = CSOutput._get_branch_currents_posneg(G, voltages, pos)
        n = G.shape[0]
        mask = G.row < G.col
        branch_currents = sparse.csr_matrix((branch_currents, (G.row[mask], G.col[mask])), shape = (n,n)) #SQUARE MATRIX, SAME DIMENSIONS AS GRAPH
        return branch_currents


    @staticmethod
    def _get_branch_currents_posneg(G, voltages, pos):
        """Calculates positive or negative node currents based on pos flag."""  
        mask = G.row < G.col
        if pos:
            vdiff = voltages[G.row[mask]]              
            vdiff -=  voltages[G.col[mask]]             
        else:
            vdiff = voltages[G.col[mask]]              
            vdiff -=  voltages[G.row[mask]]             

        conductances = np.where(G.data[mask] < 0, -G.data[mask], 0)
        del mask
        
        branch_currents = np.asarray(np.multiply(conductances, vdiff.T)).flatten()
        maxcur = max(branch_currents)
        branch_currents = np.where(np.absolute(branch_currents/maxcur) < 1e-8, 0, branch_currents) #Delete very small branch currents to save memory
        return branch_currents

    @print_timing
    def _create_voltage_map(self, node_map, voltages):
        """Creates raster map of voltages given node voltage vector."""
        voltage_map = np.zeros((self.state.nrows, self.state.ncols), dtype = 'float64')
        ind = node_map > 0
        voltage_map[np.where(ind)] = np.asarray(voltages[node_map[ind]-1]).flatten()
        return voltage_map

    @staticmethod
    def _convert_graph_to_3_col(graph, node_names): 
        """Converts a sparse adjacency matrix to 3-column format."""  
        Gcoo =  graph.tocoo()
        mask = Gcoo.data > 0
        
        graph_n_col = np.zeros((Gcoo.row[mask].size, 3), dtype="float64") #Fixme: this may result in zero-current nodes being left out.  Needed to make change so dimensions would match Gcoo.data[mask]
        
        if node_names == None:
            graph_n_col[:,0] = Gcoo.row[mask]
            graph_n_col[:,1] = Gcoo.col[mask]
        else:
            graph_n_col[:,0] = node_names[Gcoo.row[mask]]
            graph_n_col[:,1] = node_names[Gcoo.col[mask]]
        graph_n_col[:,2] = Gcoo.data[mask]
        return graph_n_col


    @staticmethod
    def _append_names_to_node_currents(node_currents, node_names):
        """Adds names of focal nodes to node current lists."""    
        output_node_currents = np.zeros((len(node_currents),2), dtype='float64')
        output_node_currents[:,0] = node_names[:]
        try:
            output_node_currents[:,1] = node_currents[:,0]
        except:
            output_node_currents[:,1] = node_currents[:]
        return output_node_currents


