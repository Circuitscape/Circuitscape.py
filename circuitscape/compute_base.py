import sys, time, logging, copy
import numpy as np
from scipy.sparse.linalg import cg
from scipy import sparse
from scipy.sparse.csgraph import connected_components, _validation #_validation is for py2exe

from cfg import CSConfig
from state import CSState
from io import CSIO
from profiler import ResourceLogger, print_rusage, gc_before, GCPreempt, LowMemRetry


class ComputeBase(object):
    logger = None
    
    """Circuitscape base class, common across all circuitscape modules"""
    def __init__(self, configFile, ext_log_handler):
        #gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_UNCOLLECTABLE | gc.DEBUG_SAVEALL)
        np.seterr(invalid='ignore')
        np.seterr(divide='ignore')

        self.state = CSState()
        self.options = CSConfig(configFile)
        self._setup_loggers(ext_log_handler)

        GCPreempt.enabled = self.options.preemptive_memory_release 
        LowMemRetry.callback = self._on_low_memory()
        LowMemRetry.max_retry = 0 if self.options.low_memory_mode else 1

        if self.options.parallelize: 
            if sys.platform.startswith('win'):
                self.options.parallelize = False
                ComputeBase.logger.warn("No support for parallelization on Windows. Option disabled.")


    @staticmethod
    def _create_formatter(log_lvl, printDashes=False):      
        if log_lvl == logging.DEBUG:
            formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s', '%m/%d/%Y %I.%M.%S.%p')
        elif printDashes==True:
            formatter = logging.Formatter('    --- %(message)s ---')
        else:
            formatter = logging.Formatter('%(message)s')
        return formatter
        

    @staticmethod
    def _create_logger(logger_name, log_lvl, log_file, screenprint_log, ext_log_handler, formatter=None):
        logger = logging.getLogger(logger_name)
        logger.setLevel(log_lvl)
        
        if formatter == None:
            formatter = ComputeBase._create_formatter(log_lvl)

        handlers = []
        
        if ext_log_handler:
            handlers.append(ext_log_handler)
            
        if log_file:
            handlers.append(logging.FileHandler(log_file))

        if screenprint_log:
            handlers.append(logging.StreamHandler())
        
        if len(handlers) == 0:  # if no loggers configured, disable logging
            handlers.append(logging.NullHandler())
        
        for handler in handlers:
            handler.setLevel(log_lvl)
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            
        return logger

    @staticmethod
    def _close_handlers(logger):
        all_handlers = list(logger.handlers)
        for handler in all_handlers:
            logger.removeHandler(handler)
            handler.flush()
            handler.close()
            
    def _setup_loggers(self, ext_log_handler):
        if None != ComputeBase.logger:  # logger has already been setup
            ComputeBase._close_handlers(ComputeBase.logger)
            ComputeBase._close_handlers(ResourceLogger.rlogger)
        
        # For backward compatibility with ArcGIS-based tools. 
        printDashes = True if (ext_log_handler == 'Screen' and self.options.version == 'unknown') else False
        if ext_log_handler == 'Screen':
            self.options.screenprint_log = True
            ext_log_handler = None

        log_lvl = getattr(logging, self.options.log_level.upper())
        
        formatter = ComputeBase._create_formatter(log_lvl, printDashes)
        
        CSIO.logger = CSState.logger = ComputeBase.logger = ComputeBase._create_logger('circuitscape', log_lvl, self.options.log_file, self.options.screenprint_log, ext_log_handler, formatter)

        if (self.options.profiler_log_file != self.options.log_file) and (self.options.profiler_log_file is not None):
            res_logger = ComputeBase._create_logger('circuitscape_profile', logging.DEBUG, self.options.profiler_log_file, self.options.screenprint_log, ext_log_handler, formatter)
        else:
            res_logger = ComputeBase.logger

        ResourceLogger.init_rusage(self.options.print_timings, self.options.print_rusages, res_logger)

        
    def log_complete_job(self):
        """Writes total time elapsed at end of run."""
        (hours,mins,secs) = self.elapsed_time(self.state.start_time)
        if hours>0:
            ComputeBase.logger.debug('Job took ' + str(hours) +' hours ' + str(mins) + ' minutes to complete.')
        else:
            ComputeBase.logger.info('Job took ' + str(mins) +' minutes ' + str(secs) + ' seconds to complete.')

    def _on_low_memory(self):
        """Runs circuitscape in low memory mode. Not incredibly helpful it seems."""
        cs = self
        @staticmethod
        def _set_low_memory_mode(etype, ex, traceback):
            cs.state.del_amg_hierarchy()
            ComputeBase.logger.exception("Circuitscape is running low on memory. Closing other programs can help free memory resources.")
            if not cs.options.low_memory_mode:
                ComputeBase.logger.warning("Switching to low memory mode, which will take somewhat longer to complete.")
                cs.options.low_memory_mode = True
                cs.options.preemptive_memory_release = True
        return _set_low_memory_mode

    
    @staticmethod
    @gc_before
    @print_rusage
    def solve_linear_system(G, rhs, solver_type, ml):
        """Solves system of equations."""
        # Solve G*x = rhs
        x = []
        if solver_type == 'cg+amg':
            G.psolve = ml.psolve
            (x, flag) = cg(G, rhs, tol = 1e-6, maxiter = 100000)
            if flag !=  0 or np.linalg.norm(G*x-rhs) > 1e-3:
                raise RuntimeError('CG did not converge. May need more iterations.') 
        elif solver_type == 'amg':
            x = ml.solve(rhs, tol = 1e-6);

        return x 
 
         
    @staticmethod
    @print_rusage
    def laplacian(G): 
        """Returns Laplacian of graph."""  
        n = G.shape[0]

        # FIXME: Potential for memory savings, if assignment is used
        G = G - sparse.spdiags(G.diagonal(), 0, n, n)
        G = -G + sparse.spdiags(G.sum(0), 0, n, n)

        return G


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


class IncludeExcludePairs:
    """Represents a set of focal points that are to be included/excluded during computation"""
    def __init__(self, filename):
        mode, point_ids, mat = CSIO.read_included_pairs(filename)
        self.is_include = (mode == "include")
        self.point_ids = point_ids
        self.mat = mat.tolil()
        self.max_id = int(np.max(point_ids)) + 1

    def is_included_pair(self, point_id1, point_id2):
        return (point_id1 < self.max_id) and (point_id2 < self.max_id) and (((self.mat[point_id1, point_id2] + self.mat[point_id2, point_id1]) > 0) == self.is_include)

    def is_included(self, point_id):
        return (point_id < self.max_id) and (((self.mat[point_id,:].sum() + self.mat[:,point_id].sum()) > 0) == self.is_include)
    
    def is_present(self, point_id):
        return (point_id in self.point_ids)
    
    def delete_point(self, point_id):
        if(point_id < self.max_id):
            self.mat[point_id, :] = 0
            self.mat[:, point_id] = 0

    def keep_only_points(self, points):
        for point_id in range(0, self.max_id):
            if point_id not in points:
                self.delete_point(point_id)

class FocalPoints:
    """Represents a set of focal points and associate logic to work with them"""
    def __init__(self, points, included_pairs, is_network):
        self.incl_pairs = copy.deepcopy(included_pairs)
        self.is_network = is_network
        if is_network:
            # network mode
            self.is_network = True
            self.points_rc = None
            self.point_ids = points
            self._prune_included_pairs()
        else:
            # raster mode
            self.is_network = False
            self.points_rc = points
            self._prune_included_pairs()
            self.point_ids = self._get_point_ids()            

    def _prune_included_pairs(self):
        """Remove excluded points from focal node list when using extra file that lists pairs to include/exclude."""
        if self.incl_pairs == None:
            return

        idx = 0
        _drop_flag = False
        
        if self.is_network:
            while idx < self.point_ids.size: #Prune out any points not in includeList
                if self.incl_pairs.is_present(self.point_ids[idx]):
                    idx += 1
                else:
                    _drop_flag = True   
                    self.point_ids = np.delete(self.point_ids, idx)
            self.incl_pairs.keep_only_points(self.point_ids)
        else:
            while idx < self.points_rc.shape[0]:
                if self.incl_pairs.is_present(self.points_rc[idx, 0]):
                    idx += 1
                else:
                    _drop_flag = True
                    self.points_rc = ComputeBase.deleterow(self.points_rc, idx)
            self.incl_pairs.keep_only_points(self.points_rc[:,0])
            
        if _drop_flag==True:
            ComputeBase.logger.info('Code to exclude pairwise calculations is activated and some entries did not match with focal node file. Some focal nodes may have been dropped.')      
            

    def _get_point_ids(self):
        """Return a list of unique focal node IDs"""
        return np.unique(np.asarray(self.points_rc[:,0]))

    def num_points(self):
        if self.is_network:
            return self.point_ids.size
        else:
            return self.points_rc.shape[0]
    
    def get_unique_coordinates(self):
        """Return a list of unique focal node IDs and x-y coordinates."""
        if self.is_network:
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
        if self.is_network:
            raise RuntimeError('Not available in network mode')
          
        if pt_idx != None:
            return (self.points_rc[pt_idx,1], self.points_rc[pt_idx,2])
            
        return self.points_rc

    def get_subset(self, idx_list):
        """Returns a subset of focal point coordinates for supplied indexes"""
        if self.is_network:
            raise RuntimeError('Not available in network mode')
        
        ncoords = len(idx_list)
        sub_coords = np.zeros((ncoords,3), int)
        for idx in range(0, ncoords):
            sub_coords[idx,:] = self.points_rc[idx_list[idx], :]
        
        return FocalPoints(sub_coords, self.incl_pairs, self.is_network)

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
        if self.is_network:
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
        if self.is_network:
            return node_map.tolist().index(self.point_ids[focal_point_idx])
        else:
            return self.grid_to_graph(self.points_rc[focal_point_idx,1], self.points_rc[focal_point_idx,2], node_map)
    
    
    def point_id(self, idx):
        if self.is_network:
            return self.point_ids[idx]
        else:
            return self.points_rc[idx,0]    
    
    def point_pair_idxs(self):
        """Returns pairs of point indices across all components
        
        Returns (x, -1) to denote end of each first node number.
        """
        if self.is_network:
            numpoints = self.point_ids.size
            for n1_idx in range(0, numpoints-1):
                for n2_idx in range(n1_idx+1, numpoints):
                    yield (n1_idx, n2_idx)
                yield(n1_idx, -1)
        else:        
            numpoints = self.point_ids.size
            
            for pt1_idx in range(0, numpoints): 
                for pt2_idx in range(pt1_idx+1, numpoints):
                    if (None != self.incl_pairs) and not self.incl_pairs.is_included_pair(self.point_ids[pt1_idx], self.point_ids[pt2_idx]):
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
        
        if self.is_network:
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
                    if (None != self.incl_pairs) and not self.incl_pairs.is_included_pair(self.points_rc[pt1_idx, 0], self.points_rc[pt2_idx, 0]):
                        continue
                    src = self.get_graph_node_idx(pt2_idx, node_map)
                    if (src >=  0 and components[src] == comp):
                        yield (pt1_idx, pt2_idx)
                yield(pt1_idx, -1)




class HabitatGraph:
    def __init__(self, g_map=None, poly_map=None, connect_using_avg_resistances=False, connect_four_neighbors_only=False, g_graph=None, node_names=None):
        if None != g_map:
            self.is_network = False
            self.g_map = g_map
            self.poly_map = poly_map
            self.connect_using_avg_resistances = connect_using_avg_resistances
            self.connect_four_neighbors_only = connect_four_neighbors_only
            
            self.node_map = HabitatGraph._construct_node_map(g_map, poly_map)
            (component_map, components) = HabitatGraph._construct_component_map(g_map, self.node_map, connect_using_avg_resistances, connect_four_neighbors_only)
            self.component_map = component_map
            self.components = components
            
            self.num_components = components.max()
            self.num_nodes = self.node_map.max()
        else:
            self.is_network = True
            self.g_graph = g_graph          # is the sparse CSR matrix 
            self.node_map = node_names    # list of node names
            
            (_num_components, C) = connected_components(g_graph)
            C += 1
            self.components = C
            
            self.num_components = C.max()
            self.num_nodes = self.node_map.size
    
    def get_graph(self):
        if self.is_network:
            return self.g_graph
        else:
            return HabitatGraph._construct_g_graph(self.g_map, self.node_map, self.connect_using_avg_resistances, self.connect_four_neighbors_only)
    
    def prune_for_component(self, keep_component, arr):
        if not self.is_network:
            raise RuntimeError('Not available in raster mode')
        indices = np.where(self.components == keep_component)
        return arr[indices]
    
    def prune_nodes_for_component(self, keep_component):
        """Removes nodes outside of component being operated on.
        
        Returns node map and adjacency matrix that only include nodes in keep_component.
        """
        if self.is_network:
            del_indices = np.where(self.components != keep_component)
            pruned_graph = ComputeBase.deleterowcol(self.g_graph, delrow=del_indices, delcol=del_indices)
            return (pruned_graph, self.prune_for_component(keep_component, self.node_map))
        else:
            selector = self.component_map == keep_component
            
            g_map_pruned = selector * self.g_map
            poly_map_pruned = []
            if self.poly_map !=  []:
                poly_map_pruned = selector * self.poly_map
    
            node_map_pruned = HabitatGraph._construct_node_map(g_map_pruned, poly_map_pruned)
            G_pruned = HabitatGraph._construct_g_graph(g_map_pruned, node_map_pruned, self.connect_using_avg_resistances, self.connect_four_neighbors_only)
             
            return (G_pruned, node_map_pruned)
            
    
    def unique_component_with_points(self, point_map):
        components_with_points = self.components_with_points(point_map)
        return None if (components_with_points.size > 1) else components_with_points[0]
        
    def components_with_points(self, point_map):
        if self.is_network:
            raise RuntimeError('Not available in network mode')
        
        sub_component_map = np.where(point_map, self.component_map, 0)
        sub_components = np.unique(sub_component_map)
        return sub_components if (sub_components[0] != 0) else sub_components[1:]

    
    @staticmethod
    @print_rusage              
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
        node_map[np.where(node_map)] = ComputeBase.relabel(node_map[np.where(node_map)], 1) #BHM 072409

        return node_map

    
    @staticmethod
    @print_rusage
    def _construct_component_map(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only):
        """Assigns component numbers to grid corresponding to pixels with non-zero conductances.
        
        Nodes with the same component number are in single, connected components.
        
        """  
        G = HabitatGraph._construct_g_graph(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only) 
        (_num_components, C) = connected_components(G)
        C += 1 # Number components from 1

        (I, J) = np.where(node_map)
        nodes = node_map[I, J].flatten()

        component_map = np.zeros(node_map.shape, dtype = 'int32')
        component_map[I, J] = C[nodes-1]

        return (component_map, C)
    
    @staticmethod
    @print_rusage
    def _construct_g_graph(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only):
        """Construct sparse adjacency matrix given raster maps of conductances and nodes."""
        numnodes = node_map.max()
        (node1, node2, conductances) = HabitatGraph._get_conductances(g_map, node_map, connect_using_avg_resistances, connect_four_neighbors_only)
        return HabitatGraph._make_sparse_csr(node1, node2, conductances, numnodes)
        
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
        (s_horiz, t_horiz) = HabitatGraph._neighbors_horiz(g_map)
        (s_vert,  t_vert)  = HabitatGraph._neighbors_vert(g_map)

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
            (s_dr, t_dr) = HabitatGraph._neighbors_diag1(g_map)
            (s_dl, t_dl) = HabitatGraph._neighbors_diag2(g_map)

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


class Output:
    """Handles output of current and voltage maps"""
    def __init__(self, options, state, report_status, node_names=None):
        self.options = options
        self.state = state
        self.report_status = report_status
        self.node_names = node_names
        self.nn_sorted = node_names[np.argsort(node_names)] if (None != node_names) else None
        
        self.is_network = (self.options.data_type == 'network')
        self.scenario = self.options.scenario
        self.voltage_maps = {}
        self.current_maps = {}


    def get_c_map(self, name, remove=False):
        """Get current map identified by name. Remove it from storage if remove is True."""
        ret = self.current_maps.get(name, None)
        if remove:
            self.rm_c_map(name)
        return ret
    
    def rm_c_map(self, name):
        """Remove current map identified by name if present"""
        if self.current_maps.has_key(name):
            del self.current_maps[name]

    def store_max_c_map(self, max_name, name, remove=False):
        """Store the maximum currents into map identified by max_name by comparing it with that identified by name.
        Implemented only for raster mode."""
        self.store_max_c_map_values(max_name, self.current_maps[name])
        if remove:
            self.rm_c_map(name)
            
    def store_max_c_map_values(self, max_name, values):
        self.current_maps[max_name] = np.maximum(self.current_maps[max_name], values)

    def accumulate_c_map_from(self, name, fromname):
        self.accumulate_c_map_with_values(name, self.current_maps[fromname])
    
    def accumulate_c_map_with_values(self, name, values):
        if self.is_network:
            _branch_currents, node_currents, branch_currents_array, node_map = values
            full_branch_currents, full_node_currents, _bca, _np = self.current_maps[name]
            pos = np.searchsorted(self.nn_sorted, node_map)
            full_node_currents[pos] += node_currents
            bc_x = np.searchsorted(self.nn_sorted, branch_currents_array[:,0])
            bc_y = np.searchsorted(self.nn_sorted, branch_currents_array[:,1])            
            full_branch_currents = full_branch_currents + sparse.csr_matrix((branch_currents_array[:,2], (bc_x, bc_y)), shape=full_branch_currents.shape)
            self.current_maps[name] = (full_branch_currents, full_node_currents, branch_currents_array, node_map)
        else:
            self.current_maps[name] += values
            
    def accumulate_c_map(self, name, voltages, G, node_map, finitegrounds, local_src, local_dst):
        self._write_store_c_map(name, False, False, True, voltages, G, node_map, finitegrounds, local_src, local_dst)

    def store_c_map(self, name, voltages, G, node_map, finitegrounds, local_src, local_dst):
        self._write_store_c_map(name, False, False, False, voltages, G, node_map, finitegrounds, local_src, local_dst)
    
    def write_c_map(self, name, remove=False, voltages=None, G=None, node_map=None, finitegrounds=None, local_src=None, local_dst=None):
        self._write_store_c_map(name, remove, True, False, voltages, G, node_map, finitegrounds, local_src, local_dst)
        
    def _write_store_c_map(self, name, remove, write, accumulate, voltages, G, node_map, finitegrounds, local_src, local_dst):
        if self.is_network:
            if voltages != None:
                (node_currents, branch_currents) = self._create_current_maps(voltages, G, finitegrounds)
                branch_currents_array = Output._convert_graph_to_3_col(branch_currents, node_map)
            else:
                branch_currents, node_currents, branch_currents_array, _node_map = self.current_maps[name]
            
            if write:                
                # Append node names and convert to array format
                node_currents_array = Output._append_names_to_node_currents(node_currents, node_map)                
                # node_currents_array = Output._append_names_to_node_currents(node_currents, self.nn_sorted)# Fixme: Need something like this for multipe components in network mode?
                CSIO.write_currents(self.options.output_file, branch_currents_array, node_currents_array, name, self.options)
                
            if remove:
                self.rm_c_map(name)
            elif accumulate:
                full_branch_currents, full_node_currents, _bca, _np = self.current_maps[name]
                pos = np.searchsorted(self.nn_sorted, node_map)
                full_node_currents[pos] += node_currents
                bc_x = np.searchsorted(self.nn_sorted, branch_currents_array[:,0])
                bc_y = np.searchsorted(self.nn_sorted, branch_currents_array[:,1])                
                full_branch_currents = full_branch_currents + sparse.csr_matrix((branch_currents_array[:,2], (bc_x, bc_y)), shape=full_branch_currents.shape)
                self.current_maps[name] = (full_branch_currents, full_node_currents, branch_currents_array, node_map)
            else:
                self.current_maps[name] = (branch_currents, node_currents, branch_currents_array, node_map)
        else:
            if voltages != None:
                while LowMemRetry.retry():
                    with LowMemRetry():
                        current_map = self._create_current_maps(voltages, G, finitegrounds, node_map)   
            else:
                current_map = self.current_maps[name]                                                                                
    
            if (self.options.set_focal_node_currents_to_zero==True) and (local_src is not None) and (local_dst is not None):
                # set source and target node currents to zero
                focal_node_pair_map = np.where(node_map == local_src+1, 0, 1)
                focal_node_pair_map = np.where(node_map == local_dst+1, 0, focal_node_pair_map)                                                
                current_map = np.multiply(focal_node_pair_map, current_map)
                del focal_node_pair_map
            
            if write:
                if name=='' and self.scenario != 'advanced': 
                    curMapName = 'cum_curmap' #For backward compatibility- ArcGIS
                else:
                    curMapName = 'curmap'
                fileadd = name if (name=='') else ('_'+name)
                CSIO.write_aaigrid(curMapName, fileadd, self._log_transform(current_map), self.options, self.state)
            if remove:
                self.rm_c_map(name)
            elif accumulate:
                self.current_maps[name] += current_map
            else:
                self.current_maps[name] = current_map


    def alloc_c_map(self, name):
        """Allocate space to store current maps identified by name"""
        if self.is_network:
            nnodes = len(self.node_names)
            branch_currents = sparse.csr_matrix((nnodes,nnodes))
            node_currents = np.zeros((nnodes, 1), dtype='float64')
            self.current_maps[name] = (branch_currents, node_currents, None, None)
        else:
            self.current_maps[name] = np.zeros((self.state.nrows, self.state.ncols), dtype='float64')


    def _log_transform(self, map_to_transform):
        if self.options.log_transform_maps:
            map_to_transform = np.where(map_to_transform > 0, np.log10(map_to_transform), self.state.nodata)
        return map_to_transform

    def get_v_map(self, name, remove=False):
        """Get voltage map identified by name. Remove it from storage if remove is True."""
        ret = self.voltage_maps.get(name, None)
        if remove:
            self.rm_v_map(name)
        return ret
    
    def rm_v_map(self, name):
        """Delete voltage map allocated with name if present"""
        if self.voltage_maps.has_key(name):
            del self.voltage_maps[name]
        
    def write_v_map(self, name, remove=False, voltages=None, node_map=None):
        """Write voltage map identified by name. Remove the map after that if remove is True.
        If voltages and node_map are provided, create space for voltage map disregarding prior allocation.
        """
        if self.report_status==True:
            ComputeBase.logger.info('writing voltage map ' + name)
            
        if self.is_network:
            if voltages == None:
                vm = self.voltage_maps[name]
                nn = self.node_names
            else:
                vm = voltages
                nn = node_map
            CSIO.write_voltages(self.options.output_file, vm, nn, name)
        else:
            fileadd = name if (name=='') else ('_'+name)
            if voltages == None:
                vm = self.voltage_maps[name]
            else:
                vm = self._create_voltage_map(node_map, voltages)
            CSIO.write_aaigrid('voltmap', fileadd, vm, self.options, self.state)
            
        if remove:
            self.rm_v_map(name)
        elif voltages != None:
            if self.is_network:
                self.alloc_v_map(name)
                self.accumulate_v_map(name, voltages, node_map)
            else:
                self.voltage_maps[name] = vm

    def accumulate_v_map(self, name, voltages, node_map):
        """Create and accumulate voltage map into the space identified by name"""
        if self.is_network:
            idxs = np.asarray([np.nonzero(self.node_names == a)[0][0] for a in node_map])
            vm = self.voltage_maps[name]
            vm[idxs] += voltages 
            self.voltage_maps[name] = vm
        else:
            self.voltage_maps[name] +=  self._create_voltage_map(node_map, voltages)
            
    def alloc_v_map(self, name):
        """Allocate space for a new voltage map with the given name"""
        if self.is_network:
            self.voltage_maps[name] = np.zeros(len(self.node_names), dtype='float64')
        else:
            self.voltage_maps[name] = np.zeros((self.state.nrows, self.state.ncols), dtype='float64')
        

    @gc_before
    @print_rusage
    def _create_current_maps(self, voltages, G, finitegrounds, node_map=None):
        """In raster mode, returns raster current map given node voltage vector, adjacency matrix, etc.
        In network mode returns node and branch currents given voltages in arbitrary graphs.
        """  
        G =  G.tocoo()
        node_currents = Output._get_node_currents(voltages, G, finitegrounds)
        
        if self.is_network:
            node_currents_col = np.zeros((node_currents.shape[0],1), dtype='float64')
            node_currents_col[:,0] = node_currents[:]
            branch_currents = Output._get_branch_currents(G, voltages, True) 
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
        node_currents_pos = Output._get_node_currents_posneg(G, voltages, finitegrounds, True) 
        node_currents_neg = Output._get_node_currents_posneg(G, voltages, finitegrounds, False)
        node_currents = np.where(node_currents_neg > node_currents_pos, node_currents_neg, node_currents_pos)

        return np.asarray(node_currents)[0]

    
    @staticmethod
    def _get_node_currents_posneg(G, voltages, finitegrounds, pos):
        """Calculates positive or negative node currents based on pos flag."""  
        branch_currents = Output._get_branch_currents(G, voltages, pos)
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
        branch_currents = Output._get_branch_currents_posneg(G, voltages, pos)
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

    @print_rusage
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
        output_node_currents[:,0] = node_names[:] # Fixme: fails in network mode when have multiple components
        try:
            output_node_currents[:,1] = node_currents[:,0]
        except:
            output_node_currents[:,1] = node_currents[:]
        return output_node_currents


