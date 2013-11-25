__version__ = '4.0.0-beta'
__author__ = 'Brad McRae, Viral B. Shah, and Tanmay Mohapatra'
__email__ = 'mcrae@circuitscape.org'

import time
from scipy import sparse
import numpy as np

from cs_base import print_rusage, CSBase, CSFocalPoints, CSHabitatGraph, CSOutput
from cs_io import CSIO
from cs_raster import CSRaster
from cs_profiler import gc_after

class circuitscape(CSRaster):
    def __init__(self, configFile, ext_log_handler):
        super(circuitscape, self).__init__(configFile, ext_log_handler)

    @print_rusage
    def compute(self):
        """Main function for Circuitscape."""  
        # Code below provides a back door to network mode, because not incorporated into GUI yet
        if self.options.polygon_file == 'NETWORK': 
            self.options.data_type='network' #can also be set in .ini file
        if self.options.data_type=='network': 
            self.options.graph_file = self.options.habitat_file
            self.options.focal_node_file = self.options.point_file

        self.state.start_time = time.time()

        #Test write privileges by writing config file to output directory
        self.options.write(self.options.output_file, True)
        
        if self.options.data_type=='network':
            result, solver_failed = self.compute_network() # Call module for solving arbitrary graphs (not raster grids)
            self.log_complete_job()
            return result, solver_failed #Fixme: add in solver failed check

        return self.compute_raster()


    @print_rusage
    def compute_network(self): 
        """Solves arbitrary graphs instead of raster grids."""
        (g_graph, node_names) = self.read_graph(self.options.graph_file)
        focal_nodes = self.read_focal_nodes(self.options.focal_node_file)
        
        if self.options.use_included_pairs==True:
            self.state.included_pairs = CSIO.read_included_pairs(self.options.included_pairs_file)
        
        fp = CSFocalPoints(focal_nodes, self.state.included_pairs, True)
        g_habitat = CSHabitatGraph(g_graph=g_graph, node_names=node_names)
        cs = CSOutput(self.options, self.state, False, (g_habitat.num_nodes, g_habitat.num_nodes))
        if self.options.write_cur_maps:
            cs.alloc_c_map('')
        
        (resistances, solver_failed) = self.single_ground_all_pair_resistances(g_habitat, fp, cs, True)
        _resistances, resistances_3col = self.write_resistances(fp.point_ids, resistances)
        if self.options.write_cur_maps:
            full_branch_currents, full_node_currents, _bca, _np = cs.get_c_map('')
            full_branch_currents = CSOutput._convert_graph_to_3_col(full_branch_currents, node_names)
            full_node_currents = CSOutput._append_names_to_node_currents(full_node_currents, node_names)

            ind = np.lexsort((full_branch_currents[:, 1], full_branch_currents[:, 0]))
            full_branch_currents = full_branch_currents[ind]

            ind = np.lexsort((full_node_currents[:, 1], full_node_currents[:, 0]))
            full_node_currents = full_node_currents[ind]

            CSIO.write_currents(self.options.output_file, full_branch_currents, full_node_currents, '')
            
        return resistances_3col, solver_failed

    @gc_after
    def read_graph(self, filename):
        """Reads arbitrary graph from disk. Returns sparse adjacency matrix and node names ."""
        graph_list = CSIO.load_graph(filename)

        try:
            zeros_in_resistance_graph = False           
            nodes = CSBase.deletecol(graph_list,2) 
            node_names = np.unique(np.asarray(nodes, dtype='int32'))
            nodes[np.where(nodes>= 0)] = CSBase.relabel(nodes[np.where(nodes>= 0)], 0)
            node1 = nodes[:,0]
            node2 = nodes[:,1]
            data = graph_list[:,2]
            
            ######################## Reclassification code
            if self.options.use_reclass_table == True:
                try:
                    reclass_table = CSIO.read_point_strengths(self.options.reclass_file)    
                except:
                    raise RuntimeError('Error reading reclass table.  Please check file format.')
                for i in range (0,reclass_table.shape[0]):
                    data = np.where(data==reclass_table[i,0], reclass_table[i,1],data)
                circuitscape.logger.debug('Reclassified habitat graph using %s'%(self.options.reclass_file,))
            ########################
            
            if self.options.habitat_map_is_resistances == True:
                zeros_in_resistance_graph = (np.where(data==0, 1, 0)).sum() > 0
                conductances = 1/data
            else:
                conductances = data
                
            numnodes = node_names.shape[0]
            G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes))

            Gdense=G.todense()
            g_graph = np.maximum(Gdense, Gdense.T) # To handle single or double entries for elements BHM 06/28/11
            g_graph = sparse.csr_matrix(g_graph)
        except:
            raise RuntimeError('Error processing graph/network file.  Please check file format')
        
        if zeros_in_resistance_graph == True:
            raise RuntimeError('Error: zero resistance values are not currently allowed in habitat network/graph input file.')
        
        return g_graph, node_names


    def read_focal_nodes(self, filename):
        """Loads list of focal nodes for arbitrary graph."""  
        focal_nodes = CSIO.load_graph(filename)
        try:    
            if filename == self.options.graph_file:#If graph was used as focal node file, then just want first two columns for focal_nodes.
                focal_nodes = CSBase.deletecol(focal_nodes, 2)
            focal_nodes = np.unique(np.asarray(focal_nodes))
        except:
            raise RuntimeError('Error processing focal node file.  Please check file format')
        return focal_nodes


