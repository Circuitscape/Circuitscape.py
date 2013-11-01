#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import time, logging
#from numpy import *
from scipy import sparse
import numpy as np

from cs_base import print_timing, CSBase, CSFocalPoints, CSHabitatGraph, CSOutput
from cs_io import CSIO
from cs_raster import CSRaster


class circuitscape(CSRaster):
    def __init__(self, configFile, logger_func):
        super(circuitscape, self).__init__(configFile, logger_func)

    @print_timing
    def compute(self):
        """Main function for Circuitscape."""  
        # Code below provides a back door to network mode, because not incorporated into GUI yet
        if self.options.polygon_file == 'NETWORK': 
            self.options.data_type='network' #can also be set in .ini file
        if self.options.data_type=='network': 
            self.options.graph_file = self.options.habitat_file
            self.options.focal_node_file = self.options.point_file

        self.state.start_time = time.time()
        self.state.last_gui_yield_time = time.time()        
        self.log('',1)
        self.log('',2)

        #Test write privileges by writing config file to output directory
        self.options.write(self.options.output_file, True)
        
        if self.options.data_type=='network':
            result, solver_failed = self.compute_network() # Call module for solving arbitrary graphs (not raster grids)
            self.log_complete_job()
            return result, solver_failed #Fixme: add in solver failed check

        return self.compute_raster()


    @print_timing
    def compute_network(self): 
        """Solves arbitrary graphs instead of raster grids."""
        (g_graph, node_names) = self.read_graph(self.options.graph_file)
        focal_nodes = self.readFocalNodes(self.options.focal_node_file)
        
        if self.options.use_included_pairs==True:
            self.state.included_pairs = CSIO.read_included_pairs(self.options.included_pairs_file)
        
        fp = CSFocalPoints(focal_nodes, self.state.included_pairs)
        g_habitat = CSHabitatGraph(g_graph=g_graph, node_names=node_names)
        cs = CSOutput(self.options, self.state, False, (g_habitat.num_nodes, g_habitat.num_nodes))
        
        #self.state.nrows = self.state.ncols = g_habitat.num_nodes
        (resistances, solver_failed) = self.single_ground_all_pair_resistances(g_habitat, fp, cs, True)
        _resistances, resistances_3col = self.writeResistances(fp.point_ids, resistances)
        return resistances_3col, solver_failed

    
    def read_graph(self, filename):
        """Reads arbitrary graph from disk. Returns sparse adjacency matrix and node names ."""
        graph_list = CSIO.load_graph(filename)

        try:
            zeros_in_resistance_graph = False           
            nodes = CSBase.deletecol(graph_list,2) 
            nodeNames = np.unique(np.asarray(nodes, dtype='int32'))
            nodes[np.where(nodes>= 0)] = CSBase.relabel(nodes[np.where(nodes>= 0)], 0)
            node1 = nodes[:,0]
            node2 = nodes[:,1]
            data = graph_list[:,2]
            
            ######################## Reclassification code
            if self.options.use_reclass_table == True:
                try:
                    reclassTable = CSIO.read_point_strengths(self.options.reclass_file)    
                except:
                    raise RuntimeError('Error reading reclass table.  Please check file format.')
                for i in range (0,reclassTable.shape[0]):
                    data = np.where(data==reclassTable[i,0], reclassTable[i,1],data)
                logging.debug('Reclassified habitat graph using %s'%(self.options.reclass_file,))
            ########################
            
            if self.options.habitat_map_is_resistances == True:
                zeros_in_resistance_graph = (np.where(data==0, 1, 0)).sum() > 0
                conductances = 1/data
            else:
                conductances = data
                
            numnodes = nodeNames.shape[0]
            G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes))

            Gdense=G.todense()
            g_graph = np.maximum(Gdense, Gdense.T) # To handle single or double entries for elements BHM 06/28/11
            g_graph = sparse.csr_matrix(g_graph)
        except:
            raise RuntimeError('Error processing graph/network file.  Please check file format')
        
        if zeros_in_resistance_graph == True:
            raise RuntimeError('Error: zero resistance values are not currently allowed in habitat network/graph input file.')
        
        return g_graph, nodeNames


    def readFocalNodes(self, filename):
        """Loads list of focal nodes for arbitrary graph."""  
        focal_nodes = CSIO.load_graph(filename)
        try:    
            if filename == self.options.graph_file:#If graph was used as focal node file, then just want first two columns for focal_nodes.
                focal_nodes = CSBase.deletecol(focal_nodes, 2)
            focal_nodes = np.unique(np.asarray(focal_nodes))
        except:
            raise RuntimeError('Error processing focal node file.  Please check file format')
        return focal_nodes


