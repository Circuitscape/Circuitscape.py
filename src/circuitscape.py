#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import sys, time, string, os, math
import ConfigParser
import gc
import traceback
import pdb

import numpy, scipy, pyamg

wx_available = True
try:
    import wx
except ImportError:
    wx_available = False

from sys import path
from string import split
from math import sqrt
from numpy import *
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from pyamg import *

from cs_util import *
from cs_cfg import CSConfig
import cs_io
#from cs_io import *
import copy

class circuitscape:
        
    def __init__(self, configFile, logger_func):
        gc.enable()
        self.state = {}
        self.state['amg_hierarchy'] = None
        self.options = CSConfig(configFile)
        if logger_func == 'Screen':
            self.options.screenprint_log = True
            logger_func = None
        else:
            self.options.screenprint_log = False
        numpy.seterr(invalid='ignore')
        numpy.seterr(divide='ignore')
        
        self.options.use_reclass_table = False
        self.options.reclass_file = './reclass.txt'        
        
        global logger
        logger = logger_func

        print_timing_enabled(self.options.print_timings)
        
    def log(self, text,col):
        """Prints updates to GUI or python window."""  
        (hours,mins,secs) = elapsed_time(self.state['startTime'])
        (hours,mins,secs1) = elapsed_time(self.state['lastUpdateTime'])
        if secs1 > 10: # Force update every 10 secs
            self.state['lastUpdateTime'] = time.time()
            try:
                if wx_available:
                    wx.SafeYield(None, True)  
                    wx.GetApp().Yield(True)
            except:
                pass
        if logger:           
            logger(text,col)
        if self.options.screenprint_log == True and len(text) > 1 and col == 1: 
            print '    --- ',text,' ---'
        sys.stdout.flush()
        return


    @print_timing
    def compute(self):
        """Main function for Circuitscape."""  

        # Code below provides a back door to network mode, because not incorporated into GUI yet
        if self.options.polygon_file == 'NETWORK': 
             self.options.data_type='network' #can also be set in .ini file
        if self.options.data_type=='network': 
            self.options.graph_file = self.options.habitat_file
            self.options.focal_node_file = self.options.point_file

        self.state['startTime'] = time.time()
        self.state['lastUpdateTime'] = time.time()        
        self.log('',1)
        self.log('',2)

        #Test write privileges by writing config file to output directory
        fileName = self.options.output_file
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        configFile = outputDir + '//' + outputBase + '.ini' 
        if os.path.isdir(outputDir):
            self.options.write(configFile)
            #cs_cfg.writeConfigFile(configFile, self.options)
        else:
            raise RuntimeError('Output directory ' + outputDir + ' does not exist!')    
        
        if self.options.data_type=='network':
            result,solver_failed = self.network_module() # Call module for solving arbitrary graphs (not raster grids)
            self.logCompleteJob()
            return result,solver_failed #Fixme: add in solver failed check
        else:
            self.load_maps()
            numNodes = (where(self.state['g_map'] > 0, 1, 0)).sum()         
            if self.options.screenprint_log == True:        
                print '    ---  Resistance/conductance map has',numNodes,' nodes.  ---'
        if self.options.scenario=='pairwise':
            resistances,solver_failed = self.pairwise_module(self.state['g_map'], self.state['poly_map'], self.state['points_rc'])
            self.logCompleteJob()
            return resistances,solver_failed     

        elif self.options.scenario == 'advanced':
            self.options.write_max_cur_maps=False
            self.log ('Calling solver module.',1)
            voltages,current_map,solver_failed = self.advanced_module(self.state['g_map'], self.state['poly_map'], self.state['source_map'], self.state['ground_map'],None,None,None,None,None)
            self.logCompleteJob()
            if solver_failed == True:
                print'Solver failed.\n'
            return voltages, solver_failed

        else:
            resistance_vector,solver_failed = self.one_to_all_module(self.state['g_map'], self.state['poly_map'], self.state['points_rc'])
            self.logCompleteJob()
            return resistance_vector,solver_failed 
            
    
    def grid_to_graph (self, x, y, node_map):
        """Returns node corresponding to x-y coordinates in input grid."""  
        return node_map[x, y] - 1
        
        
    def logCompleteJob(self):
        """Writes total time elapsed at end of run."""  
        (hours,mins,secs) = elapsed_time(self.state['startTime'])
        if hours>0:
            self.log('Job took ' + str(hours) +' hours ' + str(mins) + ' minutes to complete.',2)
        else:
            self.log('Job took ' + str(mins) +' minutes ' + str(secs) + ' seconds to complete.',2)
 
  
    def get_overlap_polymap(self,point,point_map,poly_map_temp,new_poly_num): 
        """Creates a map of polygons (aka short-circuit or zero resistance regions) overlapping a focal node."""  
        point_poly = where(point_map == point, 1, 0) 
        poly_point_overlap = multiply(point_poly,poly_map_temp)
        overlap_vals = unique(asarray(poly_point_overlap))
        (rows) = where(overlap_vals>0)
        overlap_vals = overlap_vals[rows] #LIST OF EXISTING POLYGONS THAT OVERLAP WITH POINT
        for a in range (0, overlap_vals.size):
            poly_map_temp = where(poly_map_temp==overlap_vals[a],new_poly_num,poly_map_temp)
        poly_map_temp = where(point_map == point, new_poly_num, poly_map_temp) 
        return poly_map_temp

    @print_timing
    def network_module(self): 
        """Solves arbitrary graphs instead of raster grids."""  
        solver_failed = False
        (g_graph,nodeNames) = self.read_graph(self.options.graph_file)
        (numComponents, C) = connected_components(g_graph)
        C += 1 # Number components from 1
        self.log('Graph has ' + str(g_graph.shape[0]) + ' nodes and '+ str(numComponents) + ' components.',2)
        if numComponents > 1:
            full_graph = g_graph
        else:
            nodesInComponent = nodeNames
        
        fullResistances = [] # initialize
        fullNodeCurrents = []
        for component in range (1,numComponents + 1):
            if numComponents > 1:
                indices = nodeNames[where(C == component)]
                delIndices = where(C != component)
                indices = where(C == component)
                nodesInComponent = nodeNames[indices]                
                g_graph = deleterowcol(full_graph, delrow = delIndices, delcol = delIndices)
                
            G = self.laplacian(g_graph)
            del g_graph
            if self.options.scenario=='advanced':
                raise RuntimeError('Advanced mode is not yet implemented in graph/network mode.')
                (sources,grounds)= self.readSourcesGroundsNetwork(G, nodeNames, self.options.source_file,self.options.ground_file)
                #FIXME: retool readSourceGroundNodes to read both files and return complete sources and grounds
                result,solver_failed = self.advanced_module_network(G,sources,grounds,nodeNames)           
                
            else:
                if self.options.use_included_pairs==True:
                    self.state['includedPairs'] = self.readincludedPairs(self.options.included_pairs_file)
                focalNodes = self.readFocalNodes(self.options.focal_node_file)
                if numComponents > 1:    #Prune out any focal nodes that are not in component
                    focalNodesInComponent = focalNodes
                    includeList = list(nodesInComponent[:])
                    numFocalNodes = focalNodesInComponent.shape[0]
                    row = 0
                    while row <numFocalNodes: 
                        if focalNodesInComponent[row] in includeList: #match
                            row = row+1
                        else:
                            n = focalNodesInComponent.shape[0]
                            keep = delete (arange(0, n), row)
                            focalNodesInComponent = focalNodesInComponent[keep]
                            numFocalNodes = numFocalNodes-1
                else:
                    focalNodesInComponent = focalNodes
                    numFocalNodes = focalNodes.shape[0]
                
                if self.options.scenario=='pairwise':
                    if numFocalNodes > 1:
                        # module returns arrays with node names
                        cumBranchCurrents,cumNodeCurrents,resistances3columns,solver_failed = self.pairwise_module_network(G,focalNodesInComponent,nodesInComponent,numComponents)
                    else: #nothing to solve
                        cumNodeCurrents = zeros((len(nodesInComponent),2),dtype = 'float64')
                        cumNodeCurrents[:,0] = nodesInComponent
                        resistances3columns = []
                        cumBranchCurrents = []
                    
                    # if first connected component solved, 
                    # then create output arrays
                    if (fullResistances == []) and (resistances3columns != []):
                        fullResistances = resistances3columns
                        if self.options.write_cur_maps == True:
                            fullBranchCurrents = cumBranchCurrents 
                            fullNodeCurrents = cumNodeCurrents
                    
                    # if ongoing solve, append results to ongoing output arrays
                    elif resistances3columns != []: 
                        fullResistances = append(fullResistances , resistances3columns, axis=0)
                        if self.options.write_cur_maps == True:
                            fullBranchCurrents = append(fullBranchCurrents , cumBranchCurrents, axis=0)
                            fullNodeCurrents = append(fullNodeCurrents , cumNodeCurrents, axis=0)
                    else: # If only node in component, just modify node current array
                        if self.options.write_cur_maps == True:
                            if fullNodeCurrents == []: #if first time populated
                                fullNodeCurrents = cumNodeCurrents
                            fullNodeCurrents = append(fullNodeCurrents , cumNodeCurrents, axis=0)
                    
                else:
                    raise RuntimeError('One-to-all and all-to-one modes are not yet implemented in graph/network mode.')
                    result,solver_failed = self.one_to_all_module_network(G,focalNodes,nodeNames)  

        if self.options.write_cur_maps == True:
            ind = lexsort((fullBranchCurrents[:, 1], fullBranchCurrents[:, 0]))
            fullBranchCurrents = fullBranchCurrents[ind]                        
                
            ind = lexsort((fullNodeCurrents[:, 1], fullNodeCurrents[:, 0]))
            fullNodeCurrents = fullNodeCurrents[ind]                        

            fileadd='cum'
            self.writeCurrentsNetwork(fullBranchCurrents, fullNodeCurrents, fileadd)

        #Make lists of focal node pairs.  Use to add disconnected pairs
        #to resistance output.
        resistancePairList = list()
        for row in range(0,fullResistances.shape[0]): 
            listEntry = str(int(fullResistances[row,0])) + "_" + str(int(fullResistances[row,1]))
            resistancePairList.append(listEntry)
            # print listEntry
        # print resistancePairList
        for sourceNode in range(0,len(focalNodes)-1):
            for targetNode in range(sourceNode+1,len(focalNodes)):                     
                pair = str(int(focalNodes[sourceNode])) + "_" + str(int(focalNodes[targetNode]))   
                # print pair
                if pair not in resistancePairList:
                    # Add in disconnected pair
                    newPair = array([[focalNodes[sourceNode],focalNodes[targetNode],-1]])
                    # print 'new',newPair
                    fullResistances = append(fullResistances , newPair, axis=0)
            
        ind = lexsort((fullResistances[:, 1], fullResistances[:, 0]))
        fullResistances = fullResistances[ind] 
        self.writeResistances3columns(fullResistances)
        
        return fullResistances,solver_failed #Fixme: need to check solver failed.

          
    def pairwise_module_network(self,G,focalNodes,nodeNames,numComponents):
        """Overhead module to solves arbitrary graphs in pairwise mode.
        
        Returns branch currents in 3-col format plus 3-column voltages,
        resistances, 2-col node currents.

        Writes currents voltages and for each pair.

        """
        solver_failed = False
        if self.options.use_included_pairs==True: #Prune points
            focalNodes = self.pruneIncludedPairsNetwork(focalNodes)          
            includedPairs = self.state['includedPairs'] 
        else:
            includedPairs = ones((focalNodes.size+1,focalNodes.size+1),dtype = 'int32')

        numpoints = focalNodes.size
        resistances = -1*ones((focalNodes.size,focalNodes.size),dtype = 'float64')
        if self.options.write_cur_maps == True:
            cumNodeCurrents = zeros((G.shape[0],1),dtype = 'float64')
            cumBranchCurrents = sparse.csr_matrix((G.shape))
        if (self.options.write_cur_maps == True) or (self.options.write_volt_maps == True) or (self.options.use_included_pairs==True):        
            useResistanceCalcShortcut = False
        else:
            useResistanceCalcShortcut = True
            # This uses a quicker way to calculate pairwise resistances.  There may be something even quicker, I'm looking into this. BHM 3-15-2013
            print 'Using shortcut'
            shortcutResistances = -1 * ones((numpoints, numpoints), dtype = 'float64') 
            voltmatrix = zeros((numpoints,numpoints),dtype = 'float64')     
        dstPoint = 0
        anchorPoint = 0            
        
        x = 0
        for i in range(0, numpoints):
            if range(i, numpoints) == []:
                break
            if (useResistanceCalcShortcut==True) and (dstPoint>0):
                break #No need to continue, we've got what we need to calculate resistances
            dstPoint = dstPoint+1            
            dst = self.nameToNode(nodeNames,focalNodes[i])
            G_dst_dst = G[dst, dst] 
            G[dst,dst] = 0
            self.state['amg_hierarchy'] = None
            gc.collect()
            self.create_amg_hierarchy(G)         
            
            for j in range(i+1, numpoints):
                x = x+1
                (hours,mins,secs) = elapsed_time(self.state['startTime'])
                if useResistanceCalcShortcut==True:
                    y = numpoints
                    self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(y) + '.',1)
                else:
                    y = numpoints*(numpoints-1)/2
                    self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)            
                if includedPairs[i+1,j+1]==1:
                    src = self.nameToNode(nodeNames,focalNodes[j])
                    try:
                        voltages = self.single_ground_solver(G, src, dst)

                        resistances[i, j] = voltages[src] - voltages[dst]
                        resistances[j, i] = voltages[src] - voltages[dst]
                    except:
                        solver_failed = True
                        resistances[i, j] = -777
                        resistances[j, i] = -777
    
                    if (useResistanceCalcShortcut==True) and (dstPoint==1) and (solver_failed==False): #this occurs for first i that is in component
                        anchorPoint = i #for use later in shortcut resistance calc
                        voltmatrix = self.getVoltmatrixNetwork(i,j,numpoints,voltages,resistances,voltmatrix,focalNodes,nodeNames) 
    
                    if (self.options.write_cur_maps == True) and (solver_failed==False):
                        finitegrounds = [-9999]
                        (node_currents,branch_currents) = self.getCurrentsNetwork(G,voltages,finitegrounds)
                        
                        
                        # Append node names and convert to array format
                        branchCurrentsArray = self.convert_graph_to_3_col(branch_currents,nodeNames)
                        nodeCurrentsArray = self.append_names_to_node_currents(node_currents, nodeNames)
                        
                        if self.options.write_cum_cur_map_only==False:                
                            fileadd=str(int(focalNodes[i])) + '_' + str(int(focalNodes[j]))
                            self.writeCurrentsNetwork(branchCurrentsArray, nodeCurrentsArray, fileadd)

                        cumNodeCurrents=cumNodeCurrents+node_currents
                        cumBranchCurrents=cumBranchCurrents+branch_currents                       
    
                    if (self.options.write_volt_maps == True) and (solver_failed==False):
                        fileadd=str(int(focalNodes[i])) + '_' + str(int(focalNodes[j]))
                        self.writeVoltagesNetwork(voltages, nodeNames, fileadd)
    
                    if solver_failed==True:
                        solver_failed = False
                        solver_failed_somewhere = True                 
            if (useResistanceCalcShortcut==True) and (i==anchorPoint): #this happens once per component. Anchorpoint is the first i in component
                shortcutResistances = self.getShortcutResistances(anchorPoint,voltmatrix,numpoints,resistances,shortcutResistances)
             
            G[dst,dst] = G_dst_dst
             
        if self.options.write_cur_maps == True:
            cumBranchCurrentsArray = self.convert_graph_to_3_col(cumBranchCurrents,nodeNames)
            cumNodeCurrentsArray =  self.append_names_to_node_currents(cumNodeCurrents, nodeNames)
        else:
            cumBranchCurrentsArray = -1
            cumNodeCurrentsArray = -1

        if useResistanceCalcShortcut==True:
            resistances = shortcutResistances
        for i in range(0,numpoints):
            resistances[i, i] = 0

        outputResistances = self.append_names_to_resistances(focalNodes, resistances)       
        resistances3columns = self.convertResistances3cols(outputResistances) 
               
        return cumBranchCurrentsArray,cumNodeCurrentsArray,resistances3columns,solver_failed  


        
    def append_names_to_node_currents(self,node_currents, nodeNames):
        """Adds names of focal nodes to node current lists."""    
        outputNodeCurrents=zeros((len(node_currents),2),dtype='float64')
        outputNodeCurrents[:,0]=nodeNames[:]
        try:
            outputNodeCurrents[:,1]=node_currents[:,0]
        except:
            outputNodeCurrents[:,1]=node_currents[:]
        return outputNodeCurrents
        
        
    def append_names_to_resistances(self, point_ids, resistances):        
        """Adds names of focal nodes to resistance matrices."""  
        focal_labels = insert(point_ids, [0], 0, axis = 0)
        resistances = insert(resistances, [0], 0, axis = 0)
        resistances = insert(resistances, [0], 0, axis = 1)
        resistances[0,:] = focal_labels
        resistances[:,0] = focal_labels
        return resistances
        
        
    def writeCurrentsNetwork(self, branch_currents, node_currents, fileadd):
        """Writes currents from network operations.
        
           Inputs are arrays with node names.

        """       
        outputDir, outputFile = os.path.split(self.options.output_file)
        outputBase, outputExtension = os.path.splitext(outputFile)
        if branch_currents!=None:
            file = outputDir + '//' + outputBase + '_branch_currents_' + fileadd + '.txt'
            savetxt(file,branch_currents)
        file = outputDir + '//' + outputBase + '_node_currents_' + fileadd + '.txt'
        savetxt(file,node_currents)      
        return 


    def writeVoltagesNetwork(self, voltages, nodeNames, fileadd):
        """Saves voltage values from solving arbitrary graphs to disk."""  
        outputVoltages=zeros((len(voltages),2),dtype='float64')
        if nodeNames != None:
             outputVoltages[:,0]=nodeNames[:]
        outputVoltages[:,1]=voltages[:]
        outputDir, outputFile = os.path.split(self.options.output_file)
        outputBase, outputExtension = os.path.splitext(outputFile)
        file = outputDir + '//' + outputBase + '_voltages_' + fileadd + '.txt'
        savetxt(file,outputVoltages)      
        return        

        
    def getCurrentsNetwork(self,G,voltages,finitegrounds):
        """Returns node and branch currents given voltages in arbitrary graphs."""  
        G =  G.tocoo()
        node_currents = self.get_node_currents(voltages, G, finitegrounds)
        node_currents_col = zeros((node_currents.shape[0],1),dtype = 'float64')
        node_currents_col[:,0] = node_currents[:]
        branch_currents = self.get_branch_currents(G,voltages,True) 
        branch_currents = absolute(branch_currents) 
        G = G.tocsr()
        return node_currents_col,branch_currents
    
    
    def getVoltmatrixNetwork(self,i,j,numpoints,voltages,resistances,voltmatrix,focalNodes,nodeNames):                                            
        """Returns a matrix of pairwise voltage differences between focal nodes when operating on arbitrary graphs.
        
        Used for shortcut calculations of effective resistance when no
        voltages or currents are mapped.
        
        """  
        voltvector = zeros((numpoints,1),dtype = 'float64')  
        for point in range(1,numpoints):
           node=self.nameToNode(nodeNames,focalNodes[point])
           voltageAtPoint = voltages[node] 
           voltageAtPoint = 1-(voltageAtPoint/resistances[i, j])
           voltvector[point] = voltageAtPoint
        voltmatrix[:,j] = voltvector[:,0] 
        return voltmatrix             
             
             
    def nameToNode(self, nodeNames,name):
        """Returns node index given node ID."""  
        nodeNames = nodeNames.tolist()
        node = nodeNames.index(name)
        return node


    def namesToNodes(self, nodeNames,names):
        """Returns node indices given node IDs."""  
        nodeNames = nodeNames.tolist()
        nodes = zeros(len(names),dtype = 'int32')

        for i in range (0,len(names)):
           nodes[i] = nodeNames.index(names[i])
        return nodes

        
    def read_graph(self, filename):
        """Reads arbitrary graph from disk. Returns sparse adjacency matrix and node names ."""  
        graphList = self.load_graph(filename,datatype='float64')

        try:
            zeros_in_resistance_graph = False           
            nodes = deletecol(graphList,2) 
            nodeNames = unique(asarray(nodes))
            nodes[where(nodes>= 0)] = relabel(nodes[where(nodes>= 0)], 0)
            node1 = nodes[:,0]
            node2 = nodes[:,1]
            data = graphList[:,2]
            ######################## Reclassification code
            if self.options.use_reclass_table==True:
                try:
                    reclassTable = self.readPointStrengths(self.options.reclass_file)    
                except:
                    raise RuntimeError('Error reading reclass table.  Please check file format.')
                for i in range (0,reclassTable.shape[0]):
                    data = where(data==reclassTable[i,0],reclassTable[i,1],data)
                print'\n***** Reclassified habitat graph using', self.options.reclass_file,'*****'
            ########################            
            if self.options.habitat_map_is_resistances == True:
                zeros_in_resistance_graph = (where(data==0, 1, 0)).sum() > 0
                conductances = 1/data
            else:
                conductances = data
            numnodes = nodeNames.shape[0]
            G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes))

            Gdense=G.todense()
            g_graph = maximum(Gdense, Gdense.T) # To handle single or double entries for elements BHM 06/28/11
            g_graph = sparse.csr_matrix(g_graph)
    
        except:
            raise RuntimeError('Error processing graph/network file.  Please check file format')
        if zeros_in_resistance_graph == True:
            raise RuntimeError('Error: zero resistance values are not currently allowed in habitat network/graph input file.')
        return g_graph, nodeNames


    def readFocalNodes(self, filename):
        """Loads list of focal nodes for arbitrary graph."""  
        focalNodes = self.load_graph(filename,datatype='int32')
        try:    
            if filename==self.options.graph_file:#If graph was used as focal node file, then just want first two columns for focalNodes.
                focalNodes = deletecol(focalNodes, 2)
            focalNodes = unique(asarray(focalNodes))
        except:
            raise RuntimeError('Error processing focal node file.  Please check file format')
        return focalNodes
        
        
    def load_graph(self,filename,datatype):
        """Returns data for arbitrary graph or focal node list from file."""  
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        f = open(filename, 'r')

        try:
            graphObject = loadtxt(filename, dtype = 'Float64', comments='#') 
        except:
            try:
                graphObject = loadtxt(filename, dtype = 'Float64', comments='#', delimiter=',')
            except:
                raise RuntimeError('Error reading',type,'file.  Please check file format')
        return graphObject
        

    @print_timing
    def one_to_all_module(self, g_map, poly_map, points_rc):
        """Overhead module for one-to-all AND all-to-one modes with raster data."""  
        lastWriteTime = time.time()

        if self.options.use_included_pairs==True: #Prune points
            points_rc = self.pruneIncludedPairs(points_rc)          
            includedPairs = self.state['includedPairs'] 
        point_ids = unique(asarray(points_rc[:,0]))
        points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)      
        
        resistance_vector = zeros((point_ids.size,2),float)
        solver_failed_somewhere = False
        if self.options.write_cur_maps == True:
            cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')         
        if self.options.write_max_cur_maps==True:
            max_current_map=cum_current_map
        oneToAllStreamline = False        
        if self.options.use_included_pairs==False: #Will do this each time later if using included pairs
            point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]
       
            #combine point map and poly map
            poly_map_temp = self.get_poly_map_temp(poly_map,point_map,point_ids,None,None)
            unique_point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            unique_point_map[points_rc_unique[:,1],points_rc_unique[:,2]] = points_rc_unique[:,0]

            (strengthMap,strengths_rc) = self.getstrengthMap(points_rc_unique,self.state['pointStrengths'])            

            # Are all points in same component?  If so, can streamline.  Create G here, just once
            # FIXME: place code below into module 
            node_map = self.construct_node_map(g_map, poly_map_temp) #******WANT POLYS BURNED IN 
            (component_map, components) = self.construct_component_map(g_map, node_map)
            pointComponents = where(unique_point_map, component_map, 0)#TODO use this to choose component to operate on below, unless all points in same component
            self.log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max()) + ' components.',2)
            del node_map
            uniquePointComponentList = unique(pointComponents)
            numComponentsWithPoints = uniquePointComponentList.size
            if uniquePointComponentList[0]==0:
                numComponentsWithPoints = numComponentsWithPoints-1
            if numComponentsWithPoints==1:
                componentWithPoints = max(uniquePointComponentList)
                (G, node_map) = self.node_pruner(g_map, poly_map_temp, component_map, componentWithPoints) #FIXME drops polymaptemp node 1... #We'll keep node_map around in this case.  Only need to create it and G once.
                oneToAllStreamline = True
            else:
                del component_map, components  
        else:
            node_map = self.construct_node_map(g_map, poly_map) #******WANT POLYS BURNED IN 
            (component_map, components) = self.construct_component_map(g_map, node_map)
            self.log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max()) + ' components.',2)
            del node_map, component_map, components

        for i in range(0, point_ids.size): #These are the 'src' nodes, i.e. the 'one' in all-to-one and one-to-all
            (hours,mins,secs) = elapsed_time(self.state['startTime'])
            self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(i+1) + ' of ' + str(point_ids.size) + '.',1)

            if self.options.use_included_pairs==True: # Done above otherwise    
                #######################   
                points_rc_unique_temp = copy.copy(points_rc_unique)
                point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
                point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]       

                for pair in range(0, point_ids.size): #loop thru exclude[point,:], delete included pairs of focal point from point_map and points_rc_unique_temp 
                    if includedPairs[i+1,pair+1]==0 and i !=  pair:
                            id = point_ids[pair]
                            point_map = where(point_map==id,0,point_map)
                            points_rc_unique_temp[pair,0] = 0 #point will not be burned in to unique_point_map

                poly_map_temp = self.get_poly_map_temp2(poly_map,point_map,points_rc_unique_temp,includedPairs,i)

                unique_point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
                unique_point_map[points_rc_unique_temp[:,1],points_rc_unique_temp[:,2]] = points_rc_unique_temp[:,0]        

                (strengthMap,strengths_rc) = self.getstrengthMap(points_rc_unique_temp,self.state['pointStrengths'])
                ###########################################
                
            src = point_ids[i]
            sumConnections = unique_point_map.sum()
            sumConnections = sumConnections.sum()
            if sumConnections!= point_ids[i]: #If there are points to connect with src point
                if self.options.scenario == 'one-to-all':
                    if self.options.use_variable_source_strengths==True:
                        strength = strengths_rc[i,0]
                    else:
                        strength = 1
                    source_map = where(unique_point_map == src,strength,0)
                    ground_map = point_map
                    ground_map =  where(unique_point_map == src,0,unique_point_map)
                    ground_map =  where(ground_map,Inf,0) 
                    self.options.remove_src_or_gnd = 'rmvgnd'
                else:
                    if self.options.use_variable_source_strengths==True:
                        source_map = where(unique_point_map == src,0,strengthMap)
                    else:
                        source_map = where(unique_point_map,1,0)
                        source_map = where(unique_point_map == src,0,source_map)
                        
                    ground_map = point_map
                    ground_map =  where(unique_point_map == src,1,0)
                    ground_map =  where(ground_map,Inf,0)   

                    self.options.remove_src_or_gnd = 'rmvsrc'
                #FIXME: right now one-to-all *might* fail if there is just one node that is not grounded (I haven't encountered this problem lately BHM Nov 2009).

                if oneToAllStreamline==False:
                    G = None
                    node_map = None
                    component_map = None
                    componentWithPoints = None
                    
                resistance,current_map,solver_failed = self.advanced_module(self.state['g_map'], poly_map_temp, source_map, ground_map, src, G, node_map, component_map, componentWithPoints)
                if solver_failed == False:
                    if self.options.write_cur_maps == True:
                        cum_current_map = cum_current_map+current_map
                        if self.options.write_max_cur_maps==True:
                            max_current_map=maximum(max_current_map,current_map)
                else:
                    print'Solver failed for at least one focal node.  \nFocal nodes with failed solves will be marked with value of -777 \nin output resistance list.\n'
    
                resistance_vector[i,0] = src
                resistance_vector[i,1] = resistance
                    
                if solver_failed==True:
                    solver_failed_somewhere = True
            else:
                resistance_vector[i,0] = src
                resistance_vector[i,1] = -1            

            (hours,mins,secs) = elapsed_time(lastWriteTime)
            if secs > 120: 
                lastWriteTime = time.time()
                self.writeResistancesOneToAll(resistance_vector,'_incomplete')
      
        if solver_failed_somewhere==False:
            if self.options.write_cur_maps == True:
                self.write_aaigrid('cum_curmap', '', cum_current_map)
                if self.options.write_max_cur_maps==True:
                    self.write_aaigrid('max_curmap', '', max_current_map)

        self.writeResistancesOneToAll(resistance_vector,'')
       
        return resistance_vector,solver_failed_somewhere 

    def get_poly_map_temp(self,poly_map,point_map,point_ids,includedPairs,point1):
        """Returns polygon map for each solve given source and destination nodes.  
        
        Used in all-to-one and one-to-all modes only.
        
        """  
        if poly_map == []:
            poly_map_temp = point_map
        else:
            poly_map_temp = poly_map
            new_poly_num = numpy.max(poly_map)
            for point2 in range(0, point_ids.shape[0]):
                if (self.options.use_included_pairs==True) and (includedPairs[point1+1,point2+1] == 1):
                   new_poly_num = new_poly_num+1
                   poly_map_temp = self.get_overlap_polymap(point_ids[point2],point_map,poly_map_temp, new_poly_num) 
                else: 
                    new_poly_num = new_poly_num+1
                    poly_map_temp = self.get_overlap_polymap(point_ids[point2],point_map,poly_map_temp, new_poly_num)                 
                
        return poly_map_temp

        
    def get_poly_map_temp2(self,poly_map,point_map,points_rc_unique_temp,includedPairs,i):
        """Returns polygon map for each solve given source and destination nodes.  
        
        Used in all-to-one and one-to-all modes when included/excluded pairs are used.
        
        """  
        if poly_map == []:
            poly_map_temp = point_map
        else:
            poly_map_temp = poly_map
            new_poly_num = numpy.max(poly_map)
            #burn in src point to polygon map
            poly_map_temp = self.get_overlap_polymap(points_rc_unique_temp[i,0],point_map,poly_map_temp, new_poly_num)                     
            for point in range(0, points_rc_unique_temp.shape[0]): #burn in dst points to polygon map
                if includedPairs[i+1,point+1] == 1:  
                   new_poly_num = new_poly_num+1
                   poly_map_temp = self.get_overlap_polymap(points_rc_unique_temp[point,0],point_map,poly_map_temp, new_poly_num) 
        return poly_map_temp


    def pairwise_module(self, g_map, poly_map, points_rc):
        """Overhead module for pairwise mode with raster data."""  

        if self.options.write_cur_maps == True:
            cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
            if self.options.write_max_cur_maps==True:
                max_current_map=cum_current_map
            else:
                max_current_map=[]
        else:
            cum_current_map = []
            max_current_map=[]

         
        # If there are no focal regions, pass all points to single_ground_all_pair_resistances,
        # otherwise, pass one point at a time.
        if self.options.point_file_contains_polygons==False:
            if points_rc.shape[0]!= (unique(asarray(points_rc[:,0]))).shape[0]:
                raise RuntimeError('At least one focal node contains multiple cells.  If this is what you really want, then choose focal REGIONS in the pull-down menu') 
            
            else:
                if self.options.use_included_pairs==True:
                    points_rc = self.pruneIncludedPairs(points_rc)
#                     includedPairs = self.state['includedPairs']
#                 else:
#                     numpoints = points_rc.shape[0]
#                     includedPairs = ones((numpoints+1,numpoints+1),dtype = 'int32')
                
                reportStatus = True
                try:
                    (resistances,cum_current_map,max_current_map,solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map, points_rc,cum_current_map,max_current_map, reportStatus)
                except MemoryError: #Give it a try, but starting again never seems to helps even from GUI.
                    self.enable_low_memory(True) #This doesn't seem to really clear out memory or truly restart.
                    if self.options.write_cur_maps == True:
                        cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
                        if self.options.write_max_cur_maps==True:
                            max_current_map=cum_current_map
                        else:
                            max_current_map=[]
                    else:
                        cum_current_map = []
                        max_current_map=[]
                        
                    #Note: This does not go through when it should.
                    (resistances,cum_current_map,max_current_map,solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map, points_rc,cum_current_map,max_current_map,reportStatus)
                if solver_failed == True:
                    print('Solver failed for at least one focal node pair. ' 
                    '\nThis can happen when input resistances differ by more than' 
                    '\n~6 orders of magnitude. Pairs with failed solves will be '
                    '\nmarked with value of -777 in output resistance matrix.\n')

                point_ids = points_rc[:,0]

        else:
            if self.options.use_included_pairs==True:
                points_rc = self.pruneIncludedPairs(points_rc)
                includedPairs = self.state['includedPairs']
            else:
                numpoints = points_rc.shape[0]
                point_ids = unique(asarray(points_rc[:,0]))
                points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)#Nov13_2010   #Fixme: can just use point ids for index size
                numUniquepoints= points_rc_unique.shape[0]#Nov13_2010
                includedPairs = ones((numUniquepoints+1,numUniquepoints+1),dtype = 'int32')#Nov13_2010

            point_map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]

            point_ids = unique(asarray(points_rc[:,0]))
            points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)

            resistances = -1*ones((point_ids.size,point_ids.size),dtype = 'float64')
            x = 0
            for i in range(0, point_ids.size-1):
                for j in range(i+1, point_ids.size):
                    if includedPairs[i+1,j+1]==1:                
                        if poly_map == []:
                            poly_map_temp = zeros((self.state['nrows'],self.state['ncols']),int)
                            new_poly_num = 1
                        else:
                            poly_map_temp = poly_map
                            new_poly_num = numpy.max(poly_map)+1
                        point = point_ids[i]
                        poly_map_temp = self.get_overlap_polymap(point_ids[i],point_map,poly_map_temp, new_poly_num) 
                        poly_map_temp = self.get_overlap_polymap(point_ids[j],point_map,poly_map_temp, new_poly_num+1) 
    
                        #Get first instance of each point in points_rc
                        points_rc_temp = zeros((2,3),int)
                        points_rc_temp[0,:] = points_rc_unique[i,:]
                        points_rc_temp[1,:] = points_rc_unique[j,:]
    
                        numpoints = point_ids.size
                        x = x+1
                        y = numpoints*(numpoints-1)/2
                        (hours,mins,secs) = elapsed_time(self.state['startTime'])
                        self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
                        reportStatus = False
                        
                        (pairwise_resistance,cum_current_map,max_current_map,solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map_temp, points_rc_temp,cum_current_map,max_current_map,reportStatus)
    
                        del poly_map_temp
                        if solver_failed == True:
                            print'Solver failed for at least one focal node pair.  \nPairs with failed solves will be marked with value of -777 \nin output resistance matrix.\n'
    
                        resistances[i,j] = pairwise_resistance[0,1]
                        resistances[j,i] = pairwise_resistance[0,1]

        for i in range(0,resistances.shape[0]): #Set diagonal to zero
            resistances[i, i] = 0

        #Add row and column headers and write resistances to disk
        resistances = self.writeResistances(point_ids, resistances)

        if self.options.write_cur_maps == True:
            if solver_failed==False:
                if self.options.log_transform_maps == True:
                    cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state['nodata']) 
                    if self.options.write_max_cur_maps==True:
                        max_current_map = where(max_current_map>0,log10(max_current_map),self.state['nodata']) 
                self.write_aaigrid('cum_curmap', '', cum_current_map)
                if self.options.write_max_cur_maps==True:      
                    self.write_aaigrid('max_curmap', '', max_current_map)

        return resistances,solver_failed
    
    
    @print_timing
    def single_ground_all_pair_resistances(self, g_map, poly_map, points_rc,cum_current_map,max_current_map,reportStatus):
        """Handles pairwise resistance/current/voltage calculations.  
        
        Called once when focal points are used, called multiple times when focal regions are used.

        """  
        lastWriteTime = time.time()
        numpoints = points_rc.shape[0]
        if (self.options.use_included_pairs==False) or (self.options.point_file_contains_polygons==True):
            includedPairs = ones((numpoints+1,numpoints+1),dtype = 'int32')
        else:
            includedPairs = self.state['includedPairs']
        
        if (self.options.point_file_contains_polygons==True) or  (self.options.write_cur_maps == True) or (self.options.write_volt_maps == True) or (self.options.use_included_pairs==True): 
           useResistanceCalcShortcut = False
        else:     
           useResistanceCalcShortcut = True # We use this when there are no focal regions.  It saves time when we are also not creating maps
           shortcutResistances = -1 * ones((numpoints, numpoints), dtype = 'float64') 
           
        solver_failed_somewhere = False
        node_map = self.construct_node_map(g_map, poly_map) # Polygons burned in to node map 
        (component_map, components) = self.construct_component_map(g_map, node_map)
        
        self.log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max())+ ' components.',2)
        resistances = -1 * ones((numpoints, numpoints), dtype = 'float64')         #Inf creates trouble in python 2.5 on Windows. Use -1 instead.
        
        x = 0
        for c in range(1,int(components.max()+1)):
            points_in_this_component = self.checkPointsInComponent(c,numpoints,components,points_rc,node_map)
                        
            if points_in_this_component:
                (G, local_node_map) = self.node_pruner(g_map, poly_map, component_map, c)
                if c==int(components.max()):
                    del component_map 
                if (useResistanceCalcShortcut==True):
                    voltmatrix = zeros((numpoints,numpoints),dtype = 'float64')     #For resistance calc shortcut
                
                dstPoint = 0
                anchorPoint = 0 #For resistance calc shortcut
                
                ##############
                for i in range(0, numpoints):
                    if range(i, numpoints) == []:
                        break

                    if (useResistanceCalcShortcut==True) and (dstPoint>0): 
                        break #No need to continue, we've got what we need to calculate resistances

                    dst = self.grid_to_graph (points_rc[i,1], points_rc[i,2], node_map)
                    local_dst = self.grid_to_graph (points_rc[i,1], points_rc[i,2], local_node_map)
                    if (dst >=  0 and components[dst] == c):
                        dstPoint = dstPoint+1
                        Gsolve = []
                    
                        G_dst_dst = G[local_dst, local_dst] 
                        G[local_dst,local_dst] = 0
    
                        self.state['amg_hierarchy'] = None
                        gc.collect()
                        self.create_amg_hierarchy(G)

                        ################    
                        for j in range(i+1, numpoints):
                            if includedPairs[i+1,j+1]==1: #Test for pair in includedPairs
                                if self.state['amg_hierarchy']==None: #Called in case of memory error in current mapping
                                    self.create_amg_hierarchy(G)
                               
                                if reportStatus==True:
                                    x = x+1
                                    (hours,mins,secs) = elapsed_time(self.state['startTime'])
                                    if useResistanceCalcShortcut==True:
                                        y = numpoints
                                        self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(y) + '.',1)
                                    else:
                                        y = numpoints*(numpoints-1)/2
                                        self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
                                src = self.grid_to_graph (points_rc[j,1], points_rc[j,2], node_map)
                                local_src = self.grid_to_graph (points_rc[j,1], points_rc[j,2], local_node_map)
                                if (src >=  0 and components[src] == c):
                                    # Solve resistive network (Kirchoff's laws)
                                    solver_failed = False
                                    try:
                                        voltages = self.single_ground_solver(G, local_src, local_dst)

                                    except:
                                        solver_failed = True
                                        solver_failed_somewhere = True
                                        resistances[i, j] = -777
                                        resistances[j, i] = -777
        
                                    if solver_failed == False:
                                        if self.options.low_memory_mode==True or self.options.point_file_contains_polygons==True:
                                            self.state['amg_hierarchy'] = None
                                            gc.collect()    
                                        
                                        resistances[i, j] = voltages[local_src] - voltages[local_dst]
                                        resistances[j, i] = voltages[local_src] - voltages[local_dst]
                                        # Write maps to files
                                        frompoint = str(points_rc[i,0])
                                        topoint = str(points_rc[j,0])
                                        
                                        if useResistanceCalcShortcut==True:
                                            if dstPoint==1: #this occurs for first i that is in component
                                                anchorPoint = i #for use later in shortcult resistance calc
                                                voltmatrix = self.getVoltmatrix(i,j,numpoints,local_node_map,voltages,points_rc,resistances,voltmatrix)                                          

                                        if self.options.write_volt_maps == True:
                                            if reportStatus==True:
                                                (hours,mins,secs) = elapsed_time(self.state['startTime'])
                                                self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min writing voltage map ' + str(x) + ' of ' + str(y) + '.',1)
                                            voltage_map = self.create_voltage_map(local_node_map,voltages) 
                                            self.write_aaigrid('voltmap', '_' + frompoint + '_' + topoint, voltage_map)
                                            del voltage_map
                                        if self.options.write_cur_maps == True:
                                            if reportStatus==True:
                                                (hours,mins,secs) = elapsed_time(self.state['startTime'])
                                                self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min writing current map ' + str(x) + ' of ' + str(y) + '.',1)
                                            finitegrounds = [-9999] #create dummy value for pairwise case
                                            
                                            try:
                                                G = G.tocoo()#Can cause memory error
                                            except MemoryError:
                                                self.enable_low_memory(False)
                                                G = G.tocoo()
                                            try:
                                                current_map = self.create_current_map(voltages, G, local_node_map, finitegrounds)   
                                            except MemoryError:
                                                self.enable_low_memory(False)
                                                current_map = self.create_current_map(voltages, G, local_node_map, finitegrounds)                                                                                    
                                            G = G.tocsr()

                                            if self.options.set_focal_node_currents_to_zero==True:
                                                # set source and target node currents to zero
                                                focal_node_pair_map = where(local_node_map == local_src+1, 0, 1)
                                                focal_node_pair_map = where(local_node_map == local_dst+1, 0, focal_node_pair_map)                                                
                                                #print'fn', focal_node_pair_map
                                                current_map = multiply(focal_node_pair_map, current_map)
                                                #print 'c',current_map
                                                del focal_node_pair_map
                                            cum_current_map = cum_current_map + current_map 
                                            if self.options.write_max_cur_maps==True:
                                                max_current_map = maximum(max_current_map, current_map) 
                                            if self.options.write_cum_cur_map_only==False:
                                                if self.options.log_transform_maps==True:
                                                   current_map = where(current_map>0,log10(current_map),self.state['nodata'])
                                                self.write_aaigrid('curmap', '_' + frompoint + '_' + topoint, current_map)
                                            del current_map    

                                        (hours,mins,secs) = elapsed_time(lastWriteTime)
                                        if secs > 120: 
                                            lastWriteTime = time.time()
                                            self.saveIncompleteResistances(resistances)# Save incomplete resistances
                        if (useResistanceCalcShortcut==True and i==anchorPoint): # This happens once per component. Anchorpoint is the first i in component
                            shortcutResistances = self.getShortcutResistances(anchorPoint,voltmatrix,numpoints,resistances,shortcutResistances)
                                                
                        G[local_dst, local_dst] = G_dst_dst

                    #End for
                    self.state['amg_hierarchy'] = None
                    gc.collect()

                #End if
                self.state['amg_hierarchy'] = None
                gc.collect()

        # Finally, resistance to self is 0.
        if useResistanceCalcShortcut==True: 
            resistances = shortcutResistances
        for i in range(0,numpoints):
            resistances[i, i] = 0

        return resistances,cum_current_map,max_current_map,solver_failed_somewhere

        
    @print_timing
    def single_ground_solver(self, G, src, dst):
        """Solver used for pairwise mode."""  
        n = G.shape[0]
        rhs = zeros(n, dtype = 'float64')
        if src==dst:
            voltages = zeros(n, dtype = 'float64')
        else:
            rhs[dst] = -1
            rhs[src] = 1
            voltages = self.solve_linear_system (G, rhs)

        return voltages

        
    @print_timing
    def advanced_module(self, g_map, poly_map, source_map, ground_map,source_id, G, node_map, component_map, componentWithPoints): 
        """Overhead module for advanced mode with raster data."""  
        if node_map==None:
            oneToAllStreamline = False
        else:
            oneToAllStreamline = True
        solver_called = False
        solver_failed = False 

        if oneToAllStreamline==False:
            node_map = self.construct_node_map(g_map, poly_map)
            (component_map, components) = self.construct_component_map(g_map, node_map)
        if self.options.scenario=='advanced':
            self.log('Graph has ' + str(node_map.max()) + ' nodes and '+ str(components.max())+ ' components.',2)
            
        if self.options.write_cur_maps == True:
            cum_current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')         
        
        if self.options.write_volt_maps == True: 
            cum_voltage_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
        elif self.options.scenario=='one-to-all':
            cum_voltage_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64') 
        
        if oneToAllStreamline==False:        
            cmin = 1
            cmax = components.max()
        else:
            cmin = componentWithPoints
            cmax = componentWithPoints
            
        for c in range(cmin,int(cmax+1)):
            c_map = where(component_map == c, 1, 0)
            local_source_map = multiply(c_map, source_map)
            local_ground_map = where(c_map, ground_map, 0) 
            del c_map
            
            source_in_component = (where(local_source_map, 1, 0)).sum() > 0
            ground_in_component = (where(local_ground_map, 1, 0)).sum() > 0
            
            if (source_in_component) & (ground_in_component):
                (rows, cols) = where(local_source_map)
                values = local_source_map[rows,cols]
                local_sources_rc = c_[values,rows,cols]
                (rows, cols) = where(local_ground_map)
                values = local_ground_map[rows,cols]
                local_grounds_rc = c_[values,rows,cols]
                del rows, cols, values, local_source_map, local_ground_map 

                if oneToAllStreamline==False:
                    (G, node_map) = self.node_pruner(g_map, poly_map, component_map, c)

                numnodes = node_map.max()
                sources = zeros(numnodes)
                grounds = zeros(numnodes)
                num_local_sources = local_sources_rc.shape[0]
                num_local_grounds = local_grounds_rc.shape[0]

                for source in range(0, num_local_sources):
                    src = self.grid_to_graph (local_sources_rc[source,1], local_sources_rc[source,2], node_map)
                    # Possible to have more than one source at a node when there are polygons
                    sources[src] = sources[src] + local_sources_rc[source,0] 

                for ground in range(0, num_local_grounds):
                    gnd = self.grid_to_graph (local_grounds_rc[ground,1], local_grounds_rc[ground,2], node_map)
                    # Possible to have more than one ground at a node when there are polygons
                    grounds[gnd] = grounds[gnd] + local_grounds_rc[ground,0] 

                (sources, grounds, finitegrounds) = self.resolve_conflicts(sources, grounds)

                solver_called = True
                try:
                    voltages = self.multiple_solver(G, sources, grounds, finitegrounds) 
                    del sources, grounds
                    if c==int(cmax):
                        del g_map, poly_map, ground_map
                        gc.collect()
                except MemoryError:
                    raise MemoryError
                except:
                    voltages = -777
                    solver_failed = True
                    
                if solver_failed==False:
                    ##Voltage and current mapping are cumulative, since there may be independent components.
                    if self.options.write_volt_maps == True: 
                        cum_voltage_map +=  self.create_voltage_map(node_map,voltages) 
                    elif self.options.scenario=='one-to-all':
                        cum_voltage_map +=  self.create_voltage_map(node_map,voltages) 
                    if self.options.write_cur_maps == True:
                        G = G.tocoo()
                        cum_current_map +=  self.create_current_map(voltages, G, node_map, finitegrounds) 
                
        if self.options.write_volt_maps == True: 
            if solver_failed==False:
                filetext = 'voltmap'
                if self.options.scenario=='advanced':
                    fileadd = ''
                else:
                    fileadd = '_'+str(source_id)
                self.write_aaigrid(filetext, fileadd, cum_voltage_map)  

        if self.options.scenario=='advanced':
            if self.options.write_cur_maps == True: 
                if solver_failed==False:
                    if self.options.log_transform_maps==True:
                        cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state['nodata']) 
                    filetext = 'curmap'
                    fileadd = ''
                    self.write_aaigrid(filetext, fileadd, cum_current_map) 
            else:
               cum_current_map = None

        else:
            if self.options.write_cur_maps == True: 
                if self.options.write_cum_cur_map_only==False:
                    if solver_failed==False:
                        if self.options.log_transform_maps==True:
                            cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state['nodata']) 
                        filetext = 'curmap'   
                        fileadd = '_'+str(source_id)   
                        self.write_aaigrid(filetext, fileadd, cum_current_map) 
            else:
                cum_current_map = None
            
        if self.options.scenario=='one-to-all':
            if solver_failed==False:
                (row, col) = where(source_map>0)
                voltages = cum_voltage_map[row,col]/source_map[row,col] #allows for variable source strength
        elif self.options.scenario=='all-to-one':
            if solver_failed==False:
                voltages = 0 #return 0 voltage/resistance for all-to-one mode           

        # Advanced mode will return voltages of the last component solved only for verification purposes.  
        if solver_called==False:
            voltages = -1

        return voltages,cum_current_map,solver_failed

        
    def resolve_conflicts(self, sources, grounds):
        """Handles conflicting grounds and sources for advanced mode according to user preferences."""  
        finitegrounds = where(grounds<Inf,grounds,0)
        if (where(finitegrounds==0, 0, 1)).sum()==0:
            finitegrounds = [-9999]
        infgrounds = where(grounds==Inf,1,0)
        
        ##Resolve conflicts bewteen sources and grounds
        conflicts = logical_and(sources,grounds)
        if self.options.remove_src_or_gnd=='rmvsrc':
            sources = where(conflicts,0,sources)
        elif self.options.remove_src_or_gnd=='rmvgnd':
            grounds = where(conflicts,0,grounds)
        elif self.options.remove_src_or_gnd=='rmvall':
            sources = where(conflicts,0,sources)
        infconflicts = logical_and(sources,infgrounds)
        grounds = where(infconflicts,0,grounds)
        if size(where(sources)) == 0:
            raise RuntimeError('All sources conflicted with grounds and were removed. There is nothing to solve.') 
        if size(where(grounds)) == 0:
            raise RuntimeError('All grounds conflicted with sources and were removed.  There is nothing to solve.') 

        return (sources, grounds, finitegrounds)

        
    @print_timing
    def multiple_solver(self, G, sources, grounds, finitegrounds):
        """Solver used for advanced mode."""  
        if finitegrounds[0]==-9999:#Fixme: no need to do this, right?
            finitegrounds = zeros(G.shape[0],dtype = 'int32') #create dummy vector for pairwise case
            Gsolve = G + sparse.spdiags(finitegrounds.T, 0, G.shape[0], G.shape[0]) 
            finitegrounds = [-9999]
        else:
            Gsolve = G + sparse.spdiags(finitegrounds.T, 0, G.shape[0], G.shape[0]) 
           
        ##remove infinite grounds from graph
        infgroundlist = where(grounds==Inf)
        infgroundlist = infgroundlist[0]
        numinfgrounds = infgroundlist.shape[0]
        
        dst_to_delete = []
        for ground in range(1, numinfgrounds+1):
            dst = infgroundlist[numinfgrounds-ground]
            dst_to_delete.append(dst)
            #Gsolve = deleterowcol(Gsolve, delrow = dst, delcol = dst)
            keep = delete (arange(0, sources.shape[0]), dst)
            sources = sources[keep]            
        Gsolve = deleterowcol(Gsolve, delrow = dst_to_delete, delcol = dst_to_delete)
        
        self.create_amg_hierarchy(Gsolve)
        voltages = self.solve_linear_system(Gsolve, sources)
        del Gsolve
        self.state['amg_hierarchy'] = None

        numinfgrounds = infgroundlist.shape[0]
        if numinfgrounds>0:
            #replace infinite grounds in voltage vector
            for ground in range(numinfgrounds,0, -1): 
                node = infgroundlist[numinfgrounds - ground] 
                voltages = asmatrix(insert(voltages,node,0)).T
        return asarray(voltages).reshape(voltages.size)
            
            
    @print_timing
    def construct_node_map(self, g_map, poly_map):
        """Creates a grid of node numbers corresponding to raster pixels with non-zero conductances."""  
        node_map = zeros(g_map.shape, dtype = 'int32')
        node_map[g_map.nonzero()] = arange(1, sum(g_map>0)+1, dtype = 'int32')

        if poly_map == []:
            return node_map

        # Remove extra points from poly_map that are not in g_map
        poly_map_pruned = zeros(g_map.shape, dtype = 'int32')
        poly_map_pruned[where(g_map)] = poly_map[where(g_map)]
        
        polynums = unique (poly_map)
   
        for i in range(0, polynums.size):
            polynum = polynums[i]
            if polynum !=  0:

                (pi, pj) = where (poly_map_pruned == polynum) #
                (pk, pl) = where (poly_map == polynum) #Added 040309 BHM                
                if len(pi)>0:  
                    node_map[pk, pl] = node_map[pi[0], pj[0]] #Modified 040309 BHM  
        node_map[where(node_map)] = relabel(node_map[where(node_map)], 1) #BHM 072409

        return node_map

    @print_timing
    def construct_component_map(self, g_map, node_map):
        """Assigns component numbers to grid corresponding to pixels with non-zero conductances.
        
        Nodes with the same component number are in single, connected components.
        
        """  
        prunedMap = False
        G = self.construct_g_graph(g_map, node_map, prunedMap) 
        (numComponents, C) = connected_components(G)
        C += 1 # Number components from 1

        (I, J) = where(node_map)
        nodes = node_map[I, J].flatten()

        component_map = zeros(node_map.shape, dtype = 'int32')
        component_map[I, J] = C[nodes-1]

        return (component_map, C)


    @print_timing
    def construct_g_graph(self, g_map, node_map,prunedMap):
        """Construct sparse adjacency matrix given raster maps of conductances and nodes."""  
        numnodes = node_map.max()
        (node1, node2, conductances) = self.get_conductances(g_map, node_map)
        G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes)) # Memory hogging operation?
        g_graph = G + G.T
        return g_graph

        
     
    def convert_graph_to_3_col(self,graph,nodeNames): 
        """Converts a sparse adjacency matrix to 3-column format."""  
        Gcoo =  graph.tocoo()
        mask = Gcoo.data > 0
        
        graphNcol = zeros((Gcoo.row[mask].size,3),dtype = "float64") #Fixme: this may result in zero-current nodes being left out.  Needed to make change so dimensions would match Gcoo.data[mask]
        
        if nodeNames==None:

            graphNcol[:,0] = Gcoo.row[mask]
            graphNcol[:,1] = Gcoo.col[mask]
        else:
            graphNcol[:,0]=nodeNames[Gcoo.row[mask]]
            graphNcol[:,1]=nodeNames[Gcoo.col[mask]]
        graphNcol[:,2] = Gcoo.data[mask]
        return graphNcol
        
        
    def get_horiz_neighbors(self, g_map):
        """Returns values of horizontal neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        g_map_l = g_map[:, 0:(n-1)]
        g_map_r = g_map[:, 1:n]
        g_map_lr = double(logical_and(g_map_l, g_map_r))
        s_horiz = where(c_[g_map_lr, zeros((m,1), dtype = 'int32')].flatten())
        t_horiz = where(c_[zeros((m,1), dtype = 'int32'), g_map_lr].flatten())

        return (s_horiz, t_horiz)

        
    def get_vert_neighbors(self, g_map):
        """Returns values of vertical neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        g_map_u = g_map[0:(m-1), :]
        g_map_d = g_map[1:m    , :]
        g_map_ud = double(logical_and(g_map_u, g_map_d))
        s_vert = where(r_[g_map_ud, zeros((1,n), dtype = 'int32')].flatten())
        t_vert = where(r_[zeros((1,n), dtype = 'int32'), g_map_ud].flatten())
        
        return (s_vert, t_vert)

        
    def get_diag1_neighbors(self, g_map):
        """Returns values of 1st diagonal neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        z1 = zeros((m-1, 1), dtype = 'int32')
        z2 = zeros((1  , n), dtype = 'int32')
        
        g_map_ul  = g_map[0:m-1, 0:n-1]
        g_map_dr  = g_map[1:m  , 1:n  ]
        g_map_udr = double(logical_and(g_map_ul, g_map_dr)) 
        s_dr      = where(r_[c_[g_map_udr, z1], z2].flatten())
        t_dr      = where(r_[z2, c_[z1, g_map_udr]].flatten())
        
        return (s_dr, t_dr)

        
    def get_diag2_neighbors(self, g_map):
        """Returns values of 2nd diagonal neighbors in conductance map."""  
        m = g_map.shape[0]
        n = g_map.shape[1]

        z1 = zeros((m-1, 1), dtype = 'int32')
        z2 = zeros((1  , n), dtype = 'int32')

        g_map_ur  = g_map[0:m-1, 1:n  ]
        g_map_dl  = g_map[1:m  , 0:n-1]
        g_map_udl = double(logical_and(g_map_ur, g_map_dl)) 
        s_dl      = where(r_[c_[z1, g_map_udl], z2].flatten())
        t_dl      = where(r_[z2, c_[g_map_udl, z1]].flatten())
                        
        return (s_dl, t_dl)
        
        
    def get_conductances(self, g_map, node_map):
        """Calculates conductances between adjacent nodes given a raster conductance map.
        
        Returns an adjacency matrix with values representing node-to-node conductance values.
        
        """  
        (s_horiz, t_horiz) = self.get_horiz_neighbors(g_map)
        (s_vert,  t_vert)  = self.get_vert_neighbors(g_map)

        s = c_[s_horiz, s_vert].flatten()
        t = c_[t_horiz, t_vert].flatten()

        # Conductances
        g1 = g_map.flatten()[s]
        g2 = g_map.flatten()[t]

        if self.options.connect_using_avg_resistances == False:
            conductances = (g1+g2)/2
        else:
            conductances = 1 /((1/g1+1/g2)/2)

        if self.options.connect_four_neighbors_only == False:
            (s_dr, t_dr) = self.get_diag1_neighbors(g_map)
            (s_dl, t_dl) = self.get_diag2_neighbors(g_map)

            sd = c_[s_dr, s_dl].flatten()
            td = c_[t_dr, t_dl].flatten()

            # Conductances
            g1 = g_map.flatten()[sd]
            g2 = g_map.flatten()[td]

            if self.options.connect_using_avg_resistances == False:
                conductances_d = (g1+g2) / (2*sqrt(2))
            else:
                conductances_d =  1 / (sqrt(2)*(1/g1 + 1/g2) / 2)

            conductances = r_[conductances, conductances_d]

            s = r_[s, sd].flatten()
            t = r_[t, td].flatten()

        # Nodes in the g_graph. 
        # Subtract 1 for Python's 0-based indexing. Node numbers start from 1
        node1 = node_map.flatten()[s]-1
        node2 = node_map.flatten()[t]-1
        
        return (node1, node2, conductances)

    @print_timing
    def node_pruner(self, g_map, poly_map, component_map, keep_component):
        """Removes nodes outside of component being operated on.
        
        Returns node map and adjacency matrix that only include nodes in keep_component.
        
        """  
        selector = component_map == keep_component
        
        g_map_pruned = selector * g_map
        poly_map_pruned = []
        if poly_map !=  []:
            poly_map_pruned = selector * poly_map

        node_map_pruned = self.construct_node_map (g_map_pruned, poly_map_pruned)
        prunedMap = True
        G_pruned = self.construct_g_graph (g_map_pruned, node_map_pruned, prunedMap) 
        G = self.laplacian(G_pruned) 
        
        return (G, node_map_pruned)

        
    @print_timing
    def laplacian(self, G): 
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
            self.state['amg_hierarchy'] = None
            # construct the MG hierarchy
            ml = []
            #  scipy.io.savemat('c:\\temp\\graph.mat',mdict={'d':G})
            ml = smoothed_aggregation_solver(G)
            self.state['amg_hierarchy'] = ml
  
        return

        
    @print_timing
    def solve_linear_system(self, G, rhs): 
        """Solves system of equations."""  
        gc.collect()
        # Solve G*x = rhs
        x = []
        if self.options.solver == 'cg+amg':
            ml = self.state['amg_hierarchy']
            G.psolve = ml.psolve
            (x, flag) = sparse.linalg.cg(G, rhs, tol = 1e-6, maxiter = 100000)
            if flag !=  0 or linalg.norm(G*x-rhs) > 1e-3:
                raise RuntimeError('CG did not converge. May need more iterations.') 

        if self.options.solver == 'amg':
            ml = self.state['amg_hierarchy']
            x = ml.solve(rhs, tol = 1e-6);

        return x 

        
    @print_timing
    def create_voltage_map(self, node_map, voltages):
        """Creates raster map of voltages given node voltage vector."""  
        voltage_map = numpy.zeros((self.state['nrows'], self.state['ncols']), dtype = 'float64')
        ind = node_map > 0
        voltage_map[where(ind)] = asarray(voltages[node_map[ind]-1]).flatten()
        return voltage_map


######################### BEGIN CURRENT MAPPING CODE ########################################
    @print_timing
    def create_current_map(self, voltages, G, node_map, finitegrounds):
        """Creates raster current map given node voltage vector, adjacency matrix, etc."""  
        gc.collect()
        node_currents = self.get_node_currents(voltages, G, finitegrounds)
        (rows, cols) = where(node_map)
        vals = node_map[rows, cols]-1
        current_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')
        current_map[rows,cols] = node_currents[vals]

        return current_map

    def get_node_currents(self, voltages, G, finitegrounds):
        """Calculates currents at nodes."""  
        node_currents_pos = self.get_node_currents_posneg (G, voltages, finitegrounds, True) 
        node_currents_neg = self.get_node_currents_posneg (G, voltages, finitegrounds, False)
        node_currents = where(node_currents_neg > node_currents_pos, node_currents_neg, node_currents_pos)

        return asarray(node_currents)[0]

        
    def get_node_currents_posneg (self, G, voltages, finitegrounds, pos):
        """Calculates positive or negative node currents based on pos flag."""  
        branch_currents = self.get_branch_currents(G,voltages,pos)
        branch_currents = branch_currents-branch_currents.T #Can cause memory error
        
        branch_currents = branch_currents.tocoo() #Can cause memory error, but this and code below more memory efficient than previous version.
        mask = branch_currents.data > 0
        row  = branch_currents.row[mask]
        col  = branch_currents.col[mask]
        data = branch_currents.data[mask]
        del mask
        n = G.shape[0]
        branch_currents = sparse.csr_matrix((data, (row, col)), shape = (n,n))
           
        if finitegrounds[0]!= -9999:  
            finiteground_currents = multiply(finitegrounds, voltages)
            if pos == True:
                finiteground_currents = where(finiteground_currents < 0, -finiteground_currents, 0)
            else:
                finiteground_currents = where(finiteground_currents > 0, finiteground_currents, 0)  
            n = G.shape[0]
            branch_currents = branch_currents + sparse.spdiags(finiteground_currents.T, 0, n, n)        

        return branch_currents.sum(0)
    
    
    def get_branch_currents(self,G,voltages,pos):    
        """Calculates branch currents."""  
        branch_currents = self.get_branch_currents_posneg(G,voltages,pos)
        n = G.shape[0]
        mask = G.row < G.col
        branch_currents = sparse.csr_matrix((branch_currents, (G.row[mask], G.col[mask])), shape = (n,n)) #SQUARE MATRIX, SAME DIMENSIONS AS GRAPH
        return branch_currents

        
    def get_branch_currents_posneg(self,G,voltages,pos):
        """Calculates positive or negative node currents based on pos flag."""  
        mask = G.row < G.col
        if pos==True:
             vdiff = voltages[G.row[mask]]              
             vdiff -=  voltages[G.col[mask]]             

        else:
             vdiff = voltages[G.col[mask]]              
             vdiff -=  voltages[G.row[mask]]             

        conductances = where(G.data[mask] < 0, -G.data[mask], 0)
        del mask
        
        branch_currents = asarray(multiply(conductances,vdiff.T)).flatten()
        maxcur = max(branch_currents)
        branch_currents = where(absolute(branch_currents/maxcur) < 1e-8,0,branch_currents) #Delete very small branch currents to save memory
        return branch_currents
######################### END CURRENT MAPPING CODE ########################################        
        
        
    ### FILE I/O ###
    def write_aaigrid(self, type, fileadd, data):
        """Writes ASCII grid.  This is main raster output format for Circuitscape."""  
        if type == 'voltmap':
            if self.options.write_volt_maps == False: 
                return
            if self.options.set_null_voltages_to_nodata==True:
                ind = self.state['g_map'] == 0
                data[where(ind)] = self.state['nodata']
                del ind                
        elif type == 'curmap' or type == 'cum_curmap' or type == 'max_curmap':
            if self.options.write_cur_maps == False: 
                return
            if self.options.set_null_currents_to_nodata == True:
                ind = self.state['g_map'] == 0
                data[where(ind)] = self.state['nodata']
                del ind                
        else:
            return

        outputDir, outputFile = os.path.split(self.options.output_file)
        outputBase, outputExtension = os.path.splitext(outputFile)

        inputBase, inputExtension = os.path.splitext(self.options.habitat_file) 
        if inputExtension == '.npy': #if read in numpy array, write one out.
            file = outputDir + '//' + outputBase + '_' + type + fileadd +'.npy'
        else:
            file = outputDir + '//' + outputBase + '_' + type + fileadd +'.asc'
        cs_io.writer(file, data, self.state, self.options.compress_grids)
        
        
    def read_cell_map(self, filename):
        """Reads resistance or conductance raster into memory, converts former to conductance format."""  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = cs_io.read_header(filename)
        self.state['ncols'] = ncols
        self.state['nrows'] = nrows
        self.state['xllcorner'] = xllcorner
        self.state['yllcorner'] = yllcorner
        self.state['cellsize'] = cellsize
        if nodata==False:
            self.state['nodata'] = -9999
        else:
            self.state['nodata'] = nodata

        cell_map = cs_io.reader(filename, 'float64')

        ######################## Reclassification code
        if self.options.use_reclass_table==True:
            try:
                reclassTable = self.readPointStrengths(self.options.reclass_file)    
            except:
                raise RuntimeError('Error reading reclass table')
            for i in range (0,reclassTable.shape[0]):
                cell_map = where(cell_map==reclassTable[i,0],reclassTable[i,1],cell_map)
            print'\n***** Reclassified habitat map using', self.options.reclass_file,'*****'
        ########################
        
        if self.options.habitat_map_is_resistances == True:
            zeros_in_resistance_map = (where(cell_map==0, 1, 0)).sum() > 0
            if zeros_in_resistance_map == True: #FIXME: Should be easy to accomodate zeros in resistance map, just treat them like polygons.
                raise RuntimeError('Error: zero resistance values are not currently supported for habitat maps.  Use a short-circuit region file instead.')
            g_map = 1 / cell_map  
            g_map = where(cell_map == -9999,0,g_map)
        else:
            g_map = where(cell_map == -9999,0,cell_map)    
        g_map = where(g_map < 0,0,g_map)    
        return g_map


    def read_point_map(self, filename):
        """Reads map or text list of focal nodes from disk.  
        
        File extension is used to determine whether format is ascii grid, numpy array, or text list.
        
        """  
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        base, extension = os.path.splitext(filename)
        
        if extension == ".txt":
            try:
                points = loadtxt(filename)
            except ValueError:
                raise RuntimeError('File "'  + filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                
            
            points_rc = zeros(points.shape,dtype = 'int32')
            try:
                points_rc[:,0] = points[:,0]
                points_rc[:,1] = ceil((self.state['nrows']-(points[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
                points_rc[:,2] = ceil(((points[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
                i = argsort(points_rc[:,0])
                points_rc = points_rc[i]
            except IndexError:
                raise RuntimeError('Error extracting focal node locations. Please check file format.')                

        elif extension == ".asc" or extension == ".npy": # We use Numpy format for quickly passing grids between ArcGIS and Circuitscape.
            readingMask = False
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = cs_io.read_header(filename)
            if cellsize!= self.state['cellsize']:
                print'\n********\nWarning: Focal node raster has different \ncell size than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif ncols!= self.state['ncols']:
                print'\n********\nWarning: Focal node raster has different \nnumber of columns than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif nrows!= self.state['nrows']:
                print'\n********\nWarning: Focal node raster has different \nnumber of rows than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif xllcorner!= self.state['xllcorner']:
                print'\n********\nWarning: Focal node raster has different \nxllcorner than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            elif yllcorner!= self.state['yllcorner']:
                print'\n********\nWarning: Focal node raster has different \nyllcorner than habitat raster. \nCircuitscape will try to crudely resample the focal node raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
                point_map = self.resampleMap(filename,readingMask)
            else:
                point_map = cs_io.reader(filename, 'int32')

            (rows, cols) = where(point_map > 0)

            values = zeros(rows.shape,dtype = 'int32') 
            for i in range(0, rows.size):
                values[i] = point_map[rows[i], cols[i]]
            points_rc = c_[values,rows,cols]
            try:            
                i = argsort(points_rc[:,0])
                points_rc = points_rc[i]
            except IndexError:
                raise RuntimeError('Error extracting focal node locations. Please check file format.')                
        else:
            raise RuntimeError('Focal node file must have a .txt or .asc extension')
        
        #Check to make sure points fall within cellmap
        if min(points_rc[:,1])<0 or min(points_rc[:,2])<0:
            raise RuntimeError('At least one focal node location falls outside of habitat map')
        elif max(points_rc[:,1])>self.state['nrows']-1:
            raise RuntimeError('At least one focal node location falls outside of habitat map')
        elif  max(points_rc[:,2])>self.state['ncols']-1:
            raise RuntimeError('At least one focal node location falls outside of habitat map')
        if (unique(asarray(points_rc[:,0]))).shape[0]<2:
            raise RuntimeError('Less than two valid focal nodes found. Please check focal node location file.')                    
        
        return points_rc
        
        
    def read_poly_map(self, filename,readingMask):  
        """Reads short-circuit region map (aka polygon map) from disk."""  
        (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = cs_io.read_header(filename)
        if cellsize!= self.state['cellsize']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \ncell size than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif ncols!= self.state['ncols']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nnumber of columns than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif nrows!= self.state['nrows']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nnumber of rows than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif xllcorner!= self.state['xllcorner']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nxllcorner than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)
        elif yllcorner!= self.state['yllcorner']:
            print'\n********\nWarning: Short-circuit region or mask raster has different \nyllcorner than habitat raster. \nCircuitscape will try to crudely resample the raster. \nWe recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.'
            map = self.resampleMap(filename,readingMask)            
        else:
            map = cs_io.reader(filename, 'int32')
            
        map = where(map == nodata,0,map)        

        if readingMask==True:
            map = where(map < 0, 0, map)        

        return map

        
    def read_source_and_ground_maps(self, source_filename, ground_filename): 
        """Reads srouce and ground raster maps from disk."""  
        #FIXME: reader does not currently handle infinite inputs for ground conductances.
        if os.path.isfile(source_filename)==False:
            raise RuntimeError('File "'  + source_filename + '" does not exist')
        base, extension = os.path.splitext(source_filename)
        if extension == ".txt":  #FIXME: probably want to roll code used for reading source, ground and point text files into single utility
            try:
                sources = loadtxt(source_filename)
            except ValueError:
                raise RuntimeError('File "'  + source_filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                

            sources_rc = zeros(sources.shape,dtype = 'int32')
            sources_rc[:,0] = sources[:,0]
            sources_rc[:,1] = ceil((self.state['nrows']-(sources[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
            sources_rc[:,2] = ceil(((sources[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
            source_map = zeros((self.state['nrows'],self.state['ncols']),dtype = 'float64')
            source_map[sources_rc[:,1],sources_rc[:,2]] = sources_rc[:,0]
        elif extension=='.asc':
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = cs_io.read_header(source_filename)
            if cellsize!= self.state['cellsize']:
                raise RuntimeError('Current source raster must have same cell size and number of rows and columns as habitat raster') 
            if ncols!= self.state['ncols']:
                raise RuntimeError('Current source raster must have same cell size and number of rows and columns as habitat raster') 
            if nrows!= self.state['nrows']:
                raise RuntimeError('Current source raster must have same cell size and number of rows and columns as habitat raster') 
            if yllcorner!= self.state['yllcorner']:
                raise RuntimeError('Current source raster must have same xllcorner and yllcorner as habitat raster') 
            if xllcorner!= self.state['xllcorner']:
                raise RuntimeError('Current source raster must have same xllcorner and yllcorner as habitat raster') 
               
            source_map = cs_io.reader(source_filename, 'float64')
            source_map = where(source_map == -9999,0,source_map)

        else:
            raise RuntimeError('Current source files must have a .txt or .asc extension')
        if self.options.use_unit_currents==True:
            source_map = where(source_map,1,0)

        if os.path.isfile(ground_filename)==False:
            raise RuntimeError('File "'  + ground_filename + '" does not exist')
        base, extension = os.path.splitext(ground_filename)
        if extension == ".txt":
            try:
                grounds = loadtxt(ground_filename)
            except ValueError:
                raise RuntimeError('File "'  + ground_filename + '" is not in correct text list format. \n If it is an ASCII grid, please use .asc extension.')                

            grounds_rc = zeros(grounds.shape,dtype = 'int32')
            grounds_rc[:,0] = grounds[:,0]
            grounds_rc[:,1] = ceil((self.state['nrows']-(grounds[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
            grounds_rc[:,2] = ceil(((grounds[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
            ground_map_raw = -9999*ones((self.state['nrows'],self.state['ncols']),dtype = 'float64')
            ground_map_raw[grounds_rc[:,1],grounds_rc[:,2]] = grounds_rc[:,0]
        elif extension=='.asc':
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = cs_io.read_header(ground_filename)
            if cellsize!= self.state['cellsize']:
                raise RuntimeError('Ground raster must have same cell size and number of rows and columns as habitat raster') 
            if ncols!= self.state['ncols']:
                raise RuntimeError('Ground raster must have same cell size and number of rows and columns as habitat raster') 
            if nrows!= self.state['nrows']:
                raise RuntimeError('Ground raster must have same cell size and number of rows and columns as habitat raster') 
            if yllcorner!= self.state['yllcorner']:
                raise RuntimeError('Ground raster must have same xllcorner and yllcorner as habitat raster') 
            if xllcorner!= self.state['xllcorner']:
                raise RuntimeError('Ground raster must have same xllcorner and yllcorner as habitat raster') 
                
            ground_map_raw = cs_io.reader(ground_filename, 'float64')
        else:
            raise RuntimeError('Ground files must have a .txt or .asc extension')
        if self.options.ground_file_is_resistances==True:
            ground_map = 1 / ground_map_raw
            ground_map = where(ground_map_raw == -9999,0,ground_map)
        else:
            ground_map = where(ground_map_raw == -9999,0,ground_map_raw)
        if self.options.use_direct_grounds==True:
            ground_map = where(ground_map,Inf,0)

        conflicts = logical_and(source_map,ground_map)
        if self.options.remove_src_or_gnd=='rmvsrc':
            source_map = where(conflicts,0,source_map)
        elif self.options.remove_src_or_gnd=='rmvgnd':
            ground_map = where(conflicts,0,ground_map)
        elif self.options.remove_src_or_gnd=='rmvall':
            source_map = where(conflicts,0,source_map)
            ground_map = where(conflicts,0,ground_map)
        if size(where(source_map)) == 0:
            raise RuntimeError('No valid sources detected. Please check source file') 
        if size(where(ground_map)) == 0:
            raise RuntimeError('No valid grounds detected. Please check ground file') 
        return source_map, ground_map


    def readincludedPairs(self, filename):
        """Reads matrix denoting node pairs to include/exclude from calculations.
        
        FIXME: matrices are an inconvenient way for users to specify pairs.  Using a 
        2- or 3-column format would be easier.
        
        """  
        
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        
        try:
            f = open(filename, 'r')
            [ign, minval] = string.split(f.readline())
            [ign, maxval] = string.split(f.readline())
            minval = float(minval)
            maxval = float(maxval)
            f.close()            
            
            includedPairs = loadtxt(filename, skiprows = 2, dtype = 'Float64')
            pointIds = includedPairs[:,0]
            includedPairs = where(includedPairs>maxval,0,includedPairs)
            includedPairs = where(includedPairs<minval,0,1)             
            includedPairs[:,0] = pointIds
            includedPairs[0,:] = pointIds
            includedPairs[0,0] = -1
            i = argsort(includedPairs[:,0])
            includedPairs = includedPairs[i]
            includedPairs = includedPairs.T            
            i = argsort(includedPairs[:,0])
            includedPairs = includedPairs[i] 
        
        except:
            raise RuntimeError('Error reading focal node include/exclude matrix. Please check file format.')                

        return includedPairs

        
    def readPointStrengths(self, filename):
        """Reads list of variable source strengths from disk.  
        
        This code also used for reading file for reclassifying input data.
        
        """  
        if os.path.isfile(filename)==False:
            raise RuntimeError('File "'  + filename + '" does not exist')
        
        try:
            pointStrengths = loadtxt(filename)
        except ValueError:
            raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
           
        try:
            pointIds = pointStrengths[:,0]
            i = argsort(pointStrengths[:,0])
            pointStrengths = pointStrengths[i]

        except:
            raise RuntimeError('Error reading focal node source strength list. Please check file format.')                
        
        return pointStrengths


    def enable_low_memory(self, restart):
        """Runs circuitscape in low memory mode.  Not incredibly helpful it seems."""  
        self.state['amg_hierarchy'] = None
        gc.collect()
        if self.options.low_memory_mode==True:
            if restart==False: #If this module has already been called
                raise MemoryError
        self.options.low_memory_mode = True
        print'\n**************\nMemory error reported.'        

        type, value, tb = sys.exc_info()
        info = traceback.extract_tb(tb)
        print'Full traceback:'
        print info
        print'***************'
        filename, lineno, function, text = info[-1] # last line only
        print"\n %s:%d: %s: %s (in %s)" %\
              (filename, lineno, type.__name__, str(value), function)

        type = value = tb = None # clean up
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


    def resampleMap(self,filename,readingMask):
        """Code to crudely resample input raster if raster headers don't match (i.e. different extents or cell sizes used)."""  
        try:
            (ncols, nrows, xllcorner, yllcorner, cellsize, nodata) = cs_io.read_header(filename)
            map = cs_io.reader(filename, 'int32')
            map = where(map == nodata,0,map)
            
            if readingMask==True:
                map = where(map>0, 1, 0)      
                map = 1-map #now zeros are areas to keep, ones are masked out
                
            (rows, cols) = where(map > 0)
            
            xcoords = (cols+0.5)*cellsize+xllcorner
            ycoords = (nrows-rows-0.5)*cellsize+yllcorner
            
            values = zeros(rows.shape,dtype = 'int32') 
            for i in range(0, rows.size):
                values[i] = map[rows[i], cols[i]]
            mapCoords = c_[values,xcoords,ycoords]
    
            i = argsort(mapCoords[:,0])
            mapCoords = mapCoords[i]
            
            #From mapCoords to mapRc:
            mapRc = zeros(mapCoords.shape,dtype = 'int32')
            mapRc[:,0] = mapCoords[:,0]
            mapRc[:,1] = ceil((self.state['nrows']-(mapCoords[:,2]-self.state['yllcorner'])/self.state['cellsize']))-1
            mapRc[:,2] = ceil(((mapCoords[:,1]-self.state['xllcorner'])/self.state['cellsize']))-1
            i = argsort(mapRc[:,0])
            mapRc = mapRc[i]
    
            rows = mapRc[:,1]
            (delrows) = asarray(where(rows<0))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = deleterow(mapRc,delrows2)
            rows = mapRc[:,1] 
            (delrows) = asarray(where(rows>self.state['nrows']-1))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = deleterow(mapRc,delrows2)
            cols = mapRc[:,2]
            (delrows) = asarray(where(cols<0))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = deleterow(mapRc,delrows2)
            cols = mapRc[:,2]
            (delrows) = asarray(where(cols>self.state['ncols']-1))
            delrows2 = zeros(delrows.shape[1]) #turn into 1-d array
            delrows2[:] = delrows[:]
            if delrows2!= []:
                mapRc = deleterow(mapRc,delrows2)
            del delrows
            del delrows2
    
            #From mapRc to map:                
            map = numpy.zeros((self.state['nrows'],self.state['ncols']),int)
            map[mapRc[:,1],mapRc[:,2]] = mapRc[:,0] 
            if readingMask==True:
                map = 1-map #now zeros are areas to mask out, ones are kept
            
        except: 
            raise RuntimeError('Error resampling focal node, mask, or short-circuit region locations to match habitat map cell size and extent.  We recommend using the "Export to Circuitscape" ArcGIS tool to create ASCII grids with compatible cell size and extent.')   

        return map

        
    def getVoltmatrix(self,i,j,numpoints,local_node_map,voltages,points_rc,resistances,voltmatrix):                                            
        """Returns a matrix of pairwise voltage differences between focal nodes.
        
        Used for shortcut calculations of effective resistance when no
        voltages or currents are mapped.
        
        """  
        voltvector = zeros((numpoints,1),dtype = 'float64')  
        voltage_map = self.create_voltage_map(local_node_map,voltages) 
        for point in range(1,numpoints):
            voltageAtPoint = voltage_map[points_rc[point,1], points_rc[point,2]]
            voltageAtPoint = 1-(voltageAtPoint/resistances[i, j])
            voltvector[point] = voltageAtPoint
        voltmatrix[:,j] = voltvector[:,0] 
        return voltmatrix


    def getShortcutResistances(self,anchorPoint,voltmatrix,numpoints,resistances,shortcutResistances): #FIXME: no solver failed capability
        """Calculates all resistances to each focal node at once.
        
        Greatly speeds up resistance calculations if not mapping currents or voltages.
        
        """  
        point1 = anchorPoint
        for pointx in range(0, numpoints): #point1 is source node, i.e. the 1 in R12.  point 2 is the dst node.
            R1x = resistances[point1,pointx]
            if R1x!= -1:
                shortcutResistances[point1,pointx] = R1x
                shortcutResistances[pointx,point1] = R1x
                for point2 in range(pointx,numpoints):
                    R12 = resistances[point1,point2] 
                    if R12!= -1:
                        shortcutResistances[point2,point1] = R12
                        shortcutResistances[point1,point2] = R12
                        Vx = voltmatrix[pointx,point2]
                        R2x = 2*R12*Vx+R1x-R12
                        shortcutResistances[pointx,point2] = R2x
                        shortcutResistances[point2,pointx] = R2x   
        return shortcutResistances                        


    def writeResistances(self, point_ids, resistances):
        """Writes resistance file to disk."""  
        outputResistances = self.append_names_to_resistances(point_ids, resistances)
        fileName = self.options.output_file
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        outputFile = outputDir + '//' + outputBase + '_resistances' + outputExtension 
        savetxt (outputFile, outputResistances)

        resistances3Columns = self.convertResistances3cols(outputResistances)
        
        self.writeResistances3columns(resistances3Columns)
        
        #remove partial result file        
        oldFile = outputDir + '//' + outputBase + '_resistances_incomplete' + outputExtension
        try:
            os.remove(oldFile)
        except:
            pass 
        return outputResistances
        
        
    def writeResistances3columns(self, resistances3Columns):    
        """Writes effective resistances to disk in 3 column format."""  
        fileName = self.options.output_file
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)       
        outputFile = outputDir + '//' + outputBase + '_resistances_3columns' + outputExtension 
        savetxt (outputFile, resistances3Columns)
        return         
        
        
    def convertResistances3cols(self, resistances):
        """Converts resistances from matrix to 3-column format."""  
        numPoints = resistances.shape[0]-1
        numEntries = numPoints*(numPoints-1)/2
        resistances3columns = zeros((numEntries,3),dtype = 'float64') 
        x = 0
        for i in range(1,numPoints):
            for j in range(i+1,numPoints+1):
                resistances3columns[x,0] = resistances[i,0]    
                resistances3columns[x,1] = resistances[0,j]
                resistances3columns[x,2] = resistances[i,j]
                x = x+1
        return resistances3columns        
        
        
    def saveIncompleteResistances(self, resistances):
        """Saves resistances from ongoing calculations.  
        
        Helpful for debugging or recovering partial results after crash or user abort.
        
        """  
        fileName = self.options.output_file
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        outputFile = outputDir + '//' + outputBase + '_resistances_incomplete' + outputExtension
        savetxt (outputFile, resistances)
        return


    def pruneIncludedPairsNetwork(self,focalNodes):
        """Remove excluded points from focal node list when using extra file that lists pairs to include/exclude in network mode."""   
        includedPairs = (self.state['includedPairs'])
        includeList = list(includedPairs[0,:])
        point = 0
        dropFlag = False
        while point < focalNodes.size: #Prune out any points not in includeList
            if focalNodes[point] in includeList: #match
                point = point+1
            else:
                dropFlag = True   
                focalNodes=delete(focalNodes,point)
         
        includeList = list(focalNodes[:])
        numConnectionRows = includedPairs.shape[0]
        row = 1
        while row <numConnectionRows: #Prune out any entries in includeList that are not in focalNodes
            if includedPairs [row,0] in includeList: #match
                row = row+1
            else:
                includedPairs = deleterowcol(includedPairs,delrow = row,delcol = row)   
                dropFlag = True
                numConnectionRows = numConnectionRows-1

        self.state['includedPairs'] = includedPairs                     
#         if dropFlag==True:
#             print'\nNOTE: Code to exclude pairwise calculations is activated and \nsome entries did not match with focal node file.  \nSome focal nodes may have been dropped.'      
        return focalNodes


    def pruneIncludedPairs(self,points_rc):
        """Remove excluded points from focal node list when using extra file that lists pairs to include/exclude."""  
        includedPairs = (self.state['includedPairs'])
        includeList = list(includedPairs[0,:])
        point = 0
        dropFlag = False
        while point < points_rc.shape[0]: #Prune out any points not in includeList
            if points_rc[point,0] in includeList: #match
                point = point+1
            else:
                dropFlag = True   
                points_rc = deleterow(points_rc,point)  
         
        includeList = list(points_rc[:,0])
        numConnectionRows = includedPairs.shape[0]
        row = 1
        while row <numConnectionRows: #Prune out any entries in includeList that are not in points_rc
            if includedPairs [row,0] in includeList: #match
                row = row+1
            else:
                includedPairs = deleterowcol(includedPairs,delrow = row,delcol = row)   
                dropFlag = True
                numConnectionRows = numConnectionRows-1

        self.state['includedPairs'] = includedPairs                     
#         if dropFlag==True:
#             print'\nNOTE: Code to exclude pairwise calculations is activated and \nsome entries did not match with focal node file.  \nSome focal nodes may have been dropped.'      
        return points_rc
    
    
    def get_points_rc_unique(self,point_ids,points_rc):
        """Return a list of unique focal node IDs and x-y coordinates."""  
        points_rc_unique = zeros((point_ids.size,3), int)
        for i in range(0, point_ids.size):
            for j in range(0, points_rc.shape[0]):
                if points_rc[j,0]==point_ids[i]:
                    points_rc_unique[i,:] = points_rc[j,:] 
                    break                    
        return points_rc_unique          
        
        
    def checkPointsInComponent(self,c,numpoints,components,points_rc,node_map):
        """Checks to see if there are focal points in a given component."""  
        points_in_this_component = False            
        for pt1 in range(0, numpoints): 
            if points_in_this_component == False:
                src = self.grid_to_graph (points_rc[pt1,1], points_rc[pt1,2], node_map)
                for pt2 in range(pt1+1, numpoints):
                    dst = self.grid_to_graph (points_rc[pt2,1], points_rc[pt2,2], node_map)
                    if (src >=  0 and components[src] == c) and (dst >=  0 and components[dst] == c):
                        points_in_this_component = True
                        break        
        return points_in_this_component    

        
    def getstrengthMap(self,points_rc_unique,pointStrengths):
        """Returns map and coordinates of point strengths when variable source strengths are used."""  
        if self.options.use_variable_source_strengths==True:
            if self.options.scenario == 'one-to-all': 
                strengths_rc = self.get_strengths_rc(self.state['pointStrengths'],points_rc_unique)
                strengthMap = None
            else:
                strengths_rc = self.get_strengths_rc(pointStrengths,points_rc_unique)
                strengthMap = numpy.zeros((self.state['nrows'],self.state['ncols']),dtype = 'Float64')
                strengthMap[points_rc_unique[:,1],points_rc_unique[:,2]] = strengths_rc[:,0]     
            return strengthMap,strengths_rc
        else:
            return None,None
        
        
    def get_strengths_rc(self,pointStrengths,points_rc_unique):
        """Returns coordinates of point strengths when variable source strengths are used."""  
        strengths_rc = zeros(points_rc_unique.shape,dtype = 'float64')
        strengths_rc[:,1] = points_rc_unique[:,1]
        strengths_rc[:,2] = points_rc_unique[:,2]
        for point in range(0,points_rc_unique.shape[0]):
            try:
                strengthIds = list(pointStrengths[:,0])
                strengthValues = list(pointStrengths[:,1])
                pointId = points_rc_unique[point,0]
                indx = strengthIds.index(pointId)
                strengths_rc[point,0] = strengthValues[indx]
            except ValueError:
                strengths_rc[point,0] = 1
        return strengths_rc        

        
    @print_timing
    def load_maps(self):
        """Loads all raster maps into self.state."""  
        self.log('Reading maps',1)
        self.log('',2)
        self.state['g_map'] = self.read_cell_map(self.options.habitat_file)
        if self.options.use_polygons:
            self.state['poly_map'] = self.read_poly_map(self.options.polygon_file,readingMask = False)
        else:
            self.state['poly_map'] = []
 
        if self.options.use_mask==True:
            mask = self.read_poly_map(self.options.mask_file,readingMask = True)
            mask = where(mask !=  0, 1, 0) 
            self.state['g_map'] = multiply(self.state['g_map'],mask)
            
            sumGmap = (self.state['g_map']).sum()
            sumGmap = sumGmap.sum()
            if sumGmap==0:
                raise RuntimeError('All entries in habitat map have been dropped after masking with the mask file.  There is nothing to solve.')             
            del mask
        else:
            self.state['mask'] = []
            

        if self.options.scenario=='advanced':
            self.state['points_rc'] = []
            (self.state['source_map'], self.state['ground_map']) = self.read_source_and_ground_maps(self.options.source_file, self.options.ground_file)

        else:        
            self.state['points_rc'] = self.read_point_map(self.options.point_file)
            self.state['source_map'] = []
            self.state['ground_map'] = []

        if self.options.use_included_pairs==True:
            self.state['includedPairs'] = self.readincludedPairs(self.options.included_pairs_file)
        
        self.state['pointStrengths'] = None
        if self.options.use_variable_source_strengths==True:
            self.state['pointStrengths'] = self.readPointStrengths(self.options.variable_source_file) 
        
        self.log('Processing maps',1)
        return 
        
        
    def writeResistancesOneToAll(self,resistances,string):
        """Saves effective resistances from one-to-all calculations to disk."""  
        fileName = self.options.output_file
        outputDir, outputFile = os.path.split(fileName)
        outputBase, outputExtension = os.path.splitext(outputFile)
        if self.options.scenario == 'one-to-all':
            outputFile = outputDir + '//' + outputBase + '_resistances' + string + outputExtension
        else:
            outputFile = outputDir + '//' + outputBase + '_results' + string + outputExtension     
        savetxt (outputFile, resistances) 
        
        #remove partial result file        
        if string=='':        
            if self.options.scenario == 'one-to-all':
                oldFile = outputDir + '//' + outputBase + '_resistances_incomplete' + string + outputExtension
            else:
                oldFile = outputDir + '//' + outputBase + '_results_incomplete' + string + outputExtension
            try:
                os.remove(oldFile)
            except:
                pass 
        return

    # Not implemented at this time
    def advanced_module_network(self,G,sources,grounds,nodeNames):
        """Overhead module for advanced mode with arbitrary graphs.
        
        Represents functionality we'll want to have in 4.0.
        
        NOT IMPLEMENTED YET
        
        """  
        (sources, grounds, finitegrounds) = self.resolve_conflicts(sources, grounds)
        try:
            voltages = self.multiple_solver(G, sources, grounds, finitegrounds)
            #fixme: need to deal with case with no valid sources or grounds in a component.
            
            solver_failed = False
        except:
            solver_failed = True
        if self.options.write_cur_maps == True:
            (node_currents,branch_currents) = self.getCurrentsNetwork(G,voltages,finitegrounds)
            fileadd=''
            outputNodeCurrents = self.writeCurrentsNetwork(branch_currents, node_currents, nodeNames, fileadd)
        if self.options.write_volt_maps == True:
            self.writeVoltagesNetwork(voltages, nodeNames, fileadd='')
        
        return voltages, solver_failed
        
        
    # Not implemented at this time        
    def one_to_all_module_network(self,G,focalNodes,nodeNames):
        """Not implemented yet.  Represents functionality we'll want to have in 4.0"""  
        solver_failed = False
        if self.options.use_included_pairs==True: #Prune points
            focalNodes = self.pruneIncludedPairsNetwork(focalNodes)          
            includedPairs = self.state['includedPairs'] 

        numpoints = focalNodes.size
        numnodes = G.shape[0]
        resistances = zeros((focalNodes.size,2),dtype = 'float64')
        resistances[:,0] = focalNodes
        if self.options.write_cur_maps == True:
            cumNodeCurrents = zeros((nodeNames.size,1),dtype = 'float64')
            cumBranchCurrents = sparse.csr_matrix((G.shape))
        
        finitegrounds = [-9999]
        focalNodeLocs = self.namesToNodes(nodeNames,focalNodes) 
        
        sources = zeros(numnodes,dtype = 'float64')
        grounds = zeros(numnodes,dtype = 'float64')
        
        variableSources = self.getVariableSources(numnodes,focalNodes)
        
        if self.options.scenario == 'one-to-all':
            grounds[focalNodeLocs] = Inf
        else:
            sources = variableSources

        x = 0
        lastWriteTime = time.time()
        for i in range(0, numpoints):
            if self.options.use_included_pairs==True: #Prune points
                groundsTemp=grounds
                sourcesTemp=sources
                for pair in range(0, focalNodes.size): #loop thru exclude[point,:], delete included pairs of focal point from point_map and points_rc_unique_temp 
                    if includedPairs[i+1,pair+1]==0 and i !=  pair:
                        dropNode = focalNodeLocs[pair]
                        groundsTemp[dropNode] = 0
                        sourcesTemp[dropNode] = 0
                        
            x = x+1
            (hours,mins,secs) = elapsed_time(self.state['startTime'])
            self.log ('At ' + str(hours) +' hr ' + str(mins) + ' min solving focal node ' + str(x) + ' of '+ str(numpoints) + '.',1)

            node = focalNodeLocs[i]
            if self.options.scenario == 'one-to-all':
                grounds[node] = 0
                sources[node] = variableSources[node]
                try:
                    voltages = self.multiple_solver(G, sources, grounds, finitegrounds)
                    resistances[i,1] = voltages[node]
                except:
                    solver_failed = True
                    resistances[i,1] = -777

                grounds[node] = Inf
                sources[node] = 0
            else:
                grounds[node] = Inf
                sources[node] = 0
                try:
                    voltages = self.multiple_solver(G, sources, grounds, finitegrounds)
                except:
                    solver_failed = True
                    resistances[i,1] = -777
                grounds[node] = 0
                sources[node] = variableSources[node]               

            (hours,mins,secs) = elapsed_time(lastWriteTime)
            if secs > 120: 
                lastWriteTime = time.time()
                self.writeResistancesOneToAll(resistances,'_incomplete')       
            
            if self.options.write_cur_maps == True:
                (node_currents,branch_currents) = self.getCurrentsNetwork(G,voltages,finitegrounds)
                cumNodeCurrents=cumNodeCurrents+node_currents
                cumBranchCurrents=cumBranchCurrents+branch_currents
                if self.options.write_cum_cur_map_only==False:     
                    outputNodeCurrents = self.writeCurrentsNetwork(branch_currents, node_currents, nodeNames, fileadd=str(int(focalNodes[i])))
            if self.options.write_volt_maps == True:
                self.writeVoltagesNetwork(voltages, nodeNames, fileadd=str(int(focalNodes[i])))
                
        if self.options.write_cur_maps == True:
            outputNodeCurrents = self.writeCurrentsNetwork(cumBranchCurrents, cumNodeCurrents, nodeNames, 'cum')
            
        self.writeResistancesOneToAll(resistances,'')
        
        #Need to add row and column headers and write currents and voltagesto disk

        print 'focal nodes'
        print focalNodes
        return resistances,solver_failed

    # Not implemented at this time   
    def graph_list_to_graph(self,graphList):
        """Converts 3-column adjacency list to sparse adjacency matrix.
        
        NOT IMPLEMENTED CURRENTLY
        
        """  
        nodes = deletecol(graphList,2) 
        nodeNames = unique(asarray(nodes))
        nodes[where(nodes>= 0)] = relabel(nodes[where(nodes>= 0)], 0)
        node1 = nodes[:,0]
        node2 = nodes[:,1]
        data = graphList[:,2] # Edge weights

        numnodes = nodeNames.shape[0]
        G = sparse.csr_matrix((data, (node1, node2)), shape = (numnodes, numnodes))       

        Gdense=G.todense()
        graph = maximum(Gdense, Gdense.T) # To handle single or double entries for elements BHM 06/28/11
        graph = sparse.csr_matrix(graph)

        return graph, nodeNames


    # Not implemented at this time   
    def getVariableSources(self,numnodes,focalNodes):       
        """Returns souce strengths assigned to focal nodes by user.        

        NOT IMPLEMENTED CURRENTLY

        """  
        variableSources = ones(numnodes,dtype = 'float64')
        if self.options.use_variable_source_strengths==True:
            pointStrengths = self.readPointStrengths(self.options.variable_source_file) 
            variableSourceNames = pointStrengths[:,0]
            variableSourceNodes = self.namesToNodes(focalNodes,variableSourceNames)
            try:
                for i in range (0,len(variableSourceNames)):
                    variableSources[variableSourceNodes[i]] = pointStrengths[variableSourceNodes[i],1]
            except IndexError:
                raise RuntimeError('Error assinging variable source strengths. Please make sure focal node names match.')                
        return variableSources

        
    # Not implemented at this time   
    def readSourcesGroundsNetwork(self, G, nodeNames, sourceFile,groundFile):
        """Reads source and ground files for advanced network mode.
        
        NOT IMPLEMENTED CURRENTLY

        """  
        if os.path.isfile(sourceFile)==False:
            raise RuntimeError('File "'  + sourceFile + '" does not exist')   
        if os.path.isfile(groundFile)==False:
            raise RuntimeError('File "'  + groundFile + '" does not exist')               
        rawSources = loadtxt(sourceFile, dtype = 'float64')
        rawGrounds = loadtxt(groundFile, dtype = 'float64')

        if self.options.ground_file_is_resistances==True:
            rawGrounds[:,1] = 1 / rawGrounds[:,1]
        if self.options.use_direct_grounds==True:
            rawGrounds[:,1] = where(rawGrounds[:,1],Inf,0)
        
        numnodes = G.shape[0]
        sourceNodes = self.namesToNodes(nodeNames,rawSources[:,0]) 
        sources = zeros((numnodes),dtype = 'float64')
        sources[sourceNodes[:]] = rawSources[:,1] 
        
        groundNodes = self.namesToNodes(nodeNames,rawGrounds[:,0]) 
        grounds = zeros((numnodes),dtype = 'float64')
        grounds[groundNodes[:]] = rawGrounds[:,1] 

        conflicts = logical_and(sources,grounds)
        if self.options.remove_src_or_gnd=='rmvsrc':
            sources = where(conflicts,0,sources)
        elif self.options.remove_src_or_gnd=='rmvgnd':
            grounds = where(conflicts,0,grounds)
        elif self.options.remove_src_or_gnd=='rmvall':
            sources = where(conflicts,0,sources)
            grounds = where(conflicts,0,grounds)
        if size(where(sources)) == 0:
            raise RuntimeError('No valid sources detected. Please check source file') 
        if size(where(grounds)) == 0:
            raise RuntimeError('No valid grounds detected. Please check ground file') 
        return sources, grounds
        
