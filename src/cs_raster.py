#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import time, gc, logging
import numpy
from numpy import *
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from cs_base import CSBase, print_timing#, elapsed_time, deleterow, deleterowcol, relabel
from cs_io import CSIO


class CSRaster(CSBase):
    def __init__(self, configFile, logger_func):
        super(CSRaster, self).__init__(configFile, logger_func)

    @print_timing
    def compute_raster(self):
        """Main function for Circuitscape."""  
        
        self.load_maps()
        if self.options.screenprint_log == True:        
            num_nodes = (where(self.state.g_map > 0, 1, 0)).sum()         
            logging.debug('Resistance/conductance map has %d nodes' % (num_nodes,))

        if self.options.scenario == 'pairwise':
            resistances, solver_failed = self.pairwise_module(self.state.g_map, self.state.poly_map, self.state.points_rc)
            self.logCompleteJob()
            return resistances,solver_failed     

        elif self.options.scenario == 'advanced':
            self.options.write_max_cur_maps = False
            self.log ('Calling solver module.', 1)
            voltages, _current_map, solver_failed = self.advanced_module(self.state.g_map, self.state.poly_map, self.state.source_map, self.state.ground_map,None,None,None,None,None)
            self.logCompleteJob()
            if solver_failed == True:
                logging.error('Solver failed')
            return voltages, solver_failed

        else:
            resistance_vector, solver_failed = self.one_to_all_module(self.state.g_map, self.state.poly_map, self.state.points_rc)
            self.logCompleteJob()
            return resistance_vector, solver_failed 
            
    
  
    def get_overlap_polymap(self, point, point_map, poly_map_temp, new_poly_num): 
        """Creates a map of polygons (aka short-circuit or zero resistance regions) overlapping a focal node."""  
        point_poly = numpy.where(point_map == point, 1, 0) 
        poly_point_overlap = numpy.multiply(point_poly, poly_map_temp)
        overlap_vals = numpy.unique(numpy.asarray(poly_point_overlap))
        rows = numpy.where(overlap_vals > 0)
        overlap_vals = overlap_vals[rows] #LIST OF EXISTING POLYGONS THAT OVERLAP WITH POINT
        for a in range (0, overlap_vals.size):
            poly_map_temp = numpy.where(poly_map_temp==overlap_vals[a], new_poly_num, poly_map_temp)
        poly_map_temp = numpy.where(point_map == point, new_poly_num, poly_map_temp)
        return poly_map_temp


    def append_names_to_resistances(self, point_ids, resistances):        
        """Adds names of focal nodes to resistance matrices."""  
        focal_labels = insert(point_ids, [0], 0, axis = 0)
        resistances = insert(resistances, [0], 0, axis = 0)
        resistances = insert(resistances, [0], 0, axis = 1)
        resistances[0,:] = focal_labels
        resistances[:,0] = focal_labels
        return resistances
        
    
    @print_timing
    def one_to_all_module(self, g_map, poly_map, points_rc):
        """Overhead module for one-to-all AND all-to-one modes with raster data."""  
        lastWriteTime = time.time()

        if self.options.use_included_pairs==True: #Prune points
            points_rc = self.pruneIncludedPairs(points_rc)          
            includedPairs = self.state.included_pairs 
        point_ids = unique(asarray(points_rc[:,0]))
        points_rc_unique = self.get_points_rc_unique(point_ids,points_rc)      
        
        resistance_vector = zeros((point_ids.size,2),float)
        solver_failed_somewhere = False
        if self.options.write_cur_maps == True:
            cum_current_map = zeros((self.state.nrows, self.state.ncols),dtype = 'float64')         
        if self.options.write_max_cur_maps==True:
            max_current_map=cum_current_map
        oneToAllStreamline = False        
        if self.options.use_included_pairs==False: #Will do this each time later if using included pairs
            point_map = numpy.zeros((self.state.nrows, self.state.ncols),int)
            point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]
       
            #combine point map and poly map
            poly_map_temp = self.get_poly_map_temp(poly_map,point_map,point_ids,None,None)
            unique_point_map = numpy.zeros((self.state.nrows, self.state.ncols),int)
            unique_point_map[points_rc_unique[:,1],points_rc_unique[:,2]] = points_rc_unique[:,0]

            (strengthMap,strengths_rc) = self.getstrengthMap(points_rc_unique,self.state.point_strengths)            

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

            self.log ('solving focal node ' + str(i+1) + ' of ' + str(point_ids.size) + '.',1)

            if self.options.use_included_pairs==True: # Done above otherwise    
                #######################   
                points_rc_unique_temp = numpy.copy(points_rc_unique)
                point_map = numpy.zeros((self.state.nrows, self.state.ncols),int)
                point_map[points_rc[:,1],points_rc[:,2]] = points_rc[:,0]       

                for pair in range(0, point_ids.size): #loop thru exclude[point,:], delete included pairs of focal point from point_map and points_rc_unique_temp 
                    if includedPairs[i+1,pair+1]==0 and i !=  pair:
                            pt_id = point_ids[pair]
                            point_map = where(point_map==pt_id,0,point_map)
                            points_rc_unique_temp[pair,0] = 0 #point will not be burned in to unique_point_map

                poly_map_temp = self.get_poly_map_temp2(poly_map,point_map,points_rc_unique_temp,includedPairs,i)

                unique_point_map = numpy.zeros((self.state.nrows, self.state.ncols),int)
                unique_point_map[points_rc_unique_temp[:,1],points_rc_unique_temp[:,2]] = points_rc_unique_temp[:,0]        

                (strengthMap,strengths_rc) = self.getstrengthMap(points_rc_unique_temp,self.state.point_strengths)
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
                    
                resistance,current_map,solver_failed = self.advanced_module(self.state.g_map, poly_map_temp, source_map, ground_map, src, G, node_map, component_map, componentWithPoints)
                if solver_failed == False:
                    if self.options.write_cur_maps == True:
                        cum_current_map = cum_current_map+current_map
                        if self.options.write_max_cur_maps==True:
                            max_current_map=maximum(max_current_map,current_map)
                else:
                    logging.warning('Solver failed for at least one focal node.  \nFocal nodes with failed solves will be marked with value of -777 \nin output resistance list.\n')
    
                resistance_vector[i,0] = src
                resistance_vector[i,1] = resistance
                    
                if solver_failed==True:
                    solver_failed_somewhere = True
            else:
                resistance_vector[i,0] = src
                resistance_vector[i,1] = -1            

            (hours,mins,_secs) = self.elapsed_time(lastWriteTime)
            if mins > 2 or hours > 0: 
                lastWriteTime = time.time()
                CSIO.write_resistances_one_to_all(self.options.output_file, resistance_vector, '_incomplete', self.options.scenario)
                #self.writeResistancesOneToAll(resistance_vector,'_incomplete')
      
        if solver_failed_somewhere==False:
            if self.options.write_cur_maps == True:
                CSIO.write_aaigrid('cum_curmap', '', cum_current_map, self.options, self.state)
                if self.options.write_max_cur_maps==True:
                    CSIO.write_aaigrid('max_curmap', '', max_current_map, self.options, self.state)

        CSIO.write_resistances_one_to_all(self.options.output_file, resistance_vector, '', self.options.scenario)
        #self.writeResistancesOneToAll(resistance_vector,'')
       
        return resistance_vector,solver_failed_somewhere 

    def get_poly_map_temp(self, poly_map, point_map, point_ids, includedPairs, point1):
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

        
    def get_poly_map_temp2(self, poly_map, point_map, points_rc_unique_temp, includedPairs, i):
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


    def _pairwise_module_alloc_current_maps(self):
        cum_current_map = max_current_map = []
        if self.options.write_cur_maps == True:
            cum_current_map = zeros((self.state.nrows, self.state.ncols), dtype='float64') 
            if self.options.write_max_cur_maps == True:
                max_current_map = cum_current_map
        return (cum_current_map, max_current_map)
    
    @print_timing   
    def pairwise_module(self, g_map, poly_map, points_rc):
        """Overhead module for pairwise mode with raster data."""  
        cum_current_map, max_current_map = self._pairwise_module_alloc_current_maps()

        # If there are no focal regions, pass all points to single_ground_all_pair_resistances,
        # otherwise, pass one point at a time.
        if self.options.point_file_contains_polygons == False:
            if points_rc.shape[0] != (unique(asarray(points_rc[:,0]))).shape[0]:
                raise RuntimeError('At least one focal node contains multiple cells.  If this is what you really want, then choose focal REGIONS in the pull-down menu') 
            
            if self.options.use_included_pairs == True:
                points_rc = self.pruneIncludedPairs(points_rc)
            
            reportStatus = True
            try:
                (resistances, cum_current_map, max_current_map, solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map, points_rc, cum_current_map, max_current_map, reportStatus)
            except MemoryError: #Give it a try, but starting again never seems to helps even from GUI.
                self.enable_low_memory(True) #This doesn't seem to really clear out memory or truly restart.
                cum_current_map, max_current_map = self._pairwise_module_alloc_current_maps()
                #Note: This does not go through when it should.
                (resistances, cum_current_map, max_current_map, solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map, points_rc,cum_current_map,max_current_map,reportStatus)
                
            if solver_failed == True:
                logging.warning('Solver failed for at least one focal node pair. ' 
                '\nThis can happen when input resistances differ by more than' 
                '\n~6 orders of magnitude. Pairs with failed solves will be '
                '\nmarked with value of -777 in output resistance matrix.\n')

            point_ids = points_rc[:,0]

        else:
            if self.options.use_included_pairs == True:
                points_rc = self.pruneIncludedPairs(points_rc)
                includedPairs = self.state.included_pairs
            else:
                numpoints = points_rc.shape[0]
                point_ids = unique(asarray(points_rc[:,0]))
                points_rc_unique = self.get_points_rc_unique(point_ids, points_rc)#Nov13_2010   #Fixme: can just use point ids for index size
                num_unique_points = points_rc_unique.shape[0]#Nov13_2010
                includedPairs = ones((num_unique_points+1, num_unique_points+1), dtype='int32')#Nov13_2010

            point_map = numpy.zeros((self.state.nrows, self.state.ncols), int)
            point_map[points_rc[:,1], points_rc[:,2]] = points_rc[:,0]

            point_ids = unique(asarray(points_rc[:,0]))
            points_rc_unique = self.get_points_rc_unique(point_ids, points_rc)

            resistances = -1 * numpy.ones((point_ids.size, point_ids.size), dtype='float64')
            x = 0
            for i in range(0, point_ids.size-1):
                for j in range(i+1, point_ids.size):
                    if includedPairs[i+1,j+1] != 1:
                        continue

                    if poly_map == []:
                        poly_map_temp = zeros((self.state.nrows, self.state.ncols),int)
                        new_poly_num = 1
                    else:
                        poly_map_temp = poly_map
                        new_poly_num = numpy.max(poly_map)+1
                    #point = point_ids[i]
                    poly_map_temp = self.get_overlap_polymap(point_ids[i], point_map, poly_map_temp, new_poly_num) 
                    poly_map_temp = self.get_overlap_polymap(point_ids[j], point_map, poly_map_temp, new_poly_num+1) 

                    #Get first instance of each point in points_rc
                    points_rc_temp = numpy.zeros((2,3), int)
                    points_rc_temp[0,:] = points_rc_unique[i,:]
                    points_rc_temp[1,:] = points_rc_unique[j,:]

                    numpoints = point_ids.size
                    x = x+1
                    y = numpoints*(numpoints-1)/2
                    self.log ('solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
                    reportStatus = False
                    
                    (pairwise_resistance, cum_current_map, max_current_map, solver_failed) = self.single_ground_all_pair_resistances(g_map, poly_map_temp, points_rc_temp,cum_current_map,max_current_map,reportStatus)

                    del poly_map_temp
                    if solver_failed == True:
                        logging.warning('Solver failed for at least one focal node pair.  \nPairs with failed solves will be marked with value of -777 \nin output resistance matrix.\n')

                    resistances[i,j] = pairwise_resistance[0,1]
                    resistances[j,i] = pairwise_resistance[0,1]

        for i in range(0,resistances.shape[0]): #Set diagonal to zero
            resistances[i, i] = 0

        #Add row and column headers and write resistances to disk
        resistances = self.writeResistances(point_ids, resistances)

        if self.options.write_cur_maps == True:
            if solver_failed == False:
                if self.options.log_transform_maps == True:
                    cum_current_map = where(cum_current_map>0, log10(cum_current_map), self.state.nodata) 
                    if self.options.write_max_cur_maps==True:
                        max_current_map = where(max_current_map>0, log10(max_current_map), self.state.nodata) 
                CSIO.write_aaigrid('cum_curmap', '', cum_current_map, self.options, self.state)
                
                if self.options.write_max_cur_maps == True:      
                    CSIO.write_aaigrid('max_curmap', '', max_current_map, self.options, self.state)

        return resistances,solver_failed
    
    
    @print_timing
    def single_ground_all_pair_resistances(self, g_map, poly_map, points_rc, cum_current_map, max_current_map, report_status):
        """Handles pairwise resistance/current/voltage calculations.  
        
        Called once when focal points are used, called multiple times when focal regions are used.
        """  
        last_write_time = time.time()
        numpoints = points_rc.shape[0]
        if (self.options.use_included_pairs==False) or (self.options.point_file_contains_polygons==True):
            included_pairs = ones((numpoints+1,numpoints+1), dtype='int32')
        else:
            included_pairs = self.state.included_pairs
        
        if (self.options.point_file_contains_polygons==True) or  (self.options.write_cur_maps == True) or (self.options.write_volt_maps == True) or (self.options.use_included_pairs==True): 
            use_resistance_calc_shortcut = False
        else:     
            use_resistance_calc_shortcut = True # We use this when there are no focal regions.  It saves time when we are also not creating maps
            shortcut_resistances = -1 * ones((numpoints, numpoints), dtype='float64') 
           
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
                if (use_resistance_calc_shortcut==True):
                    voltmatrix = zeros((numpoints,numpoints),dtype = 'float64')     #For resistance calc shortcut
                
                dstPoint = 0
                anchorPoint = 0 #For resistance calc shortcut
                
                ##############
                for i in range(0, numpoints):
                    if range(i, numpoints) == []:
                        break

                    if (use_resistance_calc_shortcut==True) and (dstPoint>0): 
                        break #No need to continue, we've got what we need to calculate resistances

                    dst = self.grid_to_graph (points_rc[i,1], points_rc[i,2], node_map)
                    local_dst = self.grid_to_graph (points_rc[i,1], points_rc[i,2], local_node_map)
                    if (dst >=  0 and components[dst] == c):
                        dstPoint = dstPoint+1
                        #Gsolve = []
                    
                        G_dst_dst = G[local_dst, local_dst] 
                        G[local_dst,local_dst] = 0
    
                        self.state.amg_hierarchy = None
                        gc.collect()
                        self.create_amg_hierarchy(G)

                        ################    
                        for j in range(i+1, numpoints):
                            if included_pairs[i+1,j+1]==1: #Test for pair in included_pairs
                                if self.state.amg_hierarchy==None: #Called in case of memory error in current mapping
                                    self.create_amg_hierarchy(G)
                               
                                # tan: parallelize here
                                if report_status==True:
                                    x = x+1
                                    if use_resistance_calc_shortcut==True:
                                        y = numpoints
                                        self.log ('solving focal node ' + str(x) + ' of '+ str(y) + '.',1)
                                    else:
                                        y = numpoints*(numpoints-1)/2
                                        self.log ('solving focal pair ' + str(x) + ' of '+ str(y) + '.',1)
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
                                            self.state.amg_hierarchy = None
                                            gc.collect()    
                                        
                                        resistances[i, j] = voltages[local_src] - voltages[local_dst]
                                        resistances[j, i] = voltages[local_src] - voltages[local_dst]
                                        # Write maps to files
                                        frompoint = str(points_rc[i,0])
                                        topoint = str(points_rc[j,0])
                                        
                                        if use_resistance_calc_shortcut==True:
                                            if dstPoint==1: #this occurs for first i that is in component
                                                anchorPoint = i #for use later in shortcult resistance calc
                                                voltmatrix = self.getVoltmatrix(i,j,numpoints,local_node_map,voltages,points_rc,resistances,voltmatrix)                                          

                                        if self.options.write_volt_maps == True:
                                            if report_status==True:
                                                self.log ('writing voltage map ' + str(x) + ' of ' + str(y) + '.',1)
                                            voltage_map = self.create_voltage_map(local_node_map,voltages) 
                                            CSIO.write_aaigrid('voltmap', '_' + frompoint + '_' + topoint, voltage_map, self.options, self.state)
                                            del voltage_map
                                        if self.options.write_cur_maps == True:
                                            if report_status==True:
                                                self.log ('writing current map ' + str(x) + ' of ' + str(y) + '.',1)
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
                                                    current_map = where(current_map>0,log10(current_map),self.state.nodata)
                                                CSIO.write_aaigrid('curmap', '_' + frompoint + '_' + topoint, current_map, self.options, self.state)
                                            del current_map    

                                        (hours,mins,_secs) = self.elapsed_time(last_write_time)
                                        if mins > 2 or hours > 0: 
                                            last_write_time = time.time()
                                            CSIO.save_incomplete_resistances(self.options.output_file, resistances)# Save incomplete resistances
                        if (use_resistance_calc_shortcut==True and i==anchorPoint): # This happens once per component. Anchorpoint is the first i in component
                            shortcut_resistances = self.getShortcutResistances(anchorPoint,voltmatrix,numpoints,resistances,shortcut_resistances)
                                                
                        G[local_dst, local_dst] = G_dst_dst

                    #End for
                    self.state.amg_hierarchy = None
                    gc.collect()

                #End if
                self.state.amg_hierarchy = None
                gc.collect()

        # Finally, resistance to self is 0.
        if use_resistance_calc_shortcut==True: 
            resistances = shortcut_resistances
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
            cum_current_map = zeros((self.state.nrows, self.state.ncols),dtype = 'float64')         
        
        if self.options.write_volt_maps == True: 
            cum_voltage_map = zeros((self.state.nrows, self.state.ncols),dtype = 'float64') 
        elif self.options.scenario=='one-to-all':
            cum_voltage_map = zeros((self.state.nrows, self.state.ncols),dtype = 'float64') 
        
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
                CSIO.write_aaigrid(filetext, fileadd, cum_voltage_map, self.options, self.state)  

        if self.options.scenario=='advanced':
            if self.options.write_cur_maps == True: 
                if solver_failed==False:
                    if self.options.log_transform_maps==True:
                        cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state.nodata) 
                    filetext = 'curmap'
                    fileadd = ''
                    CSIO.write_aaigrid(filetext, fileadd, cum_current_map, self.options, self.state) 
            else:
                cum_current_map = None

        else:
            if self.options.write_cur_maps == True: 
                if self.options.write_cum_cur_map_only==False:
                    if solver_failed==False:
                        if self.options.log_transform_maps==True:
                            cum_current_map = where(cum_current_map>0,log10(cum_current_map),self.state.nodata) 
                        filetext = 'curmap'   
                        fileadd = '_'+str(source_id)   
                        CSIO.write_aaigrid(filetext, fileadd, cum_current_map, self.options, self.state) 
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
        Gsolve = self.deleterowcol(Gsolve, delrow = dst_to_delete, delcol = dst_to_delete)
        
        self.create_amg_hierarchy(Gsolve)
        voltages = self.solve_linear_system(Gsolve, sources)
        del Gsolve
        self.state.amg_hierarchy = None

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
        node_map = numpy.zeros(g_map.shape, dtype = 'int32')
        node_map[g_map.nonzero()] = numpy.arange(1, numpy.sum(g_map>0)+1, dtype='int32')

        if poly_map == []:
            return node_map

        # Remove extra points from poly_map that are not in g_map
        poly_map_pruned = numpy.zeros(g_map.shape, dtype='int32')
        poly_map_pruned[numpy.where(g_map)] = poly_map[numpy.where(g_map)]
        
        polynums = numpy.unique(poly_map)
   
        for i in range(0, polynums.size):
            polynum = polynums[i]
            if polynum !=  0:
                (pi, pj) = numpy.where(poly_map_pruned == polynum) #
                (pk, pl) = numpy.where(poly_map == polynum) #Added 040309 BHM                
                if len(pi) > 0:  
                    node_map[pk, pl] = node_map[pi[0], pj[0]] #Modified 040309 BHM  
        node_map[numpy.where(node_map)] = self.relabel(node_map[numpy.where(node_map)], 1) #BHM 072409

        #print "point_file = %s"%(self.options.point_file,)
        #print "node_map ="
        #print node_map
        return node_map

    @print_timing
    def construct_component_map(self, g_map, node_map):
        """Assigns component numbers to grid corresponding to pixels with non-zero conductances.
        
        Nodes with the same component number are in single, connected components.
        
        """  
        prunedMap = False
        G = self.construct_g_graph(g_map, node_map, prunedMap) 
        (_numComponents, C) = connected_components(G)
        C += 1 # Number components from 1

        (I, J) = where(node_map)
        nodes = node_map[I, J].flatten()

        component_map = zeros(node_map.shape, dtype = 'int32')
        component_map[I, J] = C[nodes-1]

        #print "component_map:"
        #print component_map
        #print "component_map C:"
        #print C
        return (component_map, C)


    @print_timing
    def construct_g_graph(self, g_map, node_map,prunedMap):
        """Construct sparse adjacency matrix given raster maps of conductances and nodes."""  
        numnodes = node_map.max()
        (node1, node2, conductances) = self.get_conductances(g_map, node_map)
        G = sparse.csr_matrix((conductances, (node1, node2)), shape = (numnodes, numnodes)) # Memory hogging operation?
        g_graph = G + G.T
        
        #print "g_graph ="
        #print g_graph
        return g_graph

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
    def create_voltage_map(self, node_map, voltages):
        """Creates raster map of voltages given node voltage vector."""  
        voltage_map = numpy.zeros((self.state.nrows, self.state.ncols), dtype = 'float64')
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
        current_map = zeros((self.state.nrows, self.state.ncols),dtype = 'float64')
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
        resistances3Columns = self.convertResistances3cols(outputResistances)
        CSIO.write_resistances(self.options.output_file, outputResistances, resistances3Columns)
        return outputResistances        
        
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
        
        
    def pruneIncludedPairs(self, points_rc):
        """Remove excluded points from focal node list when using extra file that lists pairs to include/exclude."""
        included_pairs = self.state.included_pairs
        include_list = list(included_pairs[0,:])
        point = 0
        _drop_flag = False
        while point < points_rc.shape[0]: #Prune out any points not in include_list
            if points_rc[point,0] in include_list: #match
                point = point+1
            else:
                _drop_flag = True   
                points_rc = self.deleterow(points_rc, point)  
         
        include_list = list(points_rc[:,0])
        num_connection_rows = included_pairs.shape[0]
        row = 1
        while row < num_connection_rows: #Prune out any entries in include_list that are not in points_rc
            if included_pairs[row,0] in include_list: #match
                row = row+1
            else:
                included_pairs = self.deleterowcol(included_pairs, delrow=row, delcol=row)   
                _drop_flag = True
                num_connection_rows = num_connection_rows-1

        self.state.included_pairs = included_pairs
#         if _drop_flag==True:
#             print'\nNOTE: Code to exclude pairwise calculations is activated and \nsome entries did not match with focal node file.  \nSome focal nodes may have been dropped.'      
        return points_rc
    
    
    def get_points_rc_unique(self, point_ids, points_rc):
        """Return a list of unique focal node IDs and x-y coordinates."""  
        points_rc_unique = numpy.zeros((point_ids.size, 3), int)
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
                strengths_rc = self.get_strengths_rc(self.state.point_strengths,points_rc_unique)
                strengthMap = None
            else:
                strengths_rc = self.get_strengths_rc(pointStrengths,points_rc_unique)
                strengthMap = numpy.zeros((self.state.nrows,self.state.ncols),dtype = 'Float64')
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
        self.log('Reading maps', 1)
        self.log('', 2)
        reclass_file = self.options.reclass_file if self.options.use_reclass_table else None
        CSIO.read_cell_map(self.options.habitat_file, self.options.habitat_map_is_resistances, reclass_file, self.state)
        
        if self.options.use_polygons:
            self.state.poly_map = CSIO.read_poly_map(self.options.polygon_file, False, 0, self.state, True, "Short-circuit region", 'int32')
        else:
            self.state.poly_map = []
 
        if self.options.use_mask==True:
            mask = CSIO.read_poly_map(self.options.mask_file, True, 0, self.state, True, "Mask", 'int32')
            mask = where(mask !=  0, 1, 0) 
            self.state.g_map = multiply(self.state.g_map, mask)
            del mask
            
            sum_gmap = (self.state.g_map).sum()
            #sum_gmap = sum_gmap.sum()
            if sum_gmap==0:
                raise RuntimeError('All entries in habitat map have been dropped after masking with the mask file.  There is nothing to solve.')             
        else:
            self.state.mask = []

        if self.options.scenario=='advanced':
            self.state.points_rc = []
            (self.state.source_map, self.state.ground_map) = CSIO.read_source_and_ground_maps(self.options.source_file, self.options.ground_file, self.state, self.options)
        else:        
            self.state.points_rc = CSIO.read_point_map(self.options.point_file, "Focal node", self.state)
            self.state.source_map = []
            self.state.ground_map = []

        if self.options.use_included_pairs==True:
            self.state.included_pairs = CSIO.read_included_pairs(self.options.included_pairs_file)
        
        self.state.point_strengths = None
        if self.options.use_variable_source_strengths==True:
            self.state.point_strengths = CSIO.read_point_strengths(self.options.variable_source_file) 
        
        self.log('Processing maps',1)
        return 
 
