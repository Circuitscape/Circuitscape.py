## !/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import os

import unittest
import numpy as np
from cs_io import CSIO
from circuitscape import circuitscape

def approxEqual(a, b):
    m = a.shape[0]
    n = a.shape[1]

    for i in range(0, m):
        for j in range(0, n):
            if (a[i,j] != b[i,j]):
                if (abs(a[i,j] - b[i,j]) > 1e-6):  
                    return False
    return True

def cs_verifyall():
    suite = unittest.TestLoader().loadTestsFromTestCase(cs_verify)
    testResult = unittest.TextTestRunner(verbosity=0).run(suite)
    return testResult
    unittest.main()

def test_sg(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    resistances_computed, _solver_failed = cs.compute()
    resistances_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_resistances.txt') 
    
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)    


    if test_name!='sgVerify12': #This module tests the resistance shortcut which is only calculated when maps are not written.
        current_map_1_2_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap_1_2.asc', 'float64') 
        current_map_1_2_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap_1_2.asc', 'float64') 
        cum_current_map_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap.asc', 'float64') 
        cum_current_map_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap.asc', 'float64')     
        voltage_map_1_2_computed=CSIO._reader('.//verify//output//' + test_name + '_voltmap_1_2.asc', 'float64') 
        voltage_map_1_2_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_voltmap_1_2.asc', 'float64') 
            
        ut.assertEquals (approxEqual(current_map_1_2_saved, current_map_1_2_computed), True)
        ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
        ut.assertEquals (approxEqual(voltage_map_1_2_saved, voltage_map_1_2_computed), True)
        if os.path.isfile('.//verify//output//' + test_name + '_curmap_max.asc'): 
            max_current_map_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap_max.asc', 'float64') 
            max_current_map_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap_max.asc', 'float64')             
            ut.assertEquals (approxEqual(max_current_map_saved, max_current_map_computed), True)        
        
        
def test_network_sg(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)

    # These baseline outputs generated using rasters, with outputs written in graph format using 'write_baseline_results' option.
    resistances_computed, _solver_failed = cs.compute()
    resistances_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_resistances_3columns.txt') 
    cum_node_currents_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_node_currents.txt', 'float64') 
    cum_node_currents_computed=np.loadtxt('.//verify//output//' + test_name + '_node_currents.txt', 'float64') 
    branch_currents_0_1_computed=np.loadtxt('.//verify//output//' + test_name + '_branch_currents_0_1.txt', 'float64') 
    branch_currents_0_1_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_branch_currents_0_1.txt', 'float64')
    voltage_map_0_1_computed=np.loadtxt('.//verify//output//' + test_name + '_voltages_0_1.txt', 'float64') 
    voltage_map_0_1_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_voltages_0_1.txt', 'float64') 

    # Baseline cumulative branch currents generated using network code.
    cum_branch_currents_computed=np.loadtxt('.//verify//output//' + test_name + '_branch_currents.txt', 'float64') 
    cum_branch_currents_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_branch_currents.txt', 'float64')     

    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)        
    ut.assertEquals (approxEqual(cum_node_currents_saved, cum_node_currents_computed), True)
    ut.assertEquals (approxEqual(voltage_map_0_1_saved, voltage_map_0_1_computed), True)
    if test_name != 'sgNetworkVerify2': #need to replace this test.  Rounding error creates non-zero branch currents on some platforms and not others.
        ut.assertEquals (approxEqual(branch_currents_0_1_saved, branch_currents_0_1_computed), True)
    
    ut.assertEquals (approxEqual(cum_branch_currents_saved, cum_branch_currents_computed), True)


def test_one_to_all(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    resistances_computed, _solver_failed = cs.compute()

    resistances_saved=np.loadtxt('.//verify//baseline_results//' + test_name + '_resistances.txt') 

    current_map_1_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap_1.asc', 'float64') 
    current_map_1_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap_1.asc', 'float64') 

    cum_current_map_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_1_computed=CSIO._reader('.//verify//output//' + test_name + '_voltmap_1.asc', 'float64') 
    voltage_map_1_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_voltmap_1.asc', 'float64') 
    
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)
    ut.assertEquals (approxEqual(current_map_1_saved, current_map_1_computed), True)
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_1_saved, voltage_map_1_computed), True)

   
    
def test_all_to_one(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    _resistances_computed, _solver_failed = cs.compute()

    current_map_1_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap_1.asc', 'float64') 
    current_map_1_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap_1.asc', 'float64') 

    cum_current_map_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_1_computed=CSIO._reader('.//verify//output//' + test_name + '_voltmap_1.asc', 'float64') 
    voltage_map_1_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_voltmap_1.asc', 'float64') 
    
    ut.assertEquals (approxEqual(current_map_1_saved, current_map_1_computed), True)
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_1_saved, voltage_map_1_computed), True)


        
def test_mg(ut, test_name):
    #print test_name
    configFile='.//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    _voltages = cs.compute()
   
    cum_current_map_computed=CSIO._reader('.//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_computed=CSIO._reader('.//verify//output//' + test_name + '_voltmap.asc', 'float64') 
    voltage_map_saved=CSIO._reader('.//verify//baseline_results//' + test_name + '_voltmap.asc', 'float64') 
    
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_saved, voltage_map_computed), True)


class cs_verify(unittest.TestCase):

    def test_network_pairwise_1(self):
        # Baseline test case generated with points8_mod.asc and cellmap8_mod.asc
        # Settings: habitat map is resistances, 8 neighbor connection
        test_network_sg(self, 'sgNetworkVerify1') 
 
    def test_network_pairwise_2(self):
        # Baseline test case generated with verify\2 test case
        # Settings: habitat map is conductances, 4 neighbors
        test_network_sg(self, 'sgNetworkVerify2') 
         
    def test_network_pairwise_3(self):
        # Baseline test case generated with verify\4\cellmap5x5 and points5x5
        # Settings: conductances/8N with resulting baseline graph converted to resistances.  
        test_network_sg(self, 'sgNetworkVerify3') 
         
    def test_single_ground_all_pairs_resistances_1(self):
        test_sg(self, 'sgVerify1') 
 
    def test_single_ground_all_pairs_resistances_2(self):
        test_sg(self, 'sgVerify2') 
 
    def test_single_ground_all_pairs_resistances_3(self):
        test_sg(self, 'sgVerify3') 
 
    def test_single_ground_all_pairs_resistances_4(self):
        test_sg(self, 'sgVerify4') 
 
    def test_single_ground_all_pairs_resistances_5(self):
        test_sg(self, 'sgVerify5') 
 
    def test_single_ground_all_pairs_resistances_6(self):
        test_sg(self, 'sgVerify6') 
 
    def test_single_ground_all_pairs_resistances_7(self):
        test_sg(self, 'sgVerify7') 
 
    def test_single_ground_all_pairs_resistances_8(self):
        test_sg(self, 'sgVerify8') 
 
    def test_single_ground_all_pairs_resistances_9(self):
        test_sg(self, 'sgVerify9') 
 
    def test_single_ground_all_pairs_resistances_10(self):
        test_sg(self, 'sgVerify10')
 
    def test_single_ground_all_pairs_resistances_11(self):
        test_sg(self, 'sgVerify11') 
     
    def test_single_ground_all_pairs_resistances_12(self):
        test_sg(self, 'sgVerify12') 
        
#     def test_single_ground_all_pairs_resistances_13(self):
#         test_sg(self, 'sgVerify13')         
    # TODO: correct the method name and result data 
    def test_single_ground_all_pairs_resistances_13(self):
        # Tests nodata output and max current options
        test_sg(self, 'sgVerify14') 
         
         
    def test_multiple_ground_module_1(self):
        test_mg(self, 'mgVerify1') 
       
    def test_multiple_ground_module_2(self):
        test_mg(self, 'mgVerify2') 
 
    def test_multiple_ground_module_3(self):
        test_mg(self, 'mgVerify3') 
         
    def test_multiple_ground_module_4(self):
        test_mg(self, 'mgVerify4')         
 
    def test_multiple_ground_module_5(self):
        test_mg(self, 'mgVerify5')         
   
    def test_one_to_all_module_1(self):
        test_one_to_all(self, 'oneToAllVerify1') 
   
    def test_one_to_all_module_2(self):
        test_one_to_all(self, 'oneToAllVerify2') 
 
    def test_one_to_all_module_3(self):
        test_one_to_all(self, 'oneToAllVerify3') 
 
    def test_one_to_all_module_4(self):
        test_one_to_all(self, 'oneToAllVerify4') 
 
    def test_one_to_all_module_5(self):
        test_one_to_all(self, 'oneToAllVerify5') 
   
    def test_one_to_all_module_6(self):
        test_one_to_all(self, 'oneToAllVerify6') 
 
    def test_one_to_all_module_7(self):
        test_one_to_all(self, 'oneToAllVerify7') 
 
    def test_one_to_all_module_8(self):
        test_one_to_all(self, 'oneToAllVerify8') 
 
    def test_one_to_all_module_9(self):
        test_one_to_all(self, 'oneToAllVerify9') 
   
    def test_one_to_all_module_10(self):
        test_one_to_all(self, 'oneToAllVerify10') 
 
    def test_one_to_all_module_11(self):
        test_one_to_all(self, 'oneToAllVerify11') 
 
    def test_one_to_all_module_12(self):
        test_one_to_all(self, 'oneToAllVerify12') 
         
    def test_all_to_one_module_1(self):
        test_all_to_one(self, 'allToOneVerify1') 
   
    def test_all_to_one_module_2(self):
        test_all_to_one(self, 'allToOneVerify2') 
 
    def test_all_to_one_module_3(self):
        test_all_to_one(self, 'allToOneVerify3') 
 
    def test_all_to_one_module_4(self):
        test_all_to_one(self, 'allToOneVerify4') 
 
    def test_all_to_one_module_5(self):
        test_all_to_one(self, 'allToOneVerify5') 
   
    def test_all_to_one_module_6(self):
        test_all_to_one(self, 'allToOneVerify6') 
 
    def test_all_to_one_module_7(self):
        test_all_to_one(self, 'allToOneVerify7') 
 
    def test_all_to_one_module_8(self):
        test_all_to_one(self, 'allToOneVerify8') 
 
    def test_all_to_one_module_9(self):
        test_all_to_one(self, 'allToOneVerify9') 
   
    def test_all_to_one_module_10(self):
        test_all_to_one(self, 'allToOneVerify10') 
 
    def test_all_to_one_module_11(self):
        test_all_to_one(self, 'allToOneVerify11')   
         
    def test_all_to_one_module_12(self):
        test_all_to_one(self, 'allToOneVerify12')           
         
            
if __name__ == '__main__':
    unittest.main()
