##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: cs_large_verify.py 545 2009-05-05 23:43:31Z mcrae $
##

import unittest
import numpy as np
from cs_io import CSIO
from circuitscape import circuitscape

print 'Verifying code with Large Test Problems.'

def approxEqual(a, b):
    m = a.shape[0]
    n = a.shape[1]

    for i in range(0, m):
        for j in range(0, n):
            if (a[i,j] != b[i,j]):
##                if a[i,j]!=0:
##                    if (abs((a[i,j] - b[i,j])/a[i,j]) > 1e-6):
##                        return False
##                else:
                if (abs(a[i,j] - b[i,j]) > 1e-3):
                    return False
    return True

def cs_verifyall():
    suite = unittest.TestLoader().loadTestsFromTestCase(cs_verify)
    testResult = unittest.TextTestRunner(verbosity=0).run(suite)
    return testResult
    unittest.main()

def test_sg(ut, test_name):
    configFile='..//Large_Test_Problems//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    resistances_computed,_solver_failed = cs.compute()

    resistances_saved=np.loadtxt('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_resistances.txt') 

    current_map_1_2_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_curmap_1_2.asc', 'float64') 
    current_map_1_2_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_curmap_1_2.asc', 'float64') 

    cum_current_map_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_1_2_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_voltmap_1_2.asc', 'float64') 
    voltage_map_1_2_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_voltmap_1_2.asc', 'float64') 
    
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)
    ut.assertEquals (approxEqual(current_map_1_2_saved, current_map_1_2_computed), True)
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_1_2_saved, voltage_map_1_2_computed), True)

def test_one_to_all(ut, test_name):
    configFile='..//Large_Test_Problems//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    resistances_computed,_solver_failed = cs.compute()

    resistances_saved=np.loadtxt('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_resistances.txt') 

    current_map_1_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_curmap_1.asc', 'float64') 
    current_map_1_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_curmap_1.asc', 'float64') 

    cum_current_map_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_1_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_voltmap_1.asc', 'float64') 
    voltage_map_1_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_voltmap_1.asc', 'float64') 
    
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)
    ut.assertEquals (approxEqual(current_map_1_saved, current_map_1_computed), True)
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_1_saved, voltage_map_1_computed), True)
        
def test_mg(ut, test_name):
    configFile='..//Large_Test_Problems//verify//config_files//' + test_name + '.ini'
    cs = circuitscape(configFile, None)
    _voltages = cs.compute()
   
    cum_current_map_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_curmap.asc', 'float64') 
    cum_current_map_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_curmap.asc', 'float64') 

    voltage_map_computed=CSIO._reader('..//Large_Test_Problems//verify//output//' + test_name + '_voltmap.asc', 'float64') 
    voltage_map_saved=CSIO._reader('..//Large_Test_Problems//verify//baseline_results//' + test_name + '_voltmap.asc', 'float64') 
    
    ut.assertEquals (approxEqual(cum_current_map_saved, cum_current_map_computed), True)
    ut.assertEquals (approxEqual(voltage_map_saved, voltage_map_computed), True)
        
class cs_verify(unittest.TestCase):
    def test_single_ground_all_pairs_resistances_1(self):
        test_sg(self, '250k') 

    def test_single_ground_all_pairs_resistances_2(self):
        test_sg(self, '750k') 
 
    def test_single_ground_all_pairs_resistances_3(self):
        test_sg(self, '1m') 

if __name__ == '__main__':
    unittest.main()
