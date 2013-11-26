#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae, Viral B. Shah. and Tanmay Mohapatra

import unittest, os
import numpy as np
from circuitscape.cs_io import CSIO
from circuitscape import circuitscape

print 'Verifying code with Large Test Problems.'

BIG_TESTS_ROOT      = '.'
BIG_TESTS_CFG       = os.path.join(BIG_TESTS_ROOT, 'verify', 'config_files')
BIG_TESTS_BASELINE  = os.path.join(BIG_TESTS_ROOT, 'verify', 'baseline_results')
BIG_TESTS_OUT       = os.path.join(BIG_TESTS_ROOT, 'verify', 'output')

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
                if (abs(a[i,j] - b[i,j]) > 1e-4):
                    return False
    return True

def compare_results(ut, test_name, result_file, compressed):
    result_name = test_name + '_' + result_file
    if compressed:
        result_name += '.gz'
    computed    = CSIO._ascii_grid_reader(os.path.join(BIG_TESTS_OUT,      result_name), 'float64') 
    saved       = CSIO._ascii_grid_reader(os.path.join(BIG_TESTS_BASELINE, result_name), 'float64')
    ut.assertEquals(approxEqual(saved, computed), True) 

def cs_verifyall():
    suite = unittest.TestLoader().loadTestsFromTestCase(cs_verify)
    testResult = unittest.TextTestRunner(verbosity=0).run(suite)
    return testResult
    unittest.main()

def test_sg(ut, test_name):
    configFile = os.path.join(BIG_TESTS_CFG, test_name + '.ini')
    cs = circuitscape(configFile, None)
    resistances_computed, _solver_failed = cs.compute()

    resistances_saved = np.loadtxt(os.path.join(BIG_TESTS_BASELINE, test_name + '_resistances.txt'))
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)

    compare_results(ut, test_name, 'curmap_1_2.asc', True)
    compare_results(ut, test_name, 'curmap.asc', True)
    compare_results(ut, test_name, 'voltmap_1_2.asc', True)

def test_one_to_all(ut, test_name):
    configFile = os.path.join(BIG_TESTS_CFG, test_name + '.ini')
    cs = circuitscape(configFile, None)
    resistances_computed,_solver_failed = cs.compute()

    resistances_saved = np.loadtxt(os.path.join(BIG_TESTS_BASELINE, test_name + '_resistances.txt'))
    ut.assertEquals (approxEqual(resistances_saved, resistances_computed), True)

    compare_results(ut, test_name, 'curmap_1.asc', True)
    compare_results(ut, test_name, 'curmap.asc', True)
    compare_results(ut, test_name, 'voltmap_1.asc', True)
        
def test_mg(ut, test_name):
    configFile = os.path.join(BIG_TESTS_CFG, test_name + '.ini')
    cs = circuitscape(configFile, None)
    _voltages = cs.compute()
   
    compare_results(ut, test_name, 'curmap.asc', True)
    compare_results(ut, test_name, 'voltmap.asc', True)
        
class cs_verify(unittest.TestCase):
    def test_single_ground_all_pairs_resistances_1(self):
        test_sg(self, '250k') 

    #def test_single_ground_all_pairs_resistances_2(self):
    #    test_sg(self, '750k') 
 
    def test_single_ground_all_pairs_resistances_3(self):
        test_sg(self, '1m') 

if __name__ == '__main__':
    unittest.main()
