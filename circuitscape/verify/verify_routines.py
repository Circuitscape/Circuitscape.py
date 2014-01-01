import os, sys, unittest
import numpy as np
import circuitscape as cscape

TESTS_ROOT      = 'circuitscape'
TESTS_CFG       = os.path.join(TESTS_ROOT,  'verify', 'config_files')
TESTS_BASELINE  = os.path.join(TESTS_ROOT,  'verify', 'baseline_results')
TESTS_OUT       = os.path.join(TESTS_ROOT,  'verify', 'output')
EXT_LOGGER      = None

def approxEqual(a, b):
    m = a.shape[0]
    n = a.shape[1]

    for i in range(0, m):
        for j in range(0, n):
            if (a[i,j] != b[i,j]):
                if (abs(a[i,j] - b[i,j]) > 1e-6):  
                    return False
    return True

def compare_results(ut, test_name, result_file, compressed):
    result_name = test_name + '_' + result_file
    if compressed:
        result_name += '.gz'
        
    if result_file.endswith("asc"):
        computed    = cscape.CSIO._ascii_grid_reader(os.path.join(TESTS_OUT,      result_name), 'float64') 
        saved       = cscape.CSIO._ascii_grid_reader(os.path.join(TESTS_BASELINE, result_name), 'float64')
    else:
        computed    = np.loadtxt(os.path.join(TESTS_OUT,      result_name), 'float64')
        saved       = np.loadtxt(os.path.join(TESTS_BASELINE, result_name), 'float64')
    ut.assertEquals(approxEqual(saved, computed), True) 

def load_config(test_name):
    #print test_name
    configFile = os.path.join(TESTS_CFG, test_name + '.ini')
    cs = cscape.Compute(configFile, EXT_LOGGER)
    cs.logger.info("Running test: " + test_name)
    _out_dir, out_file = os.path.split(cs.options.output_file)
    cs.options.output_file = os.path.join(TESTS_OUT, out_file)
    return cs

def set_paths(root_path=None, out_path=None):
    global TESTS_OUT, TESTS_ROOT, TESTS_CFG, TESTS_BASELINE
    TESTS_OUT = out_path
    TESTS_ROOT = root_path if root_path else 'circuitscape'
    
    TESTS_CFG       = os.path.join(TESTS_ROOT,  'verify', 'config_files')
    TESTS_BASELINE  = os.path.join(TESTS_ROOT,  'verify', 'baseline_results')
    if None == TESTS_OUT:
        TESTS_OUT   = os.path.join(TESTS_ROOT,  'verify', 'output')

def cs_verifyall(root_path=None, out_path=None, ext_logger=None, stream=None):
    global EXT_LOGGER
    EXT_LOGGER = ext_logger
    set_paths(root_path, out_path)
    suite = unittest.TestLoader().loadTestsFromTestCase(cs_verify)
    if stream == None:
        stream = sys.stderr
    testResult = unittest.TextTestRunner(verbosity=0, stream=stream).run(suite)
    EXT_LOGGER = None
    return testResult

def test_sg(ut, test_name):
    cs = load_config(test_name)
    resistances_computed, _solver_failed = cs.compute()
    
    resistances_saved = np.loadtxt(os.path.join(TESTS_BASELINE, test_name + '_resistances.txt'))
    ut.assertEquals(approxEqual(resistances_saved, resistances_computed), True)    

    if test_name != 'sgVerify12': #This module tests the resistance shortcut which is only calculated when maps are not written.
        compare_results(ut, test_name, 'curmap_1_2.asc', False)
        compare_results(ut, test_name, 'cum_curmap.asc', False)
        compare_results(ut, test_name, 'voltmap_1_2.asc', False)
            
        if os.path.isfile('.//verify//output//' + test_name + '_curmap_max.asc'):
            compare_results(ut, test_name, 'curmap_max.asc', False) 
        
        
def test_network_sg(ut, test_name):
    cs = load_config(test_name)

    # These baseline outputs generated using rasters, with outputs written in graph format using 'write_baseline_results' option.
    resistances_computed, _solver_failed = cs.compute()
    
    resistances_saved = np.loadtxt(os.path.join(TESTS_BASELINE, test_name + '_resistances_3columns.txt'))
    ut.assertEquals(approxEqual(resistances_saved, resistances_computed), True)
    
    compare_results(ut, test_name, 'node_currents_cum.txt', False)
    compare_results(ut, test_name, 'voltages_0_1.txt', False)
    compare_results(ut, test_name, 'branch_currents_cum.txt', False)

    if test_name != 'sgNetworkVerify2': #need to replace this test.  Rounding error creates non-zero branch currents on some platforms and not others.
        compare_results(ut, test_name, 'branch_currents_0_1.txt', False)


def test_network_mg(ut, test_name):
    cs = load_config(test_name)
        
    _voltages = cs.compute()

    compare_results(ut, test_name, 'node_currents.txt', False)
    compare_results(ut, test_name, 'voltages.txt', False)
    compare_results(ut, test_name, 'branch_currents.txt', False)


def test_one_to_all(ut, test_name):
    cs = load_config(test_name)

    resistances_computed, _solver_failed = cs.compute()

    resistances_saved = np.loadtxt(os.path.join(TESTS_BASELINE, test_name + '_resistances.txt'))
    ut.assertEquals(approxEqual(resistances_saved, resistances_computed), True)

    compare_results(ut, test_name, 'curmap_1.asc', False)
    compare_results(ut, test_name, 'cum_curmap.asc', False)
    compare_results(ut, test_name, 'voltmap_1.asc', False)
   
    
def test_all_to_one(ut, test_name):
    cs = load_config(test_name)
    
    _resistances_computed, _solver_failed = cs.compute()

    compare_results(ut, test_name, 'curmap_1.asc', False)
    compare_results(ut, test_name, 'cum_curmap.asc', False)
    compare_results(ut, test_name, 'voltmap_1.asc', False)

        
def test_mg(ut, test_name):
    cs = load_config(test_name)
        
    _voltages = cs.compute()

    compare_results(ut, test_name, 'curmap.asc', False)
    compare_results(ut, test_name, 'voltmap.asc', False)


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

    def test_network_multiple_ground_1(self):
        # Simple test with one source and one ground
        test_network_mg(self, 'mgNetworkVerify1') 
        
    def test_network_multiple_ground_2(self):
        # Simple test with one source and five grounds. Same as Fig 2 in McRae et al. 2008.
        test_network_mg(self, 'mgNetworkVerify2') 

    def test_network_multiple_ground_3(self):
        # Modification of mgNetworkVerify2 with 4 current sources instead of 1.
        test_network_mg(self, 'mgNetworkVerify3') 
        
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
        
    def test_single_ground_all_pairs_resistances_13(self):
        test_sg(self, 'sgVerify13')         

    def test_single_ground_all_pairs_resistances_14(self):
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
         
