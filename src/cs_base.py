##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import sys, time, gc, traceback, logging, inspect
import numpy as np
from scipy.sparse.linalg import cg
from scipy import sparse
from pyamg import smoothed_aggregation_solver

from cs_cfg import CSConfig
from cs_state import CSState

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
    
        
    def logCompleteJob(self):
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
