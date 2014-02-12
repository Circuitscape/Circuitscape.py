from multiprocessing.pool import ThreadPool
import multiprocessing, os, logging, pickle
from profiler import gc_after, print_rusage
from pyamg import smoothed_aggregation_solver

class CSState:
    logger = None
    
    def __init__(self):
        self.worker_pool = None
        self.amg_hierarchy = None
        self.cellsize = None
        self.g_map = None
        self.ground_map = None
        self.included_pairs = None
        self.mask = None
        self.ncols = None
        self.nodata = None
        self.nrows = None
        self.point_strengths = None
        self.points_rc = None                           # focal node map
        self.point_file_contains_polygons = False
        self.poly_map = None                            # habitat short circuit region map
        self.source_map = None
        self.start_time = None
        self.version = None
        self.xllcorner = None
        self.yllcorner = None

    def worker_pool_create(self, max_workers, reuse=False):
        if(None != self.worker_pool):
            if reuse:
                return
            else:
                raise RuntimeError("Existing worker pool not closed")
        #CSState.logger.debug("creating worker pool")
        max_workers = min(multiprocessing.cpu_count(), max_workers)
        self.worker_pool = ThreadPool(max_workers if (max_workers > 0) else None)


    def worker_pool_submit(self, function, callback, *args):
        if(None != self.worker_pool):
            self.worker_pool.apply_async(self.async(function), args=args, callback=callback)
        else:
            try:
                result = function(*args)
            except:
                result = None
            callback(result)
            

    def worker_pool_wait(self):
        if self.worker_pool != None:
            self.worker_pool.close()
            self.worker_pool.join()
            self.worker_pool = None
            #CSState.logger.debug("closing worker pool")
        
    @gc_after
    def del_amg_hierarchy(self):
        if self.amg_hierarchy != None:
            self.amg_hierarchy = None

    @print_rusage
    def create_amg_hierarchy(self, G, solver): 
        """Creates AMG hierarchy."""  
        if solver in ['amg', 'cg+amg']:
            # construct the MG hierarchy
            #  scipy.io.savemat('c:\\temp\\graph.mat',mdict={'d':G})
            self.amg_hierarchy = smoothed_aggregation_solver(G)


    def async(self, func):
        def wrapper(*args):
            read_fd, write_fd = os.pipe()
            child_pid = os.fork()
            if (0 == child_pid):
                pid = str(os.getpid())
                logging.disable(logging.CRITICAL) # disable logging in child to avoid python bug : http://bugs.python.org/issue6721
                try:
                    os.close(read_fd)
                    write_file = os.fdopen(write_fd, 'w')
                    result = func(*args)
                except Exception as e:
                    result = e
                write_file.write(pickle.dumps(result))
                write_file.close()
                os._exit(os.EX_OK)
            elif (child_pid > 0):
                pid = str(child_pid)
                try:
                    #self.logger.debug("parallel: waiting for " + pid)
                    os.close(write_fd)
                    read_file = os.fdopen(read_fd)
                    result = read_file.read()
                    results = pickle.loads(result)
                    #self.logger.debug("parallel: got results from " + pid)
                    if isinstance(results, Exception):
                        self.logger.exception("parallel: got error from " + pid + ": " + str(results))
                        results = None
                except:
                    self.logger.exception("parallel: exception waiting for results from " + pid)
                    results = None
                finally:
                    read_file.close()
                    #self.logger.debug("parallel: waiting for termination of " + pid)
                    os.waitpid(child_pid, 0)
                    #self.logger.debug("parallel: terminated " + pid)
            else:
                self.logger.error("parallel: unable to create new processes")
                results = None
            return results
        wrapper.__name__ = 'async_' + func.func_name
        return wrapper
