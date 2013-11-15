from multiprocessing.pool import ThreadPool
import logging

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
        self.points_rc = None
        self.poly_map = None
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
        #logging.debug("creating worker pool")
        self.worker_pool = ThreadPool(max_workers if (max_workers > 0) else None)
    
    def worker_pool_wait(self):
        if self.worker_pool != None:
            self.worker_pool.close()
            self.worker_pool.join()
            self.worker_pool = None
            #logging.debug("closing worker pool")
