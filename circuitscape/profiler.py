import time, os, gc


class ResourceLogger:
    print_timings = False
    print_rusages = False    
    print_res_spaces = 0
    psutil_available = False
    resource_available = False
    print_nz_only = False

    rlogger = None
    proc = None
    proc_has_io_counters = False
    
    t1 = []
    mem1 = []
    io1 = []
    rusage1 = []

    @staticmethod
    def init_rusage(print_t=False, print_r=False, logger=None):
        ResourceLogger.print_timing_enabled(print_t)
        ResourceLogger.print_rusage_enabled(print_r)
        GCPreempt.logger = ResourceLogger.rlogger = logger
        
        if ResourceLogger.psutil_available and ResourceLogger.print_rusages:
            ResourceLogger.proc = psutil.Process(os.getpid())
            ResourceLogger.proc_has_io_counters = hasattr(ResourceLogger.proc, 'get_io_counters')

    @staticmethod
    def print_timing_enabled(is_enabled):
        """Enables or disables the print_timings decorator."""
        ResourceLogger.print_timings = is_enabled
    
    @staticmethod
    def print_rusage_enabled(is_enabled):
        """Enables or disables the print_timings decorator."""
        ResourceLogger.print_rusages = is_enabled and (ResourceLogger.psutil_available or ResourceLogger.resource_available)
    
    @staticmethod
    def is_enabled():
        return ResourceLogger.print_rusages or ResourceLogger.print_timings
    
    @staticmethod
    def do_pre():
        ResourceLogger.print_res_spaces +=  2
        ResourceLogger.t1.append(time.time())
        if ResourceLogger.print_rusages:
            if ResourceLogger.resource_available:
                ResourceLogger.rusage1.append(resource.getrusage(resource.RUSAGE_SELF))
            elif ResourceLogger.proc:
                ResourceLogger.rusage1.append(ResourceLogger.proc.get_cpu_times())
                
            if ResourceLogger.proc:
                ResourceLogger.mem1.append(ResourceLogger.proc.get_ext_memory_info())
                ResourceLogger.io1.append(ResourceLogger.proc.get_io_counters() if ResourceLogger.proc_has_io_counters else psutil.disk_io_counters())
    
    @staticmethod
    def do_post(func_name):
        t2 = time.time()
        ResourceLogger.print_res_spaces -=  2
        
        if ResourceLogger.print_rusages: # resource not available for Windows
            cpu_diffs = []
            mem_diffs = []
            io_diffs = []
            
            if ResourceLogger.resource_available:
                rusage1 = ResourceLogger.rusage1.pop()
                rusage2 = resource.getrusage(resource.RUSAGE_SELF)
                ResourceLogger.append_diff(cpu_diffs, 'user',       rusage2.ru_utime - rusage1.ru_utime)
                ResourceLogger.append_diff(cpu_diffs, 'system',     rusage2.ru_stime - rusage1.ru_stime)
                ResourceLogger.append_diff(cpu_diffs, 'elapsed',    t2-ResourceLogger.t1.pop())
            elif ResourceLogger.proc:
                rusage1 = ResourceLogger.rusage1.pop()
                rusage2 = ResourceLogger.proc.get_cpu_times()
                for attr in ['user', 'system']:
                    if hasattr(rusage1, attr):
                        ResourceLogger.append_diff(cpu_diffs, attr, getattr(rusage2, attr) - getattr(rusage1, attr))
                ResourceLogger.append_diff(cpu_diffs, 'elapsed',    t2-ResourceLogger.t1.pop())                
            
            if ResourceLogger.proc:
                mem1 = ResourceLogger.mem1.pop()
                io1 = ResourceLogger.io1.pop()
                mem2 = ResourceLogger.proc.get_ext_memory_info()
                io2 = ResourceLogger.proc.get_io_counters() if ResourceLogger.proc_has_io_counters else psutil.disk_io_counters()
                for attr in ['rss', 'vms', 'shared', 'text', 'data']:
                    if hasattr(mem1, attr):
                        ResourceLogger.append_diff(mem_diffs, attr, getattr(mem2, attr) - getattr(mem1, attr))
                for attr in ['read_count', 'write_count', 'read_bytes', 'write_bytes']:
                    if hasattr(io1, attr):
                        ResourceLogger.append_diff(io_diffs, attr, getattr(io2, attr) - getattr(io1, attr))
            elif ResourceLogger.resource_available:
                ResourceLogger.append_diff(mem_diffs, 'maxrss', rusage2.ru_maxrss - rusage1.ru_maxrss)
                ResourceLogger.append_diff(mem_diffs, 'shared', rusage2.ru_ixrss - rusage1.ru_ixrss)
                ResourceLogger.append_diff(mem_diffs, 'heap',   rusage2.ru_idrss - rusage1.ru_idrss)
                ResourceLogger.append_diff(mem_diffs, 'stack',  rusage2.ru_isrss - rusage1.ru_isrss)
                            
                ResourceLogger.append_diff(io_diffs, 'read_count',  rusage2.ru_inblock - rusage1.ru_inblock)
                ResourceLogger.append_diff(io_diffs, 'write_count', rusage2.ru_oublock - rusage1.ru_oublock)

            mem_str = 'mem(' + ' '.join(mem_diffs) + ')'
            io_str = 'io(' + ' '.join(io_diffs) + ')'
            cpu_str = 'cpu(' + ' '.join(cpu_diffs) + ')'
            log_str = ' '.join([cpu_str, mem_str, io_str])
        else:
            log_str = 'cpu(elapsed=%d)'%(t2-ResourceLogger.t1.pop())
        ResourceLogger.rlogger.info("%s%s %s"%(" "*ResourceLogger.print_res_spaces, func_name, log_str))

    @staticmethod
    def append_diff(stor, attr_name, attr_val):
        attr_val = int(attr_val)
        if ResourceLogger.print_nz_only and (attr_val == 0):
            return
        stor.append(attr_name+'='+str(attr_val))
        
    def __init__(self, call_ctx):
        self.call_ctx = call_ctx
    
    def __enter__(self):
        ResourceLogger.do_pre()
    
    def __exit__(self, *args):
        ResourceLogger.do_post(self.call_ctx)

try:
    import psutil
    ResourceLogger.psutil_available = True
except:
    ResourceLogger.psutil_available = False

try:
    import resource
    ResourceLogger.resource_available = True
except:
    ResourceLogger.resource_available = False


def print_rusage(func):
    """Prints CPU and memory usage for functions with print_resources decorator."""
    def wrapper(*args):
        if not ResourceLogger.is_enabled():
            return func(*args)

        ResourceLogger.do_pre()  
        res = func(*args)
        ResourceLogger.do_post(func.func_name)
        return res
    wrapper.__name__ = 'print_rusage_' + func.func_name
    return wrapper



#---------------------------------------------
# GC helper classes begin
#---------------------------------------------

class GCPreempt:
    enabled = False
    BEFORE = 1
    AFTER = 2
    logger = None
    WHEN_NAMES = ['', 'before', 'after']
    
    def __init__(self, ctx_name="", mode=3, indent=False):
        self.mode = mode
        self.ctx_name = ctx_name
        self.indent = indent

    def do_gc(self, mode, indent):
        if ResourceLogger.is_enabled():
            c1 = gc.get_count()
            gc.collect()
            c2 = gc.get_count()
            ResourceLogger.print_res_spaces +=  indent
            GCPreempt.logger.info('%scollected (%d, %d, %d) objects %s %s' % (" "*ResourceLogger.print_res_spaces, c1[0]-c2[0], c1[1]-c2[1], c1[2]-c2[2], GCPreempt.WHEN_NAMES[mode], self.ctx_name))
        else:
            gc.collect()

    def __enter__(self):
        if GCPreempt.enabled:
            if(self.mode & 1):
                self.do_gc(1, 2*int(self.indent))
    
    def __exit__(self, *args):
        if GCPreempt.enabled:
            if(self.mode & 2):
                self.do_gc(2, -2*int(self.indent))


def gc_wrap(func):
    """Preemptively call gc around functions that require memory and create garbage."""
    def wrapper(*args):
        with GCPreempt(func.func_name, 3, True):
            return func(*args)
    wrapper.__name__ = 'gc_wrap_' + func.func_name
    return wrapper

def gc_before(func):
    """Preemptively call gc before function that requires lots of memory."""
    def wrapper(*args):
        with GCPreempt(func.func_name, 1, True):
            return func(*args)
    wrapper.__name__ = 'gc_before_' + func.func_name
    return wrapper

def gc_after(func):
    """Preemptively call gc after function that creates heavy garbage."""
    def wrapper(*args):
        with GCPreempt(func.func_name, 2, True):
            return func(*args)
    wrapper.__name__ = 'gc_after_' + func.func_name
    return wrapper



class LowMemRetry:
    callback = None
    max_retry = 1
    retry_count = 0
    can_retry = True
    
    def __enter__(self):
        LowMemRetry.can_retry = True
        return self
    
    def __exit__(self, etype, value, traceback):
        if isinstance(value, MemoryError):
            LowMemRetry.is_memory_error = True
            LowMemRetry.retry_count += 1
            if (None != LowMemRetry.callback):
                LowMemRetry.callback(etype, value, traceback)

            if LowMemRetry.retry():
                return True
        LowMemRetry.can_retry = False
        return False

    @staticmethod
    def retry():
        if LowMemRetry.can_retry:
            return (LowMemRetry.retry_count <= LowMemRetry.max_retry)
        LowMemRetry.can_retry = True
        return False

def lowmem_retry(func):
    """On encountering MemoryError, retry if possible"""
    def wrapper(*args):
        while LowMemRetry.retry():
            with LowMemRetry():        
                return func(*args)
    wrapper.__name__ = 'lowmem_retry_' + func.func_name
    return wrapper
