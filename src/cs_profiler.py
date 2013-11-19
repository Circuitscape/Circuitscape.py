import logging, time, os

class ResourceLogger:
    print_timings = False
    print_rusages = False    
    print_res_spaces = 0
    psutil_available = True
    resource_available = True

    rlogger = None
    proc = None
    proc_has_io_counters = False
    
    t1 = []
    mem1 = []
    io1 = []

    @staticmethod
    def init_rusage(print_t=False, print_r=False, logger=None):
        ResourceLogger.print_timing_enabled(print_t)
        ResourceLogger.print_rusage_enabled(print_r)
        ResourceLogger.rlogger = logger
        
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
            ResourceLogger.rusage1 = resource.getrusage(resource.RUSAGE_SELF)
            if ResourceLogger.proc:
                ResourceLogger.mem1.append(ResourceLogger.proc.get_ext_memory_info())
                ResourceLogger.io1.append(ResourceLogger.proc.get_io_counters() if ResourceLogger.proc_has_io_counters else psutil.disk_io_counters())
    
    @staticmethod
    def do_post(func_name):
        t2 = time.time()
        ResourceLogger.print_res_spaces -=  2
        
        if ResourceLogger.print_rusages:
            rusage = resource.getrusage(resource.RUSAGE_SELF)
            utime = rusage.ru_utime - ResourceLogger.rusage1.ru_utime
            stime = rusage.ru_stime - ResourceLogger.rusage1.ru_stime
            cpu_str = 'cpu(user=%d, system=%d, elapsed=%d)'%(utime,stime, (t2-ResourceLogger.t1.pop()))
            mem_diffs = []
            io_diffs = []
            if ResourceLogger.proc:
                mem1 = ResourceLogger.mem1.pop()
                io1 = ResourceLogger.io1.pop()
                mem2 = ResourceLogger.proc.get_ext_memory_info()
                io2 = ResourceLogger.proc.get_io_counters() if ResourceLogger.proc_has_io_counters else psutil.disk_io_counters()
                for attr in ['rss', 'vms', 'shared', 'text', 'data']:
                    if hasattr(mem1, attr):
                        mem_diffs.append(attr+'='+str(getattr(mem2, attr) - getattr(mem1, attr)))
                for attr in ['read_count', 'write_count', 'read_bytes', 'write_bytes']:
                    if hasattr(io1, attr):
                        io_diffs.append(attr+'='+str(getattr(io2, attr) - getattr(io1, attr)))
            else:
                mem_diffs.append('maxrss=' + str(rusage.ru_maxrss - ResourceLogger.rusage1.ru_maxrss))
                mem_diffs.append('shared=' + str(rusage.ru_ixrss - ResourceLogger.rusage1.ru_ixrss))
                mem_diffs.append('heap=' + str(rusage.ru_idrss - ResourceLogger.rusage1.ru_idrss))
                mem_diffs.append('stack=' + str(rusage.ru_isrss - ResourceLogger.rusage1.ru_isrss))
                io_diffs.append('read_count=' + str(rusage.ru_inblock - ResourceLogger.rusage1.ru_inblock))
                io_diffs.append('write_count=' + str(rusage.ru_oublock - ResourceLogger.rusage1.ru_oublock))
            mem_str = 'mem(' + ' '.join(mem_diffs) + ')'
            io_str = 'io(' + ' '.join(io_diffs) + ')'
            log_str = ' '.join([cpu_str, mem_str, io_str])
        else:
            log_str = 'cpu(elapsed=%d)'%(t2-ResourceLogger.t1.pop())
        ResourceLogger.rlogger.info("%s%s %s"%(" "*ResourceLogger.print_res_spaces, func_name, log_str))
        
try:
    import psutil
except:
    ResourceLogger.psutil_available = False

try:
    import resource
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
    return wrapper


