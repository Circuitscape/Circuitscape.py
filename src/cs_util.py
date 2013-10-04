##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 
##

import os, string, sys
import time
from numpy import *

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


def elapsed_time(startTime): 
    """Returns elapsed time given a start time."""    
    now=time.time()
    elapsed=now-startTime
    secs=int(elapsed)
    mins=int(elapsed/60)
    hours=int(mins/60)
    mins=mins-hours*60
    secs=secs-mins*60-hours*3600
    return hours,mins,secs

def deleterow(A, delrow):
    m = A.shape[0]
    n = A.shape[1]
    keeprows = delete (arange(0, m), delrow)
    keepcols = arange(0, n)
    return A[keeprows][:,keepcols]
        
def deletecol(A, delcol):
    m = A.shape[0]
    n = A.shape[1]
    keeprows = arange(0, m)
    keepcols = delete (arange(0, n), delcol)
    return A[keeprows][:,keepcols]

def deleterowcol(A, delrow, delcol):
    m = A.shape[0]
    n = A.shape[1]

    keeprows = delete (arange(0, m), delrow)
    keepcols = delete (arange(0, n), delcol)

    return A[keeprows][:,keepcols]

def relabel(oldlabel, offset=0):
    newlabel = zeros(size(oldlabel), dtype='int32')
    s = sort(oldlabel)
    perm = argsort(oldlabel)
    f = where(diff(concatenate(([s[0]-1], s))))
    newlabel[f] = 1
    newlabel = cumsum(newlabel)
    newlabel[perm] = copy(newlabel)
    return newlabel-1+offset
