#!/usr/bin/python
##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: cs_run.py 795 2012-01-04 18:14:31Z mcrae $
##


import sys
from cs_compute import *
if len(sys.argv) == 1:
    print 'Error: Circuitscape configuration (.ini) file required.'
else:
    configFile = sys.argv[1]
    print 'Calling Circuitscape...'
    cs = cs_compute(configFile, 'Screen')
    resistances = cs.compute()
    print resistances
    print '\nCircuitscape calculations completed.'
    
