#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae, Viral B. Shah. and Tanmay Mohapatra 

import sys
from circuitscape.compute import Compute

if len(sys.argv) == 1:
    print 'Error: Circuitscape configuration (.ini) file required.'
else:
    configFile = sys.argv[1]
    cs = Compute(configFile, 'Screen')
    resistances = cs.compute()
    print resistances

