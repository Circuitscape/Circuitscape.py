#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae and Viral B. Shah. 

import sys
from circuitscape import *

if len(sys.argv) == 1:
    print 'Error: Circuitscape configuration (.ini) file required.'
else:
    configFile = sys.argv[1]
    cs = circuitscape(configFile, 'Screen')
    resistances = cs.compute()
    print resistances

