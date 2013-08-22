#!/usr/bin/python
##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: cs_run.py 804 2012-07-30 23:05:05Z mcrae $
##


import sys
from circuitscape import *

if len(sys.argv) == 1:
    print 'Error: Circuitscape configuration (.ini) file required.'
else:
    configFile = sys.argv[1]
    cs = circuitscape(configFile, 'Screen')
    resistances = cs.compute()
    print resistances

