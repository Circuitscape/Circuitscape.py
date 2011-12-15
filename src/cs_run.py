#!/usr/bin/python
##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: cs_run.py 440 2009-03-13 01:05:07Z viral $
##

import sys
from cs_compute import *

configFile = sys.argv[1]
cs = cs_compute(configFile, None)
resistances = cs.compute()
print resistances
