#!/usr/bin/python
##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 

import sys
from numpy import *
from circuitscape import *

nrows = int(sys.argv[1])
ncols = int(sys.argv[2])
f = open('./verify/stress/cellmap.asc', 'w')
f.write('ncols         ' + str(ncols) + '\n')
f.write('nrows         ' + str(nrows) + '\n')
f.write('xllcorner     ' + str(0) + '\n')
f.write('yllcorner     ' + str(0) + '\n')
f.write('cellsize      ' + str(1) + '\n')
f.write('NODATA_value  ' + str(9999) + '\n')

delimiter = ''
fmt = ['%.6f ']*ncols
format = delimiter.join(fmt)
data = ones((nrows, ncols), dtype='float64')
for row in data:
    f.write(format % tuple(row) + '\n')
f.close()
    
configFile = './verify/stress/stress.ini'
cs = circuitscape(configFile, 'Screen')
resistances = cs.compute()
print resistances
