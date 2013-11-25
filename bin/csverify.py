#!/usr/bin/python
##
## Circuitscape (C) 2013, Brad McRae, Viral B. Shah. and Tanmay Mohapatra

from circuitscape.verify import cs_verifyall
import sys

if __name__ == '__main__':
    out_path = None
    root_path = 'circuitscape'
    
    if len(sys.argv) > 1:
        if len(sys.argv) == 3:
            out_path = sys.argv.pop()
        elif len(sys.argv) > 3:
            raise RuntimeError("invalid number of options")
            
        root_path = sys.argv.pop()

    cs_verifyall(root_path, out_path)
