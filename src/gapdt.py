##
## Circuitscape (C) 2008, Brad McRae and Viral B. Shah. 
##
## $Id: gapdt.py 740 2011-04-13 02:50:14Z viral $
##

from numpy import *
from scipy import sparse

class gapdt:

    def relabel(self, oldlabel, offset=0):
        newlabel = zeros(size(oldlabel), dtype='int32')
        s = sort(oldlabel)
        perm = argsort(oldlabel)
        f = where(diff(concatenate(([s[0]-1], s))))
        newlabel[f] = 1
        newlabel = cumsum(newlabel)
        newlabel[perm] = copy(newlabel)
        return newlabel-1+offset

    def subsref(self, A, I, J):
        B = A[:, J][I, :]
        
        return B

    def deleterowcol(self, A, delrow, delcol):
        m = A.shape[0]
        n = A.shape[1]

        keeprows = delete (arange(0, m), delrow)
        keepcols = delete (arange(0, n), delcol)

        return A[keeprows][:,keepcols]
