'''
Created on Jun 26, 2018

@author: pmackenz
'''
from numpy import array
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix

class matrixDataType(object):
    '''
    variables:
        self.smat
    
    methods:
        def __init__(self, data=[], idxi=[], idxj=[])
        def __str__(self, *args, **kwargs)
        def __repr__(self, *args, **kwargs)
        def add(self, val, i, j)
        def toCSCmatrix(self)
        def toCSRmatrix(self)
    '''
    
    def __init__(self, ndof):
        self.smat = dok_matrix((ndof,ndof))
        
    def __str__(self, *args, **kwargs):
        return str(self.smat)
    
    def __repr__(self, *args, **kwargs):
        return repr(self.smat)
        
    def add(self, val, i, j):
        self.smat[i,j] += val
        
    def toCSCmatrix(self):
        return csc_matrix(self.smat, copy=False)
        
    def toCSRmatrix(self):
        return csr_matrix(self.smat, copy=False)
        

