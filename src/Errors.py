'''
Created on Jun 26, 2018

@author: pmackenz
'''

class CellIndexError(Exception):
    '''
    classdocs
    '''

    def __init__(self, e):
        '''
        Constructor
        '''
        super().__init__()
        self.i = e[0]
        self.j = e[1]
        self.k = e[2]
        self.pos = e[3]
        
    def __str__(self):
        return "position {}/{} resulted in i={}, j={} and k={}".format(*self.pos, self.i, self.j, self.k)
    