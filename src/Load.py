'''
Created on Nov 21, 2015

@author: pmackenz
'''

from Vector import *

class Load(object):
    '''
    classdocs
    '''

    def __init__(self, force=[0,0,0], id=-1):
        '''
        Constructor
        '''
        self.force = Vector(force)
        self.node = id
        
    def atNode(self,ID):
        self.node = ID
        
    def getNodeID(self):
        return self.node
    
    def getForce(self):
        return self.force
    