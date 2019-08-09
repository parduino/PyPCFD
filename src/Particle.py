'''
Created on Jun 10, 2018

@author: pmackenz
'''

from numpy import array, ones, zeros, identity
import globalCounter as GC

class Particle(object):
    '''
    variables:
        self.id    = ParticleID
        self.mass  = mp
        self.pos   = xp
        self.vel   = vp
        self.accel = zeros(2)
        self.mu = 0.0;
        self.stress = zeros(3)
        self.strain = zeros(3)
        self.p      = 0.0
        self.strainRate = zeros(3)
        self.deformationGradient = identity(2)
    
    methods:
        def __init__(self, mp=1.0, xp=zeros(2), vp=zeros(2))
        def setViscosity(self, mu)
        def setVelocity(self, v)
        def addToVelocity(self, dv)
        def addToPosition(self, dx)
        def velocity(self)      # return particle velocity
        def position(self)      # return particle position
        def mass(self)          # return particle mass
        def strain(self)        # return particle strain
        def strainRate(self)    # return rate of deformation tensor
        def stress(self)        # return particle stress
    '''

    def __init__(self, mp=1.0, xp=zeros(2), vp=zeros(2)):
        '''
        Constructor
        '''
        # assign a unique ID
        GC.ParticleID += 1
        self.id    = GC.ParticleID
        
        self.mass  = mp
        self.pos   = xp
        self.vel   = vp
        self.accel = zeros(2)
        
        self.mu = 0.0;
        
        self.stress = zeros(3)
        self.strain = zeros(3)
        self.p      = 0.0
        self.strainRate = zeros(3)
        self.deformationGradient = identity(2)
        
    def setViscosity(self, mu):
        self.mu = mu;
        
    def addToVelocity(self, dv):
        self.vel += dv
        
    def setVelocity(self, v):
        self.vel = v
        
    def velocity(self):
        return self.vel.copy()
    
    def addToPosition(self, dx):
        self.pos += dx
        
    def position(self):
        return self.pos.copy()
    
    def mass(self):
        return self.mass
    
    def strain(self):
        return self.strain.copy()
    
    def strainRate(self):
        return self.strainRate.copy()
    
    def stress(self):
        stress = array([
                        2.*self.mu*self.strainRate[0] - self.p,
                        2.*self.mu*self.strainRate[1] - self.p,
                        self.mu*self.strainRate[2]
                        ])
        return stress

    def setDeformationGradient(self, newValue):
        self.deformationGradient = newValue

    def getDeformationGradient(self):
        return self.deformationGradient
        