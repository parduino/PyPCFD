'''
Created on Jun 13, 2018
Modified on June 24, 2018 for separate images

@author: pmackenz
'''
import numpy as np
import os


class Writer(object):
    '''
    variables:
        self.particlesPresent     # boolean
        self.height
        self.width
        self.nNodesX = nCellsX+1
        self.nNodesY = nCellsY+1
        self.Y
        self.X
        self.Vx
        self.Vy
        self.speed
        self.gs
        self.tracerPoints = [[],[]]
    
    methods:
        def __init__(self)
        
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.DATA_COUNTER  = -1
        self.particlesPresent = False
        
        self.width  = 1.0
        self.height = 1.0

        if not os.path.isdir("data"):
            os.mkdir("data")
        
    def setGrid(self, width, height, nCellsX, nCellsY):
        self.height = height
        self.width  = width
        self.nNodesX = nCellsX+1
        self.nNodesY = nCellsY+1
        
        x = np.linspace(0,width ,(nCellsX+1))
        y = np.linspace(0,height,(nCellsY+1))
        
        self.X, self.Y = np.meshgrid(x, y)
        self.Vx = np.zeros_like(self.X)
        self.Vy = np.zeros_like(self.X)
        
        # define tracer points
        
        xp = self.width*(0.5+np.arange(nCellsX)) / nCellsX
        yp = self.height*(0.5+np.arange(nCellsY))/ nCellsY
        refPtsX, refPtsY = np.meshgrid(xp,yp)
        
        self.tracerPoints = [refPtsX,refPtsY]
        
    def setData(self, nodes):
        self.particlesPresent = False
        
        self.P  = np.zeros_like(self.X)
        self.Vx = np.zeros_like(self.X)
        self.Vy = np.zeros_like(self.X)
        self.Fx = np.zeros_like(self.X)
        self.Fy = np.zeros_like(self.X)
        
        for i in range(self.nNodesX):
            for j in range(self.nNodesY):
                node = nodes[i][j]
                self.P[j,i] = node.getPressure()
                vel = node.getVelocity()
                self.Vx[j,i] = vel[0]
                self.Vy[j,i] = vel[1]
                force = node.getForce()
                self.Fx[j,i] = force[0]
                self.Fy[j,i] = force[1]
        
        self.speed = np.sqrt(self.Vx*self.Vx + self.Vy*self.Vy)
        
    def setParticleData(self, particles):
        if (len(particles) == 0):
            return
        
        self.particlesPresent = True
        
        self.ParticleX  = []
        self.ParticleY  = []
        self.ParticleVx = []
        self.ParticleVy = []
        
        for p in particles:
            pos = p.position()
            vel = p.velocity()
            self.ParticleX.append(pos[0])
            self.ParticleY.append(pos[1])
            self.ParticleVx.append(vel[0])
            self.ParticleVy.append(vel[1])
            
    def writeData(self, time):

        # store nodal co-ordinates only once at t=0
        if self.DATA_COUNTER == -1:
            fname = os.path.join('data','nodeXcoordinates.txt')
            np.savetxt(fname, self.X)

            fname = os.path.join('data', 'nodeYcoordinates.txt')
            np.savetxt(fname, self.Y)

        # store velocities and pressures at each time-step
        self.DATA_COUNTER += 1
        hdr = 't={:08.8f}s'.format(time)
        fname = os.path.join('data',"vx{:03d}.txt".format(self.DATA_COUNTER))
        np.savetxt(fname, self.Vx, header=hdr)

        fname = os.path.join('data', "vy{:03d}.txt".format(self.DATA_COUNTER))
        np.savetxt(fname, self.Vy, header=hdr)

        fname = os.path.join('data', "pressure{:03d}.txt".format(self.DATA_COUNTER))
        np.savetxt(fname, self.P, header=hdr)


