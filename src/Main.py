'''
Created on Nov 21, 2015

@author: pmackenz
'''

from Domain import *

def Main():
    # defne the Reynolds number
    Re = 100
    
    # set sliding velocity
    velocity = 1.0
    
    # mass density of the fluid
    density = 1000.
    
    # set side-length of the analysis domain
    edgeDomain      = 1.
    # set the number of cells per edge
    numCellsPerEdge = 16
    
    # viscosity of the fluid
    viscosity = density * velocity * edgeDomain / Re
    
    # create an analysis domain
    domain = Domain(edgeDomain,edgeDomain,numCellsPerEdge,numCellsPerEdge)
    
    # configure the analysis type
    doInit         = False
    solveVstar     = True
    solveP         = True
    solveVtilde    = True
    solveVenhanced = True
    updatePosition = True 
    updateStress   = False
    addTransient   = True
    
    domain.setAnalysis(doInit, solveVstar, solveP, solveVtilde, solveVenhanced, updatePosition, updateStress, addTransient)
    domain.setParameters(Re, density, velocity)
    domain.setInitialState()
    
    CFL = 1.0
    dt = domain.getTimeStep(CFL)
    
    print(dt)
    
    #print(domain)
    
    time = 0.0
    
    dt = 0.5
    
    while (time+dt <= 10.0):
        time += dt
        domain.runAnalysis(time)
    
    dt = 5.0
    
    while (time+dt <= 100.0):
        time += dt
        domain.runAnalysis(time)
    
    

if __name__ == '__main__':
    Main()