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
    
    print(u"CFL=1 equals to \u0394t={:f}".format(dt))
    
    #print(domain)
    
    # define load history and print interval
    
    dt1 = 0.5
    target1 = 10.0
    
    dt2 = 0.5
    target2 = 10.0

# ************* don't mess with stuff below *************

    time = 0.0
    
    dt = dt1
    
    while (time+dt <= target2):
        time += dt
        domain.runAnalysis(time)
    
    dt = dt2
    
    while (time+dt <= target2):
        time += dt
        domain.runAnalysis(time)
    
    

if __name__ == '__main__':
    Main()