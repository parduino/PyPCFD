'''
Created on Nov 21, 2015

@author: pmackenz
'''

from Domain import *

def Main():
    # defne the Reynolds number
    Re = 1000
    
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
    
    domain.runAnalysis(0.5)
    domain.runAnalysis(1.0)
    domain.runAnalysis(1.5)
    domain.runAnalysis(2.0)
    domain.runAnalysis(2.5)
    domain.runAnalysis(3.0)
    domain.runAnalysis(3.5)
    domain.runAnalysis(4.0)
    domain.runAnalysis(4.5)
    domain.runAnalysis(5.0)
    domain.runAnalysis(5.5)
    domain.runAnalysis(6.0)
    domain.runAnalysis(6.5)
    domain.runAnalysis(7.0)
    domain.runAnalysis(7.5)
    domain.runAnalysis(8.0)
    domain.runAnalysis(8.5)
    domain.runAnalysis(9.0)
    domain.runAnalysis(9.5)
    domain.runAnalysis(10.0)
    domain.runAnalysis(15.0)
    domain.runAnalysis(20.0)
    domain.runAnalysis(25.0)
    domain.runAnalysis(30.0)
    domain.runAnalysis(35.0)
    domain.runAnalysis(40.0)
    domain.runAnalysis(45.0)
    domain.runAnalysis(50.0)
    domain.runAnalysis(55.0)
    domain.runAnalysis(60.0)
    domain.runAnalysis(65.0)
    domain.runAnalysis(70.0)
    domain.runAnalysis(75.0)
    domain.runAnalysis(80.0)
    domain.runAnalysis(85.0)
    domain.runAnalysis(90.0)
    domain.runAnalysis(95.0)
    domain.runAnalysis(100.0)
    
    
    

if __name__ == '__main__':
    Main()