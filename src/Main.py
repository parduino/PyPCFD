'''
Created on Nov 21, 2015

@author: pmackenz
'''

from Domain import *
import subprocess
import ButcherTableau as integrator
from math import floor
from Mappings import *

def Main():
    # defne the Reynolds number
    Re = 1000
    Re = 1
    
    # set sliding velocity
    velocity = 1.0
    
    # mass density of the fluid
    density = 1000.
    
    # set side-length of the analysis domain
    edgeDomain      = 1.
    # set the number of cells per edge
    numCellsPerEdge = 16
    numCellsPerEdge = 8
    numCellsPerEdge = 4
    #numCellsPerEdge = 2
    
    # viscosity of the fluid
    viscosity = density * velocity * edgeDomain / Re
    
    # create an analysis domain
    #domain = Domain(edgeDomain, edgeDomain, numCellsPerEdge, numCellsPerEdge, mappingFunction=IdentityMap())
    domain = Domain(edgeDomain, edgeDomain, numCellsPerEdge, numCellsPerEdge, mappingFunction=FineEdgeMap())

    domain.createParticles(2,2)
    
    # configure the analysis type
    doInit         = False
    solveVstar     = True
    solveP         = True
    solveVtilde    = True
    solveVenhanced = False
    updatePosition = True
    updateStress   = False
    addTransient   = True
    
    domain.setAnalysis(doInit,
                       solveVstar,
                       solveP,
                       solveVtilde,
                       solveVenhanced,
                       updatePosition,
                       updateStress,
                       addTransient)
    domain.setParameters(Re, density, velocity)
    domain.setInitialState()
    
    CFL = 1.0
    dt = domain.getTimeStep(CFL)
    
    print(u"CFL=1 equals to \u0394t={:f}".format(dt))


    print(domain)

    # define load history and print interval

    dt1 = 0.01000
    target1 = 0.1

    dt2 = 0.1
    target2 = 1.0

# ************* don't mess with stuff below *************

    domain.particleTrace(True)

    # defining plot settings
    domain.setPlotInterval(dt1)

    # defining output settings
    domain.setWriteInterval(-1)

    # initializing starting time
    time = 0.0
    
    # run first segment
    #domain.setTimeIntegrator(integrator.ExplicitEuler())
    domain.setTimeIntegrator(integrator.RungeKutta4())

    domain.plotParticleTrace('tracePlot{:04d}.png'.format(floor(time*100)))

    dt = dt1
    while (time+dt <= target1 + 0.1 * dt):
        time += dt
        domain.runAnalysis(time)

    domain.plotParticleTrace('tracePlot{:04d}.png'.format(floor(time*100)))

    # run second segment
    domain.setTimeIntegrator(integrator.RungeKutta4())

    dt = dt2
    while (time + dt <= target2 + 0.1 * dt):
        time += dt
        domain.runAnalysis(time)

        if (time % 1.0 < 0.5*dt):   # write 1.0 sec duration trace plots
            domain.plotParticleTrace('tracePlot{:04d}.png'.format(floor(time*100)))
            domain.particleTrace(False)   # this wipes old trace
            domain.particleTrace(True)    # this restarts trace

    
    # generate the animation
    subprocess.run('./makeAnim.sh')
    

if __name__ == '__main__':
    Main()
