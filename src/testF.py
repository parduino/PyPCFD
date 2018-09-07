from Domain import *
import matplotlib.pyplot as plt
from math import log


def Main():
    # configure the analysis type
    doInit         = False
    solveVstar     = False
    solveP         = False
    solveVtilde    = False
    solveVenhanced = False
    updatePosition = True
    updateStress   = False
    addTransient   = False
    plotFigures    = True
    writeOutput    = False

    errorList = []
    dtList = []
    dt = 1.0
    while ( dt>1.0e-10 ):
        dtList.append(dt)
        domain = Domain(width=2, height=2, nCellsX=8, nCellsY=8)
        domain.setAnalysis(doInit, solveVstar, solveP, solveVtilde, solveVenhanced, updatePosition, updateStress, addTransient, plotFigures, writeOutput)
        domain.setState(0)
        # domain.plotData()
        domain.setState(dt=dt)
        errorList.append(domain.updateParticleMotion(dt)[5]) # returns error for all particles. Use [0] for first [1] for second etc..
        print('dt = {:.2E}, error = {:.3E}'.format(dt, errorList[-1]))
        dt /= 10.0
        # domain.plotData()

    plt.figure(1)
    plt.loglog(dtList, errorList, 'k-')
    plt.scatter(dtList, errorList)
    plt.xlabel('$\Delta t$ (s)')
    plt.ylabel('$|| F_{numerical} - F_{analytical} ||_{2}$')
    plt.grid(True)
    plt.axis('equal')
    plt.savefig("convergence.pdf", bbox_inches='tight')

    slope = log(errorList[0] / errorList[-1]) / log(dtList[0] / dtList[-1])
    print('slope = {:.2E}'.format(slope))


if __name__ == '__main__':
    Main()