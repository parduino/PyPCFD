from Domain import *
import matplotlib.pyplot as plt
import matplotlib
from math import log
import os
from Motion import *
from ButcherTableau import *

def runAnalysis(numAlg, motion):
    # configure the analysis type
    doInit = False
    solveVstar = False
    solveP = False
    solveVtilde = False
    solveVenhanced = False
    updatePosition = True
    updateStress = False
    addTransient = False
    plotFigures = True
    writeOutput = False

    Ferrors = []
    positionErrors = []
    dtList = []
    dt = 1.0
    while (dt > 1.0e-8):
        dtList.append(dt)
        domain = Domain(width=1., height=1., nCellsX=1, nCellsY=1,
                        motion=motion,
                        particleUpdateScheme=numAlg)

        domain.setAnalysis(doInit, solveVstar, solveP,
                           solveVtilde, solveVenhanced,
                           updatePosition, updateStress,
                           addTransient, plotFigures, writeOutput)

        # you need to set the velocity field to the initial velocity field
        # (or to any fixed time throughout the test !!!)
        #domain.setState(dt)    # this line is wrong
        domain.setState(0)

        # returns error for all particles. Use [0] for first [1] for second etc..
        FerrorList, positionErrorList = domain.updateParticleMotion(dt)

        Ferrors.append(FerrorList[0])
        positionErrors.append(positionErrorList[0])

        print('dt = {:.2E}, Position error = {:.3E}, F error = {:.3E}'.format(dt, positionErrors[-1], Ferrors[-1]))
        # print('dt = {:.2E}, Position error = {:.3E}'.format(dt, positionErrors[-1]))
        # print("\n")
        if (Ferrors[-1] < 1.e-16):
            break
        dt /= 2.0

    # create folder to store images
    if not os.path.isdir("images"):
        os.mkdir("images")

    # Plots for deformation gradient errors
    fig = plt.figure()
    matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
    matplotlib.rcParams['font.size'] = 15
    ax1 = fig.gca()
    ax1.loglog(dtList, Ferrors, 'k-o', linewidth=2, label="simulation")

    x = array([dtList[0], dtList[-1]])
    y = array([Ferrors[0], Ferrors[0] * (dtList[-1] / dtList[0]) ** (1.)])
    ax1.loglog(x, y, 'y--', linewidth=2, label="1st order")

    x = array([dtList[0], dtList[-1]])
    y = array([Ferrors[0], Ferrors[0] * (dtList[-1] / dtList[0]) ** (2.)])
    ax1.loglog(x, y, 'b:', linewidth=2, label="2nd order")

    x = array([dtList[0], dtList[-1]])
    y = array([Ferrors[0], Ferrors[0] * (dtList[-1] / dtList[0]) ** (3.)])
    ax1.loglog(x, y, 'g-.', linewidth=2, label="3rd order")

    x = array([dtList[0], dtList[-1]])
    y = array([Ferrors[0], Ferrors[0] * (dtList[-1] / dtList[0]) ** (4.)])
    ax1.loglog(x, y, 'm--', linewidth=2, label="4th order")

    x = array([dtList[0], dtList[-1]])
    y = array([Ferrors[0], Ferrors[0] * (dtList[-1] / dtList[0]) ** (5.)])
    ax1.loglog(x, y, 'r--', linewidth=2, label="5th order")

    # ax1.set_ylim(1e-16, 1e1)
    ax1.set_xlabel('$\Delta t$ (s)')
    ax1.set_ylabel('$|| F_{numerical} - F_{analytical} ||_{2}$')

    plt.ylim([1e-16,1.])

    ax1.legend(loc="best")

    ax1.grid(True)

    # ax1.axis('tight')
    plt.savefig(os.path.join("images", "{}_{}_F_convergence.pdf".format(numAlg, motion)), pad_inches=0,
                bbox_inches='tight')
    plt.savefig(os.path.join("images", "{}_{}_F_convergence.png".format(numAlg, motion)), pad_inches=0,
                bbox_inches='tight')

    slope = log(Ferrors[0] / Ferrors[-1]) / log(dtList[0] / dtList[-1])
    print('{} {} Deformation Gradient convergence slope = {:.2f}'.format(numAlg, motion, slope))

    # Plots for position errors
    fig = plt.figure()
    matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
    matplotlib.rcParams['font.size'] = 15
    ax2 = fig.gca()
    ax2.loglog(dtList, positionErrors, 'k-o', linewidth=2, label="simulation")

    x = array([dtList[0], dtList[-1]])
    y = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (1.)])
    ax2.loglog(x, y, 'y--', linewidth=2, label="1st order")

    x = array([dtList[0], dtList[-1]])
    y = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (2.)])
    ax2.loglog(x, y, 'b:', linewidth=2, label="2nd order")

    x = array([dtList[0], dtList[-1]])
    y = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (3.)])
    ax2.loglog(x, y, 'g-.', linewidth=2, label="3rd order")

    x = array([dtList[0], dtList[-1]])
    y = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (4.)])
    ax2.loglog(x, y, 'm--', linewidth=2, label="4th order")

    x = array([dtList[0], dtList[-1]])
    y = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (5.)])
    ax2.loglog(x, y, 'r--', linewidth=2, label="5th order")

    # ax2.set_ylim(1e-16, 1e3)
    ax2.set_xlabel('$\Delta t$ (s)')
    ax2.set_ylabel('$|| x_{numerical} - x_{analytical} ||_{2}$')

    ax2.legend(loc="best")

    plt.ylim([1e-16,1.])

    ax2.grid(True)
    # ax2.axis('tight')
    plt.savefig(os.path.join("images", "{}_{}_Position_convergence.pdf".format(numAlg, motion)), pad_inches=0,
                bbox_inches='tight')
    plt.savefig(os.path.join("images", "{}_{}_Position_convergence.png".format(numAlg, motion)), pad_inches=0,
                bbox_inches='tight')

    slope = log(positionErrors[0] / positionErrors[-1]) / log(dtList[0] / dtList[-1])
    print('{} {} Position convergence slope = {:.2f}'.format(numAlg, motion, slope))
    # print(positionErrors)
    # print(dtList)



def Main():
    numAlgorithms = {"ee": ExplicitEuler(), "rk4": RungeKutta4(), "midpt":MidPointRule(), "heun": HeunsMethod()}
    motionDict = {"m1": Motion1(), "m2": Motion2()}

    # filenames = []
    # for algName, numalg in numAlgorithms.items():
    #     for mName, motion in motionDict.items():
    #         runAnalysis(numalg, motion)
    #         filenames.append("{}_{}_Position_convergence.pdf".format(numalg, motion))
    #         filenames.append("{}_{}_F_convergence.pdf".format(numalg, motion))
    #         print("\n")
    #
    # for fn in filenames:
    #     print(fn)

    #runAnalysis(RungeKutta4(), Motion1())
    runAnalysis(numAlgorithms['heun'], motionDict['m1'])

if __name__ == '__main__':
    Main()