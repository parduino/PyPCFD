from Domain import *
import matplotlib.pyplot as plt
import matplotlib
from math import log
import os
from Motion import *
from ButcherTableau import *


class LocalConvergenceTest(object):

    def __init__(self, fileType='png'):
        self.numAlgorithms = (ExplicitEuler(), RungeKutta4(), MidPointRule(), MidPointMethod2())
        self.motionList = (Motion1(), Motion2())
        self.fileType = fileType

    def runAnalysis(self):
        for numalg in self.numAlgorithms:
            for motion in self.motionList:
                self.runCase(numalg, motion)
                print()

    def runCase(self, numAlg, motion):
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
            # domain.setState(dt)    # this line is wrong
            domain.setState(0)
            x0 = domain.getParticles()[0].position() # save original position of particle

            # update particle
            domain.updateParticleMotion(dt)
            # calculate and store errors from updated particle position
            particle = domain.getParticles()[0]
            posError = norm( particle.position() - motion.getAnalyticalPosition(x0, dt) )
            FError = norm(particle.getDeformationGradient() - motion.getAnalyticalF(x0, dt) )
            Ferrors.append(FError)
            positionErrors.append(posError)

            print('dt = {:.2E}, Position error = {:.3E}, F error = {:.3E}'.format(dt, positionErrors[-1], Ferrors[-1]))
            if (Ferrors[-1] < 1.e-14 or positionErrors[-1] < 1.e-14):
                break
            dt /= 10.

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

        ax1.set_ylim(1e-17, 1e2)
        ax1.set_xlabel('$\Delta t$ (s)')
        ax1.set_ylabel('$|| F_{numerical} - F_{analytical} ||_{2}$')

        plt.ylim([1e-17, 1e2])

        ax1.legend(loc="best")

        ax1.grid(True)

        fileName = "{}_{}_F_convergence.{}".format(numAlg, motion, self.fileType)
        fileNameWithPath = os.path.join("images", fileName)

        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

        slope = log(Ferrors[0] / Ferrors[3]) / log(dtList[0] / dtList[-1])
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

        ax2.set_ylim([1e-17, 1e1])

        ax2.grid(True)
        # ax2.axis('tight')
        fileName = "{}_{}_Position_convergence.{}".format(numAlg, motion, self.fileType)
        fileNameWithPath = os.path.join("images", fileName)

        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

        plt.close()

        slope = log(positionErrors[0] / positionErrors[3]) / log(dtList[0] / dtList[-1])
        print('{} {} Position convergence slope = {:.2f}'.format(numAlg, motion, slope))
        # print(positionErrors)
        # print(dtList)