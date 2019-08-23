import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['font.size'] = 15

from math import log, floor
from numpy import array
from numpy.linalg import norm
import os

from Domain import *
from Motion import *
from ButcherTableau import *


class LocalConvergenceTest(object):

    def __init__(self, motion, algorithm, fileType='png', nCells=1):
        self.numAlgorithms = (algorithm,)
        self.motionList = (motion,)
        self.fileType = fileType
        if nCells >= 1.:
            self.nCells= floor(nCells)
        else:
            self.nCells = 1

    def runAnalysis(self):
        for motion in self.motionList:
            # create POSITION plots for each motion
            # fig, ax1 = plt.subplots()
            for j, numalg in enumerate(self.numAlgorithms):
                dtList, positionErrors, Ferrors = self.runCase(numalg, motion)

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

        Ferrors = []
        positionErrors = []
        dtList = []
        dt = 1.0
        while (dt > 1.0e-10):
            dtList.append(dt)
            domain = Domain(width=1., height=1., nCellsX=self.nCells, nCellsY=self.nCells)
            domain.setMotion(motion)
            domain.setTimeIntegrator(numAlg)

            domain.setAnalysis(doInit, solveVstar, solveP,
                               solveVtilde, solveVenhanced,
                               updatePosition, updateStress,
                               addTransient)

            domain.setPlotInterval(dt)        # plot at the end of each time step
            domain.setWriteInterval(-1)       # no recorder output

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

            mask = 'dt = {:.2E}, Position error = {:.3E}, F error = {:.3E}'
            print(mask.format(dt, positionErrors[-1], Ferrors[-1]))
            if (Ferrors[-1] < 1.e-15 or positionErrors[-1] < 1.e-15):
                break
            # dt /= 10.
            dt /= 2.  # allows for a more precise identification of numeric truncation error

        # create folder to store images
        if not os.path.isdir("images"):
            os.mkdir("images")

        self.plotPositionErrors(dtList, positionErrors, numAlg, motion)

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

        plt.close()

        slope = log(Ferrors[0] / Ferrors[-1]) / log(dtList[0] / dtList[-1])
        print('{} {} Deformation Gradient convergence slope = {:.2f}'.format(numAlg, motion, slope))

        return (dtList, positionErrors, Ferrors)

    def plotPositionErrors(self, dtList, positionErrors, numAlg, motion):
        # Plots for position errors
        fig, ax2 = plt.subplots()
        matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
        matplotlib.rcParams['font.size'] = 15
        ax2 = fig.gca()
        ax2.loglog(dtList, positionErrors, 'k-o', linewidth=2, label="simulation")

        x = array([dtList[0], dtList[-1]])
        y1 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (1.)])
        y2 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (2.)])
        y3 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (3.)])
        y4 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (4.)])
        y5 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (5.)])

        ax2.loglog(x, y1, 'y--', linewidth=2, label="1st order")
        ax2.loglog(x, y2, 'b:',  linewidth=2, label="2nd order")
        ax2.loglog(x, y3, 'g-.', linewidth=2, label="3rd order")
        ax2.loglog(x, y4, 'm--', linewidth=2, label="4th order")
        ax2.loglog(x, y5, 'r--', linewidth=2, label="5th order")

        ax2.set_xlabel('$\Delta t$ (s)')
        ax2.set_ylabel('$|| x_{numerical} - x_{analytical} ||_{2}$')
        ax2.legend(loc="best")
        ax2.set_ylim([1e-17, 1e1])

        ax2.grid(True)
        fileName = "{}_{}_Position_convergence.{}".format(numAlg, motion, self.fileType)
        fileNameWithPath = os.path.join("images", fileName)

        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

        plt.close()

        slope = log(positionErrors[0] / positionErrors[-1]) / log(dtList[0] / dtList[-1])
        print('{} {} Position convergence slope = {:.2f}'.format(numAlg, motion, slope))

    def finalizePositionPlot(self, ax, dtList, positionErrors):
        x = array([dtList[0], dtList[-1]])
        y1 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (1.)])
        y2 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (3.)])
        y3 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (3.)])
        y4 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (4.)])
        y5 = array([positionErrors[0], positionErrors[0] * (dtList[-1] / dtList[0]) ** (5.)])

        ax.loglog(x, y1, 'y--', linewidth=2, label="1st order")
        ax.loglog(x, y2, 'b:',  linewidth=2, label="2nd order")
        ax.loglog(x, y3, 'g-.', linewidth=2, label="3rd order")
        ax.loglog(x, y4, 'm--', linewidth=2, label="4th order")
        ax.loglog(x, y5, 'r--', linewidth=2, label="5th order")

        ax.set_xlabel('$\Delta t$ (s)')
        ax.set_ylabel('$|| x_{numerical} - x_{analytical} ||_{2}$')
        ax.legend(loc="best")
        ax.set_ylim([1e-17, 1e1])

        ax.grid(True)