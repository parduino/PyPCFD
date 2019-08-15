from Domain import *
import matplotlib
import matplotlib.pyplot as plt
from math import log
import os
from Motion import *
from ButcherTableau import *


class GlobalConvergenceTest(object):

    def __init__(self, fileType='png'):
        self.numAlgorithms = (ExplicitEuler(), MidPointRule(), RungeKutta4())
        self.motionList = (Motion1(), Motion2())
        # self.numAlgorithms = (MidPointRule(),)
        # self.motionList = (Motion2(),)
        self.fileType = fileType

        # configure the analysis type
        self.doInit = False
        self.solveVstar = False
        self.solveP = False
        self.solveVtilde = False
        self.solveVenhanced = False
        self.updatePosition = True
        self.updateStress = False
        self.addTransient = False

        # Domain parameters
        self.domainHeight = 1.0
        self.domainWidth = 1.0

        # create folder to store images
        if not os.path.isdir("images"):
            os.mkdir("images")

    def runAnalysis(self):
        for numalg in self.numAlgorithms:
            for motion in self.motionList:
                self.runCase(numalg, motion)

    def runCase(self, numAlg, motion):
        maxTime = 1.0
        Ferrors = []
        positionErrors = []
        NList = []
        N = 1
        dt = maxTime/N
        while (N<=1e3):
            NList.append(N)

            domain = Domain(width=1., height=1., nCellsX=1, nCellsY=1)
            domain.setMotion(motion)
            domain.setTimeIntegrator(numAlg)

            domain.setAnalysis(self.doInit, self.solveVstar, self.solveP,
                               self.solveVtilde, self.solveVenhanced,
                               self.updatePosition, self.updateStress,
                               self.addTransient)

            domain.setPlotInterval(maxTime)   # plot only at the end
            domain.setWriteInterval(-1)       # no recorder output

            # Set the velocity field to the initial velocity field
            x0 = domain.getParticles()[0].position()  # save original position of particle for comparison later

            # update particle
            for j in range(N):
                domain.setState(j*dt)
                domain.updateParticleMotion(dt)

            # calculate and store errors from updated particle position
            particle = domain.getParticles()[0]
            posError = norm( particle.position() - motion.getAnalyticalPosition(x0, dt*N))
            FError = norm(particle.getDeformationGradient() - motion.getAnalyticalF(x0, dt*N))
            Ferrors.append(FError)
            positionErrors.append(posError)

            mask = 'N = {}, dt = {:.2E}, Position error = {:.3E}, F error = {:.3E}'
            print(mask.format(N, dt, positionErrors[-1], Ferrors[-1]))
            if (Ferrors[-1] < 1.e-14 or positionErrors[-1] < 1.e-14):
                break
            N *= 10
            dt = maxTime/N

        # Plots for deformation gradient errors
        fig = plt.figure()
        matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
        matplotlib.rcParams['font.size'] = 15
        ax1 = fig.gca()
        ax1.loglog(NList, Ferrors, 'k-o', linewidth=2, label="simulation")

        x = array([NList[0], NList[-1]])
        y = array([Ferrors[0], Ferrors[0] * (NList[0] / NList[-1]) ** (1.)])
        ax1.loglog(x, y, 'y--', linewidth=2, label="1st order")

        x = array([NList[0], NList[-1]])
        y = array([Ferrors[0], Ferrors[0] * (NList[0] / NList[-1]) ** (2.)])
        ax1.loglog(x, y, 'b:', linewidth=2, label="2nd order")

        x = array([NList[0], NList[-1]])
        y = array([Ferrors[0], Ferrors[0] * (NList[0] / NList[-1]) ** (3.)])
        ax1.loglog(x, y, 'g-.', linewidth=2, label="3rd order")

        x = array([NList[0], NList[-1]])
        y = array([Ferrors[0], Ferrors[0] * (NList[0] / NList[-1]) ** (4.)])
        ax1.loglog(x, y, 'm--', linewidth=2, label="4th order")

        x = array([NList[0], NList[-1]])
        y = array([Ferrors[0], Ferrors[0] * (NList[0] / NList[-1]) ** (5.)])
        ax1.loglog(x, y, 'r--', linewidth=2, label="5th order")

        # ax1.set_ylim(1e-17, 1e2)
        ax1.set_xlabel('$N$')
        ax1.set_ylabel('$|| F_{numerical} - F_{analytical} ||_{2}$')

        plt.ylim([1e-17, 1e1])

        ax1.legend(loc="best")

        ax1.grid(True)

        fileName = "{}_{}_Global_F_convergence.{}".format(numAlg, motion, self.fileType)
        fileNameWithPath = os.path.join("images", fileName)
        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

        slope = log(Ferrors[0] / Ferrors[3]) / log(NList[0] / NList[-1])
        print('{} {} Deformation Gradient convergence slope = {:.2f}'.format(numAlg, motion, slope))

        # Plots for position errors
        fig = plt.figure()
        matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
        matplotlib.rcParams['font.size'] = 15
        ax2 = fig.gca()
        ax2.loglog(NList, positionErrors, 'k-o', linewidth=2, label="simulation")

        x = array([NList[0], NList[-1]])
        y = array([positionErrors[0], positionErrors[0] * (NList[0] / NList[-1]) ** (1.)])
        ax2.loglog(x, y, 'y--', linewidth=2, label="1st order")

        x = array([NList[0], NList[-1]])
        y = array([positionErrors[0], positionErrors[0] * (NList[0] / NList[-1]) ** (2.)])
        ax2.loglog(x, y, 'b:', linewidth=2, label="2nd order")

        x = array([NList[0], NList[-1]])
        y = array([positionErrors[0], positionErrors[0] * (NList[0] / NList[-1]) ** (3.)])
        ax2.loglog(x, y, 'g-.', linewidth=2, label="3rd order")

        x = array([NList[0], NList[-1]])
        y = array([positionErrors[0], positionErrors[0] * (NList[0] / NList[-1]) ** (4.)])
        ax2.loglog(x, y, 'm--', linewidth=2, label="4th order")

        x = array([NList[0], NList[-1]])
        y = array([positionErrors[0], positionErrors[0] * (NList[0] / NList[-1]) ** (5.)])
        ax2.loglog(x, y, 'r--', linewidth=2, label="5th order")

        # ax2.set_ylim(1e-16, 1e3)
        ax2.set_xlabel('$N$')
        ax2.set_ylabel('$|| x_{numerical} - x_{analytical} ||_{2}$')

        ax2.legend(loc="best")

        # ax2.set_ylim([1e-17, 1e2])

        ax2.grid(True)
        # ax2.axis('tight')
        fileName = "{}_{}_Global_Position_convergence.{}".format(numAlg, motion, self.fileType)
        fileNameWithPath = os.path.join("images", fileName)

        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

        plt.close()

        slope = log(positionErrors[0] / positionErrors[3]) / log(NList[0] / NList[-1])
        print('{} {} Position convergence slope = {:.2f}'.format(numAlg, motion, slope))