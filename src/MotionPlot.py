from sys import platform

import matplotlib
if "win" in platform.lower():
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')

#matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
#matplotlib.rcParams['font.size'] = 15

import matplotlib.pyplot as plt

import os

from math import ceil
import numpy as np

from Node import *


class MotionPlot(object):
    """
    !class doc

    variables:
        self.motion
        self.particles = [ array([0.5, 0.5]) ]
        self.X1 = array([0.0,0.0])
        self.size = array([1.0, 1.0])
        self.npts    ... points per trace
        self.maxTime ... max time for trace
        self.xNDs = np.array(x)
        self.yNDs = np.array(y)

    methods:
        def __init__(self, motion)
        def __str__(self)
        def __repr__(self)
        def setMaxTime(self, T)
        def setPointsPerSecond(self, n)
        def setGridNodes(self, nodes):
        def exportImage(self, filename)            ... this is the user interface
        def setDomain(self, x0, y0, width, height) ... overwrite default domain size
        def setTracers(self, particleList)         ... overwrite default tracer particles
        def createTraces(self, ax, tracers)        ... for internal use only!
    """

    def __init__(self, motion):
        self.motion = motion
        self.particles = [ np.array([0.5, 0.5]) ]
        self.X1 = np.array([0.0,0.0])
        self.size = np.array([1.0, 1.0])

        self.npts    = 30
        self.maxtime = 1.0

        self.xNDs = np.array([])
        self.yNDs = np.array([])

    def __str__(self):
        if isinstance(self.motion, Motion):
            s = "plotter for " + str(self.motion)
        else:
            s = "plotter for no motion - deactivated"
        return s

    def __repr__(self):
        return "MotionPlot({})".format(repr(self.motion))

    def setMaxTime(self, T):
        self.maxtime = T

    def setPointsPerSecond(self, n):
        if n>=1:
            self.npts = n

    def exportImage(self, filename):

        # plt.ion()   # WHY?  The initial setup selects non-interactive to begin with!
        fig = plt.figure()
        ax = fig.gca()

        ax.axis("equal")
        ax.grid(False)

        if (len(self.xNDs)>0):
            ax.scatter(self.xNDs, self.yNDs, marker='+', s=30, c="b")

        self.createTraces(ax, self.particles)

        # create folder to store images
        if not os.path.isdir("images"):
            os.mkdir("images")
        fileNameWithPath = os.path.join("images", filename)

        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')
        plt.close()

    def setDomain(self, x0, y0, width, height):
        self.X1   = np.array([x0,y0])
        self.size = np.array([width, height])

    def setTracers(self, particleList):
        self.particles = particleList

    def createTraces(self, ax, tracers):

        t = np.linspace(0.0, self.maxtime, num=ceil(self.npts*self.maxtime))

        for pos in tracers:

            x0 = np.array(pos)

            xlocs = []
            ylocs = []

            for time in t:
                xp = self.motion.getAnalyticalPosition(x0, time)
                xlocs.append(xp[0])
                ylocs.append(xp[1])

            # ax.scatter(xlocs, ylocs, s=2, c="b")
            ax.plot(xlocs, ylocs, '--', c="b")
            ax.scatter(x0[0], x0[1], s=30, c="r")

    def setGrid(self, nCellsX, nCellsY):
        theNodes = []
        id = 0

        x = np.linspace(0, self.size[0], (nCellsX + 1))
        y = np.linspace(0, self.size[1], (nCellsY + 1))

        for i in range(nCellsX+1):
            for j in range(nCellsY+1):
                id += 1
                theNode = Node(id,x[i],y[j])
                theNode.setGridCoordinates(i,j)
                theNodes.append(theNode)

        self.setGridNodes(theNodes)

    def setGridNodes(self, nodes):
        x = []
        y = []
        for node in nodes:
            x.append(node.pos[0])
            y.append(node.pos[1])
        self.xNDs = np.array(x)
        self.yNDs = np.array(y)