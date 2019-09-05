from sys import platform

import matplotlib
if "win" in platform.lower():
    matplotlib.use('TkAgg')
else:
    matplotlib.use('Agg')
matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
matplotlib.rcParams['font.size'] = 15

import matplotlib.pyplot as plt

import os

from math import ceil
import numpy as np


class ParticleTracePlot(object):
    """
    !class doc

    variables:
        self.fig
        self.size = array([1.0, 1.0])
        self.xNDs = array([])
        self.yNDs = array([])

    methods:
        def __init__(self)
        def __str__(self)
        def __repr__(self)
        def cla(self)
        def setDomain(self, x0, y0, width, height)           ... overwrite default domain size
        def addTraces(self, particleTraceList, timeList=[])  ... add entire traces for a list of particles
        def setGridNodes(self, nodes)
        def exportImage(self, filename)                      ... this is the user interface

    """

    def __init__(self):
        self.size = np.array([1.0, 1.0])
        self.fig = plt.figure()
        self.xNDs = np.array([])
        self.yNDs = np.array([])

    def __str__(self):
        return "Particle trace plotter object"

    def __repr__(self):
        return "ParticleTracePlot()"

    def cla(self):
        ax = self.fig.gca()
        ax.clear()
        self.setGridNodes([])

    def exportImage(self, filename):

        ax = self.fig.gca()

        ax.axis("equal")
        ax.grid(False)

        # create folder to store images
        if not os.path.isdir("images"):
            os.mkdir("images")
        fileNameWithPath = os.path.join("images", filename)

        plt.savefig(fileNameWithPath, pad_inches=0, bbox_inches='tight')

    def setDomain(self, x0, y0, width, height):
        self.X1   = np.array([x0,y0])
        self.size = np.array([width, height])

    def setGridNodes(self, nodes):
        x = []
        y = []
        for rowOfNodes in nodes:
            for node in rowOfNodes:
                x.append(node.pos[0])
                y.append(node.pos[1])
        self.xNDs = np.array(x)
        self.yNDs = np.array(y)

        ax = self.fig.gca()
        ax.scatter(self.xNDs, self.yNDs, marker='+', s=30, c="b")

    def addTraces(self, particleTraceList, timeList=[]):
        """
        add entire traces for a list of particles

        particleTraceList = [{'node':id, 'path':array([[x,y], ...])}, ...]
        """

        ax = self.fig.gca()

        for trace in particleTraceList:

            name = trace['node']
            x = trace['path'][:,0]
            y = trace['path'][:,1]

            ax.plot(x, y, '--', c="b")
            ax.scatter(x[0], y[0], s=30, c="r")

