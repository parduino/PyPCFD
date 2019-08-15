from LocalConvergenceTest import *
from GlobalConvergenceTest import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from numpy import array, linspace, zeros, zeros_like

def Main():
    fileType = 'pdf'
    # fileType = 'png'
    # LocalConvergenceTest(fileType).runAnalysis()
    # GlobalConvergenceTest(fileType).runAnalysis()
    m = Motion2()
    plotMotions(m)


def plotMotions(m):
    n = 1000
    t = linspace(0, 200, num=n)
    x0 = array([0.5, 1./3.])
    lims = 0.5
    xlocs = zeros_like(t)
    xlocs[0] = x0[0]
    ylocs = zeros_like(t)
    ylocs[0] = x0[1]

    # You probably won't need this if you're embedding things in a tkinter plot...
    plt.ion()
    fig = plt.figure()
    matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
    matplotlib.rcParams['font.size'] = 15
    ax3 = fig.gca()
    ax3.axis("equal")
    ax3.grid(True)


    for j in range(1, n):
        xp = m.getAnalyticalPosition(x0, t[j])
        xlocs[j] = xp[0]
        ylocs[j] = xp[1]
        ax3.scatter(xlocs[0:j], ylocs[0:j], s=10, c="b")
        ax3.scatter(xlocs[0], ylocs[0], s=30, c="r")
        ax3.scatter(m.X1[0], m.X1[1], s=30, c="k")
        ax3.scatter(m.X2[0], m.X2[1], s=30, c="m")
        # ax3.set_xlim(-lims, lims)
        # ax3.set_ylim(-lims, lims)

        plt.draw()
        plt.pause(0.0001)
        ax3.clear()

    # ax3.scatter(xlocs, ylocs, s=10)

    plt.show()


if __name__ == '__main__':
    Main()
