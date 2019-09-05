# ====== settings ================

PLOT_MOTIONS           = True
PLOT_SINGLE_STEP_TESTS = False
PLOT_MULTI_STEP_TESTS  = False

COLLATE_PLOTS = False

MOTION1 = True
MOTION2 = True
MOTION3 = True
MOTION4 = True

ALGORITHM_EXPLICIT    = True
ALGORITHM_MIDPOINT    = True
ALGORITHM_RUNGE_KUTTA = True

OUTPUT_FILE_TYPE = 'png'

NUM_CELLS = 8

# ====== the test function =======

from LocalConvergenceTest import *
from GlobalConvergenceTest import *

from Motion import *
from MotionPlot import *
from ErrorPlotter import *


def testSingleStep(theMotion):

    if COLLATE_PLOTS:
        localErrorPlot = ErrorPlotter(NUM_CELLS)
        if ALGORITHM_EXPLICIT:
            localErrorPlot.addTestData(LocalConvergenceTest(theMotion, ExplicitEuler(), OUTPUT_FILE_TYPE, NUM_CELLS))

        if ALGORITHM_MIDPOINT:
            localErrorPlot.addTestData(LocalConvergenceTest(theMotion, MidPointRule(), OUTPUT_FILE_TYPE, NUM_CELLS))

        if ALGORITHM_RUNGE_KUTTA:
            localErrorPlot.addTestData(LocalConvergenceTest(theMotion, RungeKutta4(), OUTPUT_FILE_TYPE, NUM_CELLS))
        localErrorPlot.savePlot(theMotion)

    else:
        if ALGORITHM_EXPLICIT:
            LocalConvergenceTest(theMotion, ExplicitEuler(), OUTPUT_FILE_TYPE, nCells=NUM_CELLS).runAnalysis()

        if ALGORITHM_MIDPOINT:
            LocalConvergenceTest(theMotion, MidPointRule(), OUTPUT_FILE_TYPE, nCells=NUM_CELLS).runAnalysis()

        if ALGORITHM_RUNGE_KUTTA:
            LocalConvergenceTest(theMotion, RungeKutta4(), OUTPUT_FILE_TYPE, nCells=NUM_CELLS).runAnalysis()



def testMultipleSteps(theMotion):

    if COLLATE_PLOTS:
        globalErrorPlot = ErrorPlotter(NUM_CELLS)
        if ALGORITHM_EXPLICIT:
            globalErrorPlot.addTestData(GlobalConvergenceTest(theMotion, ExplicitEuler(), OUTPUT_FILE_TYPE))

        if ALGORITHM_MIDPOINT:
            globalErrorPlot.addTestData(GlobalConvergenceTest(theMotion, MidPointRule(), OUTPUT_FILE_TYPE))

        if ALGORITHM_RUNGE_KUTTA:
            globalErrorPlot.addTestData(GlobalConvergenceTest(theMotion, RungeKutta4(), OUTPUT_FILE_TYPE))

        globalErrorPlot.savePlot(theMotion)

    else:
        if ALGORITHM_EXPLICIT:
            GlobalConvergenceTest(theMotion, ExplicitEuler(), OUTPUT_FILE_TYPE, nCells=NUM_CELLS).runAnalysis()

        if ALGORITHM_MIDPOINT:
            GlobalConvergenceTest(theMotion, MidPointRule(), OUTPUT_FILE_TYPE, nCells=NUM_CELLS).runAnalysis()

        if ALGORITHM_RUNGE_KUTTA:
            GlobalConvergenceTest(theMotion, RungeKutta4(), OUTPUT_FILE_TYPE, nCells=NUM_CELLS).runAnalysis()


def plotMotionTraces():

    if MOTION1:
        m = MotionPlot(Motion1())
        #m.setGrid(NUM_CELLS,NUM_CELLS)
        m.setMaxTime(4.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m1.png")

    if MOTION2:
        m = MotionPlot(Motion2())
        #m.setGrid(NUM_CELLS,NUM_CELLS)
        m.setMaxTime(10.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m2.png")

    if MOTION3:
        m = MotionPlot(Motion3())
        #m.setGrid(NUM_CELLS,NUM_CELLS)
        m.setMaxTime(8.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m3.png")

    if MOTION4:
        m = MotionPlot(Motion4())
        #m.setGrid(NUM_CELLS,NUM_CELLS)
        m.setMaxTime(20.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.5, .1], [0.5, .3], [0.5, .5], [0.5, .7], [0.5, .9]))
        m.exportImage("m4.png")

        m = MotionPlot(Motion4())
        #m.setGrid(NUM_CELLS,NUM_CELLS)
        m.setMaxTime(20.)
        m.setPointsPerSecond(10)
        m.setTracers(([0.1, .5], [0.3, .5], [0.5, .5], [0.7, .5], [0.9, .5]))
        m.exportImage("m4b.png")


def Main():

    if PLOT_SINGLE_STEP_TESTS:

        if MOTION1:
            testSingleStep(Motion1())

        if MOTION2:
            testSingleStep(Motion2())

        if MOTION3:
            testSingleStep(Motion3())

        if MOTION4:
            testSingleStep(Motion4())


    if PLOT_MULTI_STEP_TESTS:

        if MOTION1:
            testMultipleSteps(Motion1())

        if MOTION2:
            testMultipleSteps(Motion2())

        if MOTION3:
            testMultipleSteps(Motion3())

        if MOTION4:
            testMultipleSteps(Motion4())


    if PLOT_MOTIONS:
        plotMotionTraces()


if __name__ == '__main__':
    Main()
