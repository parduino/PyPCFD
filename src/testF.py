# ====== settings ================

PLOT_MOTIONS = True

PLOT_SINGLE_STEP_TESTS = False

PLOT_MULTI_STEP_TESTS  = True

MOTION1 = False
MOTION2 = False
MOTION3 = True

ALGORITHM_EXPLICIT    = False
ALGORITHM_MIDPOINT    = True
ALGORITHM_RUNGE_KUTTA = True

OUTPUT_FILE_TYPE = 'png'

# ====== the test function =======

from LocalConvergenceTest import *
from GlobalConvergenceTest import *

from Motion import *
from MotionPlot import *

def Main():
    fileType = OUTPUT_FILE_TYPE

    if PLOT_SINGLE_STEP_TESTS:

        if MOTION1:
            if ALGORITHM_EXPLICIT:
                LocalConvergenceTest(Motion1(), ExplicitEuler(), fileType).runAnalysis()

            if ALGORITHM_MIDPOINT:
                LocalConvergenceTest(Motion1(), MidPointRule(), fileType).runAnalysis()

            if ALGORITHM_RUNGE_KUTTA:
                LocalConvergenceTest(Motion1(), RungeKutta4(), fileType).runAnalysis()

        if MOTION2:
            if ALGORITHM_EXPLICIT:
                LocalConvergenceTest(Motion2(), ExplicitEuler(), fileType).runAnalysis()

            if ALGORITHM_MIDPOINT:
                LocalConvergenceTest(Motion2(), MidPointRule(), fileType).runAnalysis()

            if ALGORITHM_RUNGE_KUTTA:
                LocalConvergenceTest(Motion2(), RungeKutta4(), fileType).runAnalysis()

        if MOTION3:
            if ALGORITHM_EXPLICIT:
                LocalConvergenceTest(Motion3(), ExplicitEuler(), fileType).runAnalysis()

            if ALGORITHM_MIDPOINT:
                LocalConvergenceTest(Motion3(), MidPointRule(), fileType).runAnalysis()

            if ALGORITHM_RUNGE_KUTTA:
                LocalConvergenceTest(Motion3(), RungeKutta4(), fileType).runAnalysis()


    if PLOT_MULTI_STEP_TESTS:

        if MOTION1:
            if ALGORITHM_EXPLICIT:
                GlobalConvergenceTest(Motion1(), ExplicitEuler(), fileType).runAnalysis()

            if ALGORITHM_MIDPOINT:
                GlobalConvergenceTest(Motion1(), MidPointRule(), fileType).runAnalysis()

            if ALGORITHM_RUNGE_KUTTA:
                GlobalConvergenceTest(Motion1(), RungeKutta4(), fileType).runAnalysis()

        if MOTION2:
            if ALGORITHM_EXPLICIT:
                GlobalConvergenceTest(Motion2(), ExplicitEuler(), fileType).runAnalysis()

            if ALGORITHM_MIDPOINT:
                GlobalConvergenceTest(Motion2(), MidPointRule(), fileType).runAnalysis()

            if ALGORITHM_RUNGE_KUTTA:
                GlobalConvergenceTest(Motion2(), RungeKutta4(), fileType).runAnalysis()

        if MOTION3:
            if ALGORITHM_EXPLICIT:
                GlobalConvergenceTest(Motion3(), ExplicitEuler(), fileType).runAnalysis()

            if ALGORITHM_MIDPOINT:
                GlobalConvergenceTest(Motion3(), MidPointRule(), fileType).runAnalysis()

            if ALGORITHM_RUNGE_KUTTA:
                GlobalConvergenceTest(Motion3(), RungeKutta4(), fileType).runAnalysis()


    if PLOT_MOTIONS:

        if MOTION1:
            m = MotionPlot(Motion1())
            m.setMaxTime(4.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m1.png")

        if MOTION2:
            m = MotionPlot(Motion2())
            m.setMaxTime(10.)
            m.setPointsPerSecond(10)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m2.png")

        if MOTION3:
            m = MotionPlot(Motion3())
            m.setMaxTime(40.)
            m.setPointsPerSecond(4)
            m.setTracers( ([0.5, .1],[0.5, .3],[0.5, .5],[0.5, .7],[0.5, .9]) )
            m.exportImage("m3.png")


if __name__ == '__main__':
    Main()
