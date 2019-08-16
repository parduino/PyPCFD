from LocalConvergenceTest import *
from GlobalConvergenceTest import *
from numpy import array, linspace, zeros, zeros_like, tensordot
import matplotlib.pyplot as plt

def Main():
    fileType = 'pdf'
    # fileType = 'png'
    LocalConvergenceTest(fileType).runAnalysis()
    # GlobalConvergenceTest(fileType).runAnalysis()
    # Motion1().plotMotion()
    # Motion2().plotMotion()
    # m = Motion3()

    # m. getLagrangianPosition(array([1., 0.]), 0.0)
    # m.plotMotion()

if __name__ == '__main__':
    Main()
