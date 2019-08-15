from LocalConvergenceTest import *
from GlobalConvergenceTest import *
from numpy import array, linspace, zeros, zeros_like

def Main():
    fileType = 'pdf'
    # fileType = 'png'
    # LocalConvergenceTest(fileType).runAnalysis()
    # GlobalConvergenceTest(fileType).runAnalysis()
    Motion2().plotMotion()



if __name__ == '__main__':
    Main()
