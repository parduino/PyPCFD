from LocalConvergenceTest import *
from GlobalConvergenceTest import *

def Main():
    # fileType = 'pdf'
    fileType = 'png'
    LocalConvergenceTest(fileType).runAnalysis()
    GlobalConvergenceTest(fileType).runAnalysis()


if __name__ == '__main__':
    Main()
