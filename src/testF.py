from LocalConvergenceTest import *
from GlobalConvergenceTest import *
from numpy import array, linspace, zeros, zeros_like
import matplotlib.pyplot as plt

def Main():
    fileType = 'pdf'
    # fileType = 'png'
    LocalConvergenceTest(fileType).runAnalysis()
    # GlobalConvergenceTest(fileType).runAnalysis()
    # Motion1().plotMotion()
    # Motion2().plotMotion()



#     ax = addData()
#     ax.scatter(2,3)
#     matplotlib.rcParams['font.sans-serif'] = "Times New Roman"
#     matplotlib.rcParams['font.size'] = 15
#     ax.axis("equal")
#     ax.grid(True)
#     plt.show()
#
# def addData():
#     fig, ax = plt.subplots()
#     ax.scatter(1, 2)
#     return ax



if __name__ == '__main__':
    Main()
