from abc import ABCMeta, abstractmethod
from scipy.sparse.linalg import expm
from numpy import array, dot, zeros_like
from numpy.linalg import inv
from math import pi


# Just an interface for a motions
class Motion(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def getVel(self, xIJ, time):
        pass

    @abstractmethod
    def getDvDt(self, xIJ, time):
        pass

    @abstractmethod
    def getAnalyticalF(self, time):
        pass

    @abstractmethod
    def getAnalyticalPosition(self, x0, time):
        pass


class Motion1(Motion):

    def __init__(self):
        super().__init__()
        # set global motion parameters
        theta = pi
        self.X0 = array([0.5, 0.5])

        self.Vel0 = array([0.1, 0.0])  # translation velocity

        # calculate rotation matrix
        self.Omega = array([[0.0, -theta],
                            [theta, 0.0]])  # skew symmetric matrix

        # ############
        # self.Vel0 = array([0.1, 0.0])  # translation velocity
        #
        # self.Omega = array([[0.0, -theta],
        #                     [theta, 0.0]])  # skew symmetric matrix
        # #############

    def __str__(self):
        return "motion_1"

    def getVel(self, xIJ, time):
        # print("time:", time)
        return dot(self.Omega, (xIJ - time * self.Vel0)) + self.Vel0

    def getDvDt(self, xIJ, time):
        accn = -self.Omega @ self.Vel0
        return accn

    def getAnalyticalF(self, time):
        return expm(time * self.Omega)  # brute force matrix exponential

    def getAnalyticalPosition(self, x0, time):
        Q = expm(time * self.Omega)
        x = Q @ (x0 - self.X0) + time * self.Vel0
        return x


class Motion2(Motion):

    def __init__(self):
        super().__init__()
        # set global motion parameters
        # gg = 0.5
        # gamma = [gg, 1.-gg]
        # omega =

        self.gamma1 = 0.1
        self.gamma2 = 1. - self.gamma1

        self.Omega1 = array([[0., -1.], [1., 0.]])
        self.Omega2 = array([[0., -1.], [1., 0.]])/2.

        self.X1 = array([0.,0.])
        self.X2 = array([10.,10.])

        self.x0 = self.gamma1 * self.X1 + self.gamma2 * self.X2

        # to be calculated later
        self.Q1 = zeros_like(self.Omega1)
        self.Q2 = zeros_like(self.Omega2)

        self.R = zeros_like(self.Omega1)
        self.Rinv = zeros_like(self.Omega1)

        self.S1 = zeros_like(self.Q1)
        self.S2 = zeros_like(self.Q2)

        self.x1 = zeros_like(self.X1)
        self.x2 = zeros_like(self.X2)

        self.xTilde = zeros_like(self.X1)

    def __str__(self):
        return "motion_2"

    def getVel(self, xIJ, time):
        self.Q1 = expm(time*self.Omega1)
        self.Q2 = expm(time*self.Omega2)
        # print("time:", time)
        # print("O1",self.Omega1)
        # print("O2", self.Omega2)

        self.R = self.gamma1 * self.Q1 + self.gamma2 * self.Q2
        self.Rinv = inv(self.R)
        # print("Q1", self.Q1)
        # print("Q2", self.Q2)

        self.S1 = self.Q1 @ self.Rinv
        self.S2 = self.Q2 @ self.Rinv

        self.x1 = self.Q1 @ self.X1
        self.x2 = self.Q1 @ self.X2

        self.xTilde = self.gamma1 * self.x1 + self.gamma2 * self.x2
        self.xTilde -= self.x0

        v1 = self.Omega1 @ (self.S1 @ (xIJ + self.xTilde) - self.x1)
        v2 = self.Omega2 @ (self.S2 @ (xIJ + self.xTilde) - self.x2)
        v = self.gamma1 * v1 + self.gamma2 * v2

        return v

    def getDvDt(self, xIJ, time):
        print("Krish has work to do")
        raise

    def getAnalyticalF(self, time):
        Q1 = expm(time*self.Omega1)
        Q2 = expm(time*self.Omega2)

        R = self.gamma1 * Q1 + self.gamma2 * Q2
        return R

    def getAnalyticalPosition(self, x0, time):
        Q1 = expm(time*self.Omega1)
        Q2 = expm(time*self.Omega2)

        x = self.gamma1 * Q1 @ (x0 - self.x1) \
            + self.gamma2 * Q2 @ (x0 - self.x2) \
            + self.x0

        return x


