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
    def getAnalyticalF(self, xIJ, time):
        pass

    @abstractmethod
    def getAnalyticalPosition(self, x0, time):
        pass


class Motion1(Motion):

    def __init__(self):
        super().__init__()
        # set global motion parameters
        theta = pi/2.0
        self.X0 = array([0.5, 0.5])

        self.Vel0 = array([0.1, 0.1])   # translation velocity

        # calculate rotation matrix
        self.Omega = array([[0.0, -theta],
                            [theta, 0.0]])  # skew symmetric matrix

    def __str__(self):
        return "motion_1"

    def getVel(self, xIJ, time):
        v = self.Omega @ (xIJ - self.X0 - time * self.Vel0) + self.Vel0
        # print(v)
        return v

    def getDvDt(self, xIJ, time):
        return -self.Omega @ self.Vel0

    def getAnalyticalF(self, xIJ, time):
        Q = expm(time * self.Omega)
        # print(Q)
        return Q  # brute force matrix exponential

    def getAnalyticalPosition(self, x0, time):
        Q = expm(time * self.Omega)
        x = Q @ (x0 - self.X0) + time * self.Vel0 + self.X0
        return x


class Motion2(Motion):

    def __init__(self):
        super().__init__()

        # self.gamma1 = 0.1
        # self.gamma2 = 1. - self.gamma1

        # self.Omega1 = array([[0., -1.], [1., 0.]])*pi/4.
        # self.Omega2 = array([[0., -1.], [1., 0.]])*pi/6.
        #
        # self.X1 = array([0.5, 0.5])
        # self.X2 = array([5.0, 5.0])

        self.gamma1 = 0.5
        self.gamma2 = 1. - self.gamma1

        self.Omega1 = array([[0., -1.], [1., 0.]])*pi/4.
        self.Omega2 = array([[0., -1.], [1., 0.]])*pi/40.

        self.X1 = array([0.5, 0.5])
        self.X2 = array([0.5, 0.5])

        # self.gamma1 = 0.1
        # self.gamma2 = 1. - self.gamma1
        #
        # self.Omega1 = array([[0., -1.], [1., 0.]])*pi
        # self.Omega2 = array([[0., -1.], [1., 0.]])*pi/3
        #
        # self.X1 = array([0.5, 0.5])
        # self.X2 = array([0.5, 0.5])

        self.x0 = self.gamma1 * self.X1 + self.gamma2 * self.X2

    def __str__(self):
        return "motion_2"

    def getVel(self, xIJ, time):
        Q1 = expm(time*self.Omega1)
        Q2 = expm(time*self.Omega2)

        R = self.gamma1 * Q1 + self.gamma2 * Q2
        Rinv = inv(R)

        S1 = Q1 @ Rinv
        S2 = Q2 @ Rinv

        x1 = Q1 @ self.X1
        x2 = Q2 @ self.X2

        xTilde = self.gamma1 * x1 + self.gamma2 * x2
        xTilde -= self.x0

        v1 = self.Omega1 @ (S1 @ (xIJ + xTilde) - x1)
        v2 = self.Omega2 @ (S2 @ (xIJ + xTilde) - x2)
        v = self.gamma1 * v1 + self.gamma2 * v2

        # print(v)
        return v

    def getDvDt(self, xIJ, time):
        Q1 = expm(time*self.Omega1)
        Q2 = expm(time*self.Omega2)

        R = self.gamma1 * Q1 + self.gamma2 * Q2
        Rinv = inv(R)

        S1 = Q1 @ Rinv
        S2 = Q2 @ Rinv

        x1 = Q1 @ self.X1
        x2 = Q2 @ self.X2

        xTilde = self.gamma1 * x1 + self.gamma2 * x2
        xTilde -= self.x0

        v1 = self.Omega1 @ (S1 @ (xIJ + xTilde) - x1)
        v2 = self.Omega2 @ (S2 @ (xIJ + xTilde) - x2)

        gradv1 = self.Omega1 @ S1
        gradv2 = self.Omega2 @ S2
        gradv  = self.gamma1 * gradv1 + self.gamma2 * gradv2

        dvdt = self.gamma1 * (self.Omega1 @ v1) + self.gamma2 * (self.Omega2 @ v2)
        dvdt -= dot(gradv, self.getVel(xIJ, time))

        # print(dvdt)
        return dvdt


    def getAnalyticalF(self, xIJ, time):
        Q1 = expm(time*self.Omega1)
        Q2 = expm(time*self.Omega2)

        R = self.gamma1 * Q1 + self.gamma2 * Q2
        # print(R)
        return R

    def getAnalyticalPosition(self, x0, time):
        Q1 = expm(time*self.Omega1)
        Q2 = expm(time*self.Omega2)

        x = self.gamma1 * Q1 @ (x0 - self.X1) \
            +self.gamma2 * Q2 @ (x0 - self.X2) \
            +self.x0

        return x


