from scipy.sparse.linalg import expm
from numpy import array, dot, zeros, zeros_like, linspace, tensordot, pi
from numpy.linalg import inv, norm

#from scipy.optimize import newton

import os


# Just an interface for a motions
class Motion(object):

    def __init__(self):
        self.id = -1

    def __str__(self):
        return "motion_{}".format(self.id)

    def __repr__(self):
        return "Motion{}()".format(self.id)

    def getVel(self, xIJ, time):
        return array([0.0, 0.0])

    def getDvDt(self, xIJ, time):
        return array([0.0, 0.0])

    def getAnalyticalF(self, x0, time):
        return zeros((2,2))

    def getAnalyticalPosition(self, x0, time):
        return array([0.0, 0.0])


class Motion1(Motion):

    def __init__(self):
        super().__init__()
        # set global motion parameters
        theta = pi
        self.X0 = array([0.5, 0.5])

        self.Vel0 = array([0.1, 0.1])   # translation velocity

        # calculate rotation matrix
        self.Omega = array([[0.0, -theta],
                            [theta, 0.0]])  # skew symmetric matrix

    def getVel(self, xIJ, time):
        v = self.Omega @ (xIJ - self.X0 - time * self.Vel0) + self.Vel0
        # print(v)
        return v

    def getDvDt(self, xIJ, time):
        return -self.Omega @ self.Vel0

    def getAnalyticalF(self, x0, time):
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

        self.gamma1 = 0.6
        self.gamma2 = 1. - self.gamma1

        self.Omega1 = array([[0., -1.], [1., 0.]])*pi/4.
        self.Omega2 = array([[0., -1.], [1., 0.]])*pi/40.

        self.X1 = array([0.5, 0.5])
        self.X2 = array([0.8, 0.8])

        self.x0 = self.gamma1 * self.X1 + self.gamma2 * self.X2

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


    def getAnalyticalF(self, x0, time):
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


class Motion3(Motion):

    def __init__(self):
        super().__init__()
        # set global motion parameters
        theta = pi/20.0
        self.X0 = array([0.5, 0.5])

        # initialize rotation matrix
        self.Omega = array([[0.0, -theta],
                            [theta, 0.0]])  # skew symmetric matrix

        self.tolerance = 1.e-14

    def getVel(self, xIJ, time):
        X = self.getLagrangianPosition(xIJ, time)
        R = self.getR(X, time)
        v = (X @ X) * self.Omega @ R @ X
        return v

    def getDvDt(self, xIJ, time):
        X = self.getLagrangianPosition(xIJ, time)
        R = self.getR(X, time)
        V = self.getVel(X, time)
        # F = self.getAnalyticalF(X, time)  ... avoid multiple computation of expensive tensor products !

        r2 = dot(X,X)
        dVdt   = r2 * self.Omega @ V

        Y = self.Omega @ R @ X
        Z = tensordot(Y, X, axes=0)
        F = R + 2.0 * time * Z
        GradV  = r2 * self.Omega @ F + 2.0 * Z

        # grad_v = dVdt - GradV @ V    # WTF ???  del v is a 2nd-order tensor, not a vector !!!!!
        delV = GradV @ inv(F)
        dvdt   = dVdt - delV @ V

        return dvdt

    def getR(self, X, t):
        r2 = dot(X,X)
        return expm(r2*t * self.Omega)

    def getAnalyticalF(self, X, time):
        R = self.getR(X, time)
        Y = self.Omega @ R @ X
        F = R + 2.0 * time * tensordot(Y, X, axes=0)
        return F

    def getAnalyticalPosition(self, X, time):
        R = self.getR(X, time)
        return R @ X

    def zeroFunc(self, X, xIJ, time):
        return xIJ - self.getAnalyticalPosition(X, time)

    def funcPrime(self, X, xIJ, time):
        return -self.getAnalyticalF(X, time)

    def getLagrangianPosition(self, xIJ, time):
        # WHY? Xk = xIJ/2.
        R = self.getR(xIJ, time)
        X = xIJ @ R
        error = xIJ - self.getAnalyticalPosition(X, time)

        cnt = 0
        while norm(error) > self.tolerance:
            self.F = self.getAnalyticalF(X, time)
            deltaX = inv(self.F) @ error
            X += deltaX

            cnt += 1

            # error = F @ deltaX
            error = xIJ - self.getAnalyticalPosition(X, time)
            print("* step {}: error = {:12.6e}".format(cnt,norm(error)))

            if cnt > 10:
                print("Newton iteration failed to converge")
                raise

        msg = "Eulerian Position = {}, Lagrangian Position = {}, Reconstructed Eulerian Position = {}"
        print(msg.format(xIJ, X, self.getAnalyticalPosition(X, time)))

        # X = newton(self.zeroFunc, xIJ, fprime=self.funcPrime, args=(xIJ,time))
        return X


