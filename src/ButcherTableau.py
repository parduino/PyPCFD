from abc import ABCMeta, abstractmethod
from numpy import array


# Just an interface for a Butcher Table
class ButcherTableau(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_a(self):
        pass

    @abstractmethod
    def get_b(self):
        pass

    @abstractmethod
    def get_c(self):
        pass


class ExplicitEuler(ButcherTableau):

    def __init__(self):
        super().__init__()

    def __str__(self):
        return "ExplicitEuler"

    def get_a(self):
        return array([0.])       # time factors

    def get_b(self):
        return array([[0.]])  # position factors

    def get_c(self):
        return array([1.])       # update factors


class MidPointRule(ButcherTableau):

    def __init__(self):
        super().__init__()

    def __str__(self):
        return "MidPointRule"

    def get_a(self):
        return array([0., 1./2.])       # time factors

    def get_b(self):
        return array([[0., 0.],
                      [1./2., 0.]])  # position factors

    def get_c(self):
        return array([0., 1.])       # update factors


class RungeKutta4(ButcherTableau):

    def __init__(self):
        super().__init__()

    def __str__(self):
        return "RK4"

    def get_a(self):
        return array([0., 1./2., 1./2., 1.])         # time factors

    def get_b(self):
        return array([[0. ,0. ,0.,0.],
                      [0.5,0. ,0.,0.],
                      [0. ,0.5,0.,0.],
                      [0. ,0. ,1.,0.]])              # position factors

    def get_c(self):
        return array([1./6., 1./3., 1./3., 1./6.])   # update factors


class MidPointMethod2(ButcherTableau):

    def __init__(self):
        super().__init__()

    def __str__(self):
        return "AltMidPointRule"

    def get_a(self):
        return array([0., 1.])        # time factors

    def get_b(self):
        return array([[0., 0.],
                       [1., 0.]])     # position factors

    def get_c(self):
        return array([1./2., 1./2.])  # update factors


# class BogackiShampine(ButcherTableau):
#
#     def __init__(self):
#         super().__init__()
#
#     def get_a(self):
#         return array([0., 1./2., 3./4., 1.])        # time factors
#
#     def get_b(self):
#         return array([[0., 0., 0., 0.],
#                       [1./2., 0., 0, 0.],
#                       [0., 3./4., 0., 0.],
#                       [2./9., 1./3., 4./9., 0.]])     # position factors
#
#     def get_c(self):
#         return array([[2./9.,	1./3.,	4./9.,	0.],
#                      [7./.24,	1./.4,	1./.3,	1./.8]])  # update factors
