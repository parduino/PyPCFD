import numpy as np

class Mappings(object):
    """
    @class: Mapping

    variables:
        self.params

    methods:
        def __init__(self, parameters={})
        def toX(self, s, t)
        def toY(self, s, t)
        def toS(self, x, y)
        def toT(self, x, y)
        def oneToX(self, s, t)
        def oneToY(self, s, t)
        def oneToS(self, x, y)
        def oneToT(self, x, y)
    """

    def __init__(self, parameters={}):
        self.params = parameters

    # universal functions
    '''
    The following functions are to recognize different types of
    input such that the mapping can be applied to single points,
    list of points, or even multi-dimensional arrays.
    '''
    def toX(self, s, t):
        if isinstance(s, np.float64) and isinstance(t, np.float64):
            return self.oneToX(s,t)

        res = []
        for (S,T) in zip(s,t):
            res.append(self.toX(S,T))

        return np.array(res)

    def toY(self, s, t):
        if isinstance(s, np.float64) and isinstance(t, np.float64):
            return self.oneToY(s,t)

        res = []
        for (S,T) in zip(s,t):
            res.append(self.toY(S,T))

        return np.array(res)

    def toS(self, x, y):
        if isinstance(x, np.float64) and isinstance(y, np.float64):
            return self.oneToS(x,y)

        res = []
        for (X,Y) in zip(x,y):
            res.append(self.toS(X,Y))

        return np.array(res)

    def toT(self, x, y):
        if isinstance(x, np.float64) and isinstance(y, np.float64):
            return self.oneToT(x,y)

        res = []
        for (X,Y) in zip(x,y):
            res.append(self.toT(X,Y))

        return np.array(res)

    # place holder functions: identities
    '''
    These functions privide the mapping of a single point
    onto a resulting mapped single point
    '''
    def oneToX(self, s, t):
        return s

    def oneToY(self, s, t):
        return t

    def oneToS(self, x, y):
        return x

    def oneToT(self, x, y):
        return y


class IdentityMap(Mappings):

    def oneToX(self, s, t):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return d * s

    def oneToY(self, s, t):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return d * t

    def oneToS(self, x, y):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return x / d

    def oneToT(self, x, y):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return y / d


class FineEdgeMap(Mappings):

    def oneToX(self, s, t):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return d * (1.0 - np.cos(np.pi * s)) / 2.

    def oneToY(self, s, t):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return d * (1.0 - np.cos(np.pi * t)) / 2.

    def oneToS(self, x, y):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return np.acos( 1.0 - 2.*x/d ) / np.pi

    def oneToT(self, x, y):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return np.acos( 1.0 - 2.*y/d ) / np.pi
