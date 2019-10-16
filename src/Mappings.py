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
    """

    def __init__(self, parameters={}):
        self.params = parameters

    def toX(self, s, t):
        return s

    def toY(self, s, t):
        return t

    def toS(self, x, y):
        return x

    def toT(self, x, y):
        return y


class IdentityMap(Mappings):

    def toX(self, s, t):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return d * s

    def toY(self, s, t):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return d * t

    def toS(self, x, y):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return x / d

    def toT(self, x, y):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return y / d


class FineEdgeMap(Mappings):

    def toX(self, s, t):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return d * (1.0 - np.cos(np.pi * s)) / 2.

    def toY(self, s, t):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return d * (1.0 - np.cos(np.pi * t)) / 2.

    def toS(self, x, y):

        if 'width' in self.params.keys():
            d = self.params['width']
        else:
            d = 1.0

        return np.acos( 1.0 - 2.*x/d ) / np.pi

    def toT(self, x, y):

        if 'height' in self.params.keys():
            d = self.params['height']
        else:
            d = 1.0

        return np.acos( 1.0 - 2.*y/d ) / np.pi
