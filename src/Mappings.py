import numpy as np

class Mappings(object):
    """
    @class: Mapping

    variables:

    methods:
        def __init__(self)
        def identity(self, idx, xi, params=dict())
        def fineEdge(self, idx, xi, params=dict())
    """

    def __init__(self):
        pass

    def identity(self, idx, xi, params=dict()):
        return xi

    def fineEdge(self, idx, xi, params=dict()):

        d = 1.

        if idx == 0:
            if 'width' in params.keys():
                d = params['width']

        if idx == 1:
            if 'height' in params.keys():
                d = params['height']

        X = d * (1.0 - np.cos(np.pi * xi)) / 2.

        return X