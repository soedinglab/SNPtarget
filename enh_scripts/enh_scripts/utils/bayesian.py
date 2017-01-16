import numpy as np


class PvalStats:

    def __init__(self, pvals, dens, eta0):
        self._pvals = pvals
        self._dens = dens
        self._eta0 = eta0

    def density(self, pval):
        return np.interp(pval, self._pvals, self._dens)

    @property
    def eta0(self):
        return self._eta0


class DistStats:

    def __init__(self, dists, dens, min_dist, max_dist, qval):
        self._dists = dists
        self._dens = dens
        self._min_dist = min_dist
        self._max_dist = max_dist
        self._qval = qval

    def density(self, dist):
        return np.interp(dist, self._dists, self._dens)

    @property
    def min_dist(self):
        return self._min_dist

    @property
    def max_dist(self):
        return self._max_dist

    @property
    def qval(self):
        return self._qval
