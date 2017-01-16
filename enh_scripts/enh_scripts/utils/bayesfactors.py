from scipy.interpolate import interp1d
import numpy as np


class PvalueBayesFactor:

    def __init__(self, gren_raw, pos_distr, bg_cut):
        self._gren_raw = gren_raw
        self._pos_distr = pos_distr
        self._bg_cut = bg_cut

    @property
    def grenander_raw(self):
        return self._gren_raw

    @property
    def positive_distribution(self):
        return self._pos_distr

    @property
    def bg_cut(self):
        return self._bg_cut

    def __call__(self, pval):
        return self._pos_distr(pval)


class DistanceBayesFactor:

    def __init__(self, pos_density, neg_density):
        self._pos_density = pos_density
        self._neg_density = neg_density

    @property
    def pos_density(self):
        return self._pos_density

    @property
    def neg_density(self):
        return self._neg_density

    def __call__(self, dist):
        return self._pos_density(dist) / self._neg_density(dist)


class HiCBayesFactor:

    def __init__(self, raw_gren_fits, true_gren_fits, dist_binning):
        self._gren_raw = raw_gren_fits
        self._pos_distrs = true_gren_fits
        self._binning = dist_binning

    def grenander_raw(self, distbin):
        return self._gren_raw[distbin]

    def positive_distribution(self, distbin):
        return self._pos_distrs[distbin]

    def __call__(self, distance, hic_pval):
        distbin = self._binning.dist_to_bin(distance)
        return self._pos_distrs[distbin](hic_pval)


class PickleSpline():

    def __init__(self, x, y, **args):
        self._x = x
        self._y = y
        self._args = args
        self._spline = interp1d(x, y, **args)

    def __call__(self, val):
        return self._spline(val)

    def __getstate__(self):
        return self._x, self._y, self._args

    def __setstate__(self, state):
        x, y, args = state
        self._x = x
        self._y = y
        self._args = args
        self._spline = interp1d(x, y, **args)


class PvalTrueGrenander:

    def __init__(self, raw_grenander, cut_point):
        self._raw = raw_grenander
        self._cut = cut_point

    def __call__(self, pval):
        return np.maximum(self._raw(pval) - self._cut, 0) / (1 - self._cut)


class DistanceTrueGrenander:

    def __init__(self, raw_grenander, cut_point, support_range):
        self._raw = raw_grenander
        self._cut = cut_point
        self._support = support_range

    def __call__(self, dist):
        return np.maximum(self._raw(dist) - self._cut, 0) / (1 - self._cut * self._support)


class HiCTrueGrenander:

    def __init__(self, raw_grenander, cut_point):
        self._raw = raw_grenander
        self._cut = cut_point

    def __call__(self, hic_val):
        if 1 - self._cut < 1e-6:
            if np.isscalar(hic_val):
                return 0
            else:
                return np.zeros(len(hic_val))
        return np.maximum(self._raw(hic_val) - self._cut, 0) / (1 - self._cut)
