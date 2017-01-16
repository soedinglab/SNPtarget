import sys

import numpy as np


class StreamReader:

    def __init__(self, path):
        self._path = path
        self._stdin = path == "-"

    def __enter__(self):
        if self._stdin:
            self._handle = sys.stdin
        else:
            self._handle = open(self._path)
        return self._handle

    def __exit__(self, ex_type, ex_val, tb):
        if not self._stdin:
            self._handle.close()


def open_or_stdin(path):
    return StreamReader(path)


class DataSink:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def write(self, x):
        pass


def open_or_sink(file, mode="w"):
    if file is not None:
        return open(file, mode)
    else:
        return DataSink()


ECDF_DIST_RFUN = """
ecdf.dists <- function (x, fdr, xmin, xmax)
{
    # compute empirical CDF as usual
    x = sort(x)
    n = length(x)
    if (n < 1)
        stop("'x' must have 1 or more non-missing values")
    vals = sort(unique(x))
    F.raw = cumsum(tabulate(match(x, vals)))/n

    # control upper bound of F:
    # make sure that the maximum slope of (Grenander) F is eta0
    F.raw = pmin(F.raw, 1 - fdr / (xmax - xmin) * (xmax - vals))

    # control lower bound of F:
    # make sure that (Grenander F) >= eta0*vals
    F.raw = pmax(F.raw, fdr / (xmax - xmin) * (vals - xmin))

    # if necessary add an atom at 1 to make it a proper CDF
    if (vals[length(vals)] != xmax)
    {
       F.raw = c(F.raw, 1)
       vals = c(vals, xmax)
    }

    # if necessary also add an atom at 0 with weight zero to get support [0,1]
    if (vals[1] != xmin)
    {
       F.raw = c(0, F.raw)
       vals = c(xmin, vals)
    }

    # finally, modify F such that the last slope of the Grenander F
    # is *exactly* eta0
    i = length(vals)-1
    F.raw[i] = 1 - fdr / (xmax - xmin) * (xmax - vals[i])

    rval <- approxfun(vals, F.raw,
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) = c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    rval
}
"""


class DistanceBinning:

    def __init__(self, bins):
        self._bins = np.array(bins)

    def dist_to_bin(self, distance):
        ins_pos = np.searchsorted(self._bins, distance, side="right")
        if ins_pos == 0 or ins_pos == len(self._bins):
            raise ValueError("%s is out of range for [%s, %s]",
                             distance, self._bins[0], self._bins[-1])
        return ins_pos - 1

    def bin_to_range(self, distbin):
        start, end = self._bins[distbin:distbin + 2]
        return start, end - 1

    def iter_bins(self):
        for i in range(len(self._bins) - 1):
            yield self.bin_to_range(i)

    @property
    def min(self):
        return self._bins[0]

    @property
    def max(self):
        return self._bins[-1]


def read_binning(handle):
    bins = []
    for line in handle:
        [tok] = line.split()
        bins.append(int(tok))
    return DistanceBinning(bins)


def write_binning(binning, handle):
    for tick in binning._bins:
        print(tick, file=handle)
