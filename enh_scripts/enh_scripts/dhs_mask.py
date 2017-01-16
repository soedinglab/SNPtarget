import argparse
import gzip
from itertools import chain

import numpy as np
from scipy.stats import zmap
from pandas.stats.moments import rolling_median

from ngsbiotools.parsers import DHS_Parser

parser = argparse.ArgumentParser("")
parser.add_argument("dhs_file")
parser.add_argument("out_file")
parser.add_argument("--zscore_cutoff", type=float, default=2)
parser.add_argument("--zscore_rank_exclude", type=int, default=5)
parser.add_argument("--window_size", type=int, default=351)
args = parser.parse_args()

rank_excl = args.zscore_rank_exclude
zscore_vecs = []
thresh = args.zscore_cutoff

with DHS_Parser(args.dhs_file) as dhs_parser:
    for dhs in dhs_parser.read():
        # calculate z-score for each accessibilty value
        order = np.argsort(dhs.access_vals)
        acc = dhs.access_vals[order]
        zscore_vec = zmap(dhs.access_vals, acc[rank_excl:-rank_excl])
        zscore_vecs.append(zscore_vec)
zscore_mat = np.vstack(zscore_vecs).T

# apply running median
result = rolling_median(zscore_mat, window=args.window_size, axis=1)

# masking
n, d = zscore_mat.shape
win_size = args.window_size
mask = np.zeros((n,d), dtype=bool)
halfwin = win_size // 2

for i in range(n):
    for j in range(halfwin, d - halfwin):
        if result[i,j] > thresh:
            mask[i,j-halfwin:j+halfwin] = 1

acc_field_start = 4
with gzip.open(args.dhs_file, "rt") as in_file,\
        gzip.open(args.out_file, "wt") as out_file:
    header = in_file.readline().split()
    print(*header, sep='\t', file=out_file)
    for i, line in enumerate(in_file):
        toks = line.split()
        access_vals = np.array(toks[acc_field_start:])
        access_vals[mask[:,i]] = np.nan
        print(*chain(toks[:acc_field_start], access_vals), file=out_file)
