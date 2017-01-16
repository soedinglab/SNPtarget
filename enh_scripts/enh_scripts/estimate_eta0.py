import argparse

import pandas as pd
import numpy as np

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rpyn


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("--min_dist", type=float, default=0)
    parser.add_argument("--max_dist", type=float, default=2e6)
    args = parser.parse_args()

    fdrtool = importr("fdrtool")
    min_dist = args.min_dist
    max_dist = args.max_dist

    print("filepath", "cor_eta0", "dist_eta0", sep="\t")

    for f in args.files:
        df = pd.read_table(f)
        pval_cor = df.pval_e
        dists = np.abs(df.distance)

        trafo_dists = (dists - min_dist) / (max_dist - min_dist)

        cor_vec = ro.FloatVector(pval_cor)
        dist_vec = ro.FloatVector(trafo_dists)

        stat = ro.StrVector(["pvalue"])
        cor_res = fdrtool.fdrtool(cor_vec, statistic=stat, plot=False,
                                  verbose=False)
        [[_, _, cor_eta0, _]] = rpyn.ri2py(cor_res.rx2("param"))
        dist_res = fdrtool.fdrtool(dist_vec, statistic=stat, plot=False,
                                   verbose=False)

        [[_, _, dist_eta0, _]] = rpyn.ri2py(dist_res.rx2("param"))

        print(f, cor_eta0, dist_eta0, sep="\t")


if __name__ == '__main__':
    main()
