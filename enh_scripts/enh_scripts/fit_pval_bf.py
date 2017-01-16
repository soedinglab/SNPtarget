import argparse
import pickle

import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rpyn

from enh_scripts.utils.bayesfactors import (
    PvalueBayesFactor,
    PickleSpline,
    PvalTrueGrenander
)
from ngsbiotools.parsers import PredFileParser


def create_parser():
    parser = argparse.ArgumentParser("fit_pval_bf")
    parser.add_argument("pred_file")
    parser.add_argument("--eta0", type=float)
    parser.add_argument("out_file")
    parser.add_argument("--plot_file")
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    pvals = []
    with open(args.pred_file) as handle:
        parser = PredFileParser(handle)
        for rec in parser.parse():
            pvals.append(rec.pval_e)

    fdrtool = importr("fdrtool")
    pval_vec = ro.FloatVector(pvals)

    eta0 = args.eta0
    if eta0 is None:
        stat = ro.StrVector(["pvalue"])
        result = fdrtool.fdrtool(pval_vec, statistic=stat, plot=False,
                                 verbose=False)
        [[cutoff, n_cens, eta0, eta0_se]] = rpyn.ri2py(result.rx2("param"))

    constr_ecdf = fdrtool.ecdf_pval(pval_vec, eta0=eta0)
    constr_fit = fdrtool.grenander(constr_ecdf, type="decreasing")

    x_knots = rpyn.ri2py(constr_fit.rx2("x.knots"))
    f_knots = rpyn.ri2py(constr_fit.rx2("f.knots"))

    raw_gren = PickleSpline(x=x_knots, y=f_knots, kind="zero",
                            bounds_error=False, fill_value=0)
    true_gren = PvalTrueGrenander(raw_gren, eta0)

    pval_bf = PvalueBayesFactor(raw_gren, true_gren, eta0)

    with open(args.out_file, "wb") as out_handle:
        pickle.dump(pval_bf, file=out_handle)

    if args.plot_file:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 8))
        ax1.hist(pvals, normed=True, alpha=0.5, bins=20)
        xx = np.linspace(0, 1, 1000)
        ax1.plot(xx, raw_gren(xx), color="red")
        ax2.plot(xx, true_gren(xx), color="red")

        plt.savefig(args.plot_file)

if __name__ == '__main__':
    main()
