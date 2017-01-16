import argparse
import pickle
import sys
import logging
import glob

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import fmin

from rpy2.robjects.packages import importr, STAP
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rpyn

from enh_scripts.utils.bayesfactors import (
    DistanceBayesFactor,
    PickleSpline,
)
from enh_scripts.utils.utils import ECDF_DIST_RFUN
from enh_scripts import add_logging_args, init_logger

logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser("fit_pval_bf")
    parser.add_argument("pred_file")
    parser.add_argument("--pval_enrich_thresh", type=float, default=0.1)
    parser.add_argument("--pval_bg_thresh", type=float, default=0.6)
    parser.add_argument("--min_dist", type=int, default=0)
    parser.add_argument("--max_dist", type=int, default=2000000)
    parser.add_argument("--plot_file")
    parser.add_argument('--lsq_pval_estim', action='store_true')
    parser.add_argument('--qval', type=float)
    parser.add_argument("out_file")
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info('invoked script via %r', ' '.join(sys.argv))
    logger.info('started reading input file')

    dfs = []
    for pred_file in glob.glob(args.pred_file):
        df = pd.read_table(pred_file, usecols=['distance', 'pval_e'])
        logger.info('successfully read %r', pred_file)
        dfs.append(df)
    df = pd.concat(dfs)

    min_dist = args.min_dist
    max_dist = args.max_dist
    df = df[(np.abs(df.distance) >= min_dist) & (np.abs(df.distance) <= max_dist)]

    def prepare_dist(dist):
        dist = np.abs(dist)
        pseudo_dist = []
        if not np.any(dist == min_dist):
            pseudo_dist.append(min_dist)
        if not np.any(dist == max_dist):
            pseudo_dist.append(max_dist)
        return pd.concat((pd.Series(pseudo_dist), dist))

    e_thresh = args.pval_enrich_thresh
    bg_thresh = args.pval_bg_thresh
    enrich_dist = prepare_dist(df.loc[df.pval_e < e_thresh, 'distance'])
    neg_dist = prepare_dist(df.loc[df.pval_e > bg_thresh, 'distance'])

    neg_dist = prepare_dist(df.loc[(df.pval_e > 0.5) & (df.pval_e < 0.8), 'distance'])
    y, bins = np.histogram(enrich_dist, bins=20, normed=True)
    min_f = y.min() * (max_dist - min_dist)

    x_uniq = pd.concat((enrich_dist, neg_dist)).sort_values().unique()
    enrich_dist_vec = ro.FloatVector(enrich_dist)
    neg_dist_vec = ro.FloatVector(neg_dist)
    x_vec = ro.FloatVector(x_uniq)

    ecdf_fun = ro.r['ecdf']
    enrich_ecdf = ecdf_fun(enrich_dist_vec)
    neg_ecdf = ecdf_fun(neg_dist_vec)

    # estimate the q-value
    logger.info('calculating q-value from p-values')
    fdrtool = importr('fdrtool')
    rpkg = STAP(ECDF_DIST_RFUN, 'enh_scripts')

    pval_vec = ro.FloatVector(df.pval_e)
    result = fdrtool.fdrtool(pval_vec, statistic='pvalue', plot=False,
                             verbose=False)

    fdr_pvals = rpyn.ri2py(result.rx2('pval'))
    qvals = rpyn.ri2py(result.rx2('qval'))

    sort_order = np.argsort(fdr_pvals)
    fdr_pvals_sort = fdr_pvals[sort_order]
    qvals_sort = qvals[sort_order]

    thresh_index = np.searchsorted(fdr_pvals_sort, e_thresh)
    qval = qvals_sort[thresh_index]

    logger.info('estimated qval is %.4f', qval)

    # use least-square fit to estimate the qval
    ind_end = len(enrich_ecdf(x_vec)) - 1
    ind_start = int(0.6 * ind_end)

    x_ind = np.linspace(ind_start, ind_end, 50).astype(int)
    y_all = rpyn.ri2py(enrich_ecdf(x_vec))
    y_neg = rpyn.ri2py(neg_ecdf(x_vec))

    dx = np.diff(x_uniq[x_ind])
    dy = np.diff(y_all[x_ind])
    dz = np.diff(y_neg[x_ind])

    # choose qval such that the gradient of qval * neg matches the gradient
    # of the curve consisting of enriched predictions
    def opt_fun(qval):
        return np.sum(((dy - qval * dz) / dx)**2)

    qval_lsq, = fmin(opt_fun, x0=qval, disp=0)
    if args.lsq_pval_estim:
        qval = qval_lsq
    if args.qval is not None:
        qval = np.float64(args.qval)

    logger.info('least-square estimated qval is %.4f', qval_lsq)
    # if we fail this assertion, there are serious problems with the data
    # and we are not going to continue

    # moving to l_bfgs_b allows to use bounds, but this should not
    # be necessary.
#    assert 0 <= qval_lsq <= 1

    raw_tp = enrich_ecdf(x_vec) - qval * neg_ecdf(x_vec)

    # the fraction of positives at all times can be at maximum 1-qval
    y_tp = np.minimum(raw_tp, 1 - qval)

    # y_tp is an estimate of the cummulative function of TP, therefore it is
    # monotonically increasing and normed to 1 on the right side.
    y_tp = np.maximum.accumulate(y_tp)
    true_ecdf = y_tp / (1 - qval)

    true_ecdf_vec = ro.FloatVector(true_ecdf)
    dist_ecdf = ro.r('''
        function (x, y)
        {
            rval <- approxfun(x, y, method = "constant", yleft = 0,
                              yright = 1, f = 0, ties = "ordered")
            class(rval) = c("ecdf", "stepfun", class(rval))
            attr(rval, "call") <- sys.call()
            rval
        }''')

    tp_ecdf = dist_ecdf(x_vec, true_ecdf_vec)
    tp_fit = fdrtool.grenander(tp_ecdf, type="decreasing")

    tp_x = rpyn.ri2py(tp_fit.rx2('x.knots'))
    tp_y = rpyn.ri2py(tp_fit.rx2('f.knots'))

    enrich_ecdf = rpkg.ecdf_dists(enrich_dist_vec, min_f, min_dist, max_dist)
    enrich_gren = fdrtool.grenander(enrich_ecdf, type="decreasing")
    enrich_x = rpyn.ri2py(enrich_gren.rx2('x.knots'))
    enrich_y = rpyn.ri2py(enrich_gren.rx2('f.knots'))

    f = interp1d(enrich_x, enrich_y, kind='zero')
    f_tp = interp1d(tp_x, tp_y, kind='zero')

    f_fp = lambda x: (f(x) - (1 - qval) * f_tp(x)) / qval
    neg_y = f_fp(x_uniq)

    pos_dens = PickleSpline(tp_x, tp_y, kind='zero')
    neg_dens = PickleSpline(x_uniq, neg_y, kind='zero')

    dist_bf = DistanceBayesFactor(pos_dens, neg_dens)

    logger.info('writing bayes factor to output file')
    with open(args.out_file, 'wb') as out_handle:
        pickle.dump(dist_bf, file=out_handle)

    if args.plot_file:
        logger.info('plotting results')
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        x_points = np.linspace(min_dist, max_dist, 1000)
        fig, axes = plt.subplots(3, 2, figsize=(14, 12))
        ((cum_pos_ax, cum_ax), (neg_ax, pos_ax), (all_ax, bf_ax)) = axes

        # cumulative positives
        cum_pos_ax.set_title('Cumul. Positives')
        cum_pos_ax.plot(x_uniq, true_ecdf, label='corrected')
        cum_pos_ax.plot(x_uniq, raw_tp / (1 - qval), label='raw')
        cum_pos_ax.legend(loc=4)

        # cummulative curves
        cum_ax.set_title('Cumul. Curves')
        x = rpyn.ri2py(x_vec)
        cum_ax.plot(x, rpyn.ri2py(enrich_ecdf(x_vec)), label='all')
        cum_ax.plot(x, qval * rpyn.ri2py(neg_ecdf(x_vec)), label='neg')
        cum_ax.legend(loc=2)
        cum_ax.set_xlim(1.8e6)
        cum_ax.set_ylim(0.8)

        # density of the positives
        pos_ax.set_title('Positives')
        pos_ax.plot(x_points, f_tp(x_points), color='red')

        # density of the negatives
        neg_ax.set_title('Negatives')
        neg_ax.hist(neg_dist, normed=True, alpha=0.5, bins=20)
        neg_ax.plot(x_points, f_fp(x_points), color='red')

        # density of the combined predictions
        all_ax.set_title('All')
        all_ax.hist(enrich_dist, normed=True, alpha=0.5, bins=20)
        all_ax.plot(x_points, f(x_points), color='red')

        # Bayes factor
        bf_ax.set_title('Bayes factor')
        bf_ax.plot(x_points, f_tp(x_points) / f_fp(x_points), color='red')

        plt.savefig(args.plot_file)

    logger.info('all done. exiting...')

if __name__ == '__main__':
    main()
