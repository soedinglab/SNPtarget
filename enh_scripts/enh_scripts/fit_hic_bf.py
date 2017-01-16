import argparse
import logging
import os
import pickle
from collections import defaultdict
import sys
import glob

import numpy as np
from scipy.integrate import quad
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri as rpyn

import matplotlib as mpl
mpl.use("AGG")
import matplotlib.pyplot as plt

from enh_scripts.utils.bayesfactors import (
    HiCBayesFactor,
    PickleSpline,
    HiCTrueGrenander
)

from enh_scripts import HUMAN_CHROMS
from enh_scripts.utils.utils import read_binning
from enh_scripts import add_logging_args, init_logger
from ngsbiotools.parsers import PredFileParser

logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser('fit_hic_bf')
    parser.add_argument('pred_db')
    parser.add_argument('dist_bf_db')
    parser.add_argument('dist_bin_file')
    parser.add_argument('out_file')
    parser.add_argument('--plot_dir')
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info('invoked script via %s', ' '.join(sys.argv))
    fdrtool = importr('fdrtool')

    logger.info("started reading distance binning")
    with open(args.dist_bin_file) as bin_handle:
        dist_binning = read_binning(bin_handle)
    dist_min, dist_max = int(dist_binning.min), int(dist_binning.max)

    binned_pvals = defaultdict(list)
    eta0_data = defaultdict(list)
    for chrom in HUMAN_CHROMS:
        pred_glob = os.path.join(args.pred_db, '%s.*' % chrom)
        pred_hits = glob.glob(pred_glob)
        if len(pred_hits) != 1:
            logger.warn('No or ambiguous prediction data for %s. Skipping.' % chrom)
            continue
        pred_file, = pred_hits

        dist_bf_glob = os.path.join(args.dist_bf_db, '%s.*' % chrom)
        dist_bf_hits = glob.glob(dist_bf_glob)
        if len(dist_bf_hits) != 1:
            logger.warn('No or ambiguous distance bf for %s. Skipping.' % chrom)
            continue
        dist_bf_file, = dist_bf_hits

        dists = []
        pvals = []
        with open(pred_file) as handle:
            parser = PredFileParser(handle)
            for rec in parser.parse():
                if dist_min <= rec.distance < dist_max:
                    pvals.append(rec.pval_e)
                    dists.append(rec.distance)
                    binned_pvals[rec.dist_bin].append(rec.hic_pval)
        distances = np.abs(np.array(dists))

        # estimate overall eta0 from the correlation p-value distribution
        pval_vec = ro.FloatVector(pvals)
        result = fdrtool.fdrtool(pval_vec, statistic='pvalue', plot=False,
                                 verbose=False)
        [[cutoff, n_cens, chrom_eta0, eta0_se]] = rpyn.ri2py(result.rx2("param"))

        with open(dist_bf_file, 'rb') as bf_pickle:
            dist_bf = pickle.load(bf_pickle)

        for start, end in dist_binning.iter_bins():
            try:
                pos_int, _ = quad(dist_bf.pos_density, start, end)
                neg_int, _ = quad(dist_bf.neg_density, start, end)

                neg_prop = neg_int * chrom_eta0
                pos_prop = pos_int * (1 - chrom_eta0)
                bin_eta0 = neg_prop / (neg_prop + pos_prop)

                N = len(distances[(distances >= start) & (distances <= end)])
                eta0_data[chrom].append((bin_eta0, N))
            except ValueError:
                logger.warn('integrating distbin (%s,%s) on chrom %s failed',
                            start, end, chrom)
                # no eta0 estimation for this distbin
                eta0_data[chrom].append((1, 0))

    logger.info('calculating hic p-value grenander densities')

    true_gren_map = {}
    raw_gren_map = {}
    for distbin, hic_pvals in binned_pvals.items():
        start, end = dist_binning.bin_to_range(distbin)

        pval_ro = ro.FloatVector(hic_pvals)
        hic_fit = fdrtool.fdrtool(pval_ro, statistic='pvalue', plot=False,
                                  verbose=False)
        [[_, _, hic_eta0, _]] = rpyn.ri2py(hic_fit.rx2('param'))
        hic_ecdf = fdrtool.ecdf_pval(pval_ro, eta0=hic_eta0)
        hic_fit = fdrtool.grenander(hic_ecdf, type='decreasing')

        hic_x_knots = rpyn.ri2py(hic_fit.rx2('x.knots'))
        hic_f_knots = rpyn.ri2py(hic_fit.rx2('f.knots'))
        hic_gren = PickleSpline(x=hic_x_knots, y=hic_f_knots, kind='zero',
                                bounds_error=False, fill_value=0)
        raw_gren_map[distbin] = hic_gren

        # calculate dist_bin eta0 as a weighted average over the estimates of all chromosomes
        n_total = 0
        eta_sum = 0
        for eta0_list in eta0_data.values():
            chr_eta0, chr_N = eta0_list[distbin]
            eta_sum += chr_eta0 * chr_N
            n_total += chr_N
        eta0 = eta_sum / n_total

        true_gren = HiCTrueGrenander(hic_gren, eta0)
        true_gren_map[distbin] = true_gren

        if args.plot_dir:
            fname = 'bin%04d_gren_hic.png' % distbin
            outfile = os.path.join(args.plot_dir, fname)
            fig, ax = plt.subplots(figsize=(8, 4))
            x = np.linspace(0, 1, 1e4)
            ax.hist(hic_pvals, bins=20, normed=True, alpha=0.3)
            ax.plot(x, hic_gren(x), color='red')
            ax.plot(x, [eta0] * len(x), color='orange')
            fig.savefig(outfile)
            plt.close(fig)

    logger.info('writing the results')
    hic_bf = HiCBayesFactor(raw_gren_map, true_gren_map, dist_binning)
    with open(args.out_file, 'wb') as out_handle:
        pickle.dump(hic_bf, file=out_handle)
    logger.info('successfully wrote bayes factor to %s', args.out_file)

    logger.info('all done')


if __name__ == '__main__':
    main()
