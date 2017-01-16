import argparse
from collections import defaultdict
import logging
import sys


import numpy as np
from ngsbiotools.parsers import PredFileParser
from enh_scripts.utils.utils import read_binning
from enh_scripts import add_logging_args, init_logger

logger = logging.getLogger()


def main():
    parser = argparse.ArgumentParser('hic_pval')
    parser.add_argument('pred_file')
    parser.add_argument('hic_resolution', type=int)
    parser.add_argument('distbin_file')
    parser.add_argument('out_file')
    parser.add_argument('--pval_thresh', type=float, default=0.3)

    add_logging_args(parser)
    args = parser.parse_args()
    init_logger(args)
    resolution = args.hic_resolution

    logger.info("starting...")
    logger.info("invoked script via '%s'", " ".join(sys.argv))

    logger.info("started reading distance binning")
    with open(args.distbin_file) as bin_handle:
        dist_binning = read_binning(bin_handle)

    logger.info("started reading input file")
    with open(args.pred_file) as pred_handle:
        bin_pval_dict = defaultdict(list)
        bin_hic_dict = {}
        parser = PredFileParser(pred_handle)
        for rec in parser.parse():
            tss_bin = rec.gene.tss // resolution
            start_bin = rec.DHS.start // resolution

            i, j = (start_bin, tss_bin) if start_bin <= tss_bin else (tss_bin, start_bin)
            bin_pval_dict[i, j].append(rec.pval_e)
            bin_hic_dict[i, j] = rec.hic_strength
    logger.info("finished phase 1: aggregating data for each bin")

    bg_ind = defaultdict(lambda: 0)
    bg_hic_dict_tmp = defaultdict(list)
    thresh = args.pval_thresh
    for (i, j), pvals in bin_pval_dict.items():
        try:
            distbin = dist_binning.dist_to_bin((j - i) * resolution)
        except ValueError:
            continue

        # count multiple times
        # if np.min(pvals) >= thresh:
        #    bg_ind[i, j] += 1
        #    bg_hic_dict_tmp[distbin].append(bin_hic_dict[i, j])
        n_insignificant = np.sum(np.array(pvals) >= thresh)
        bg_ind[i, j] += n_insignificant
        bg_hic_dict_tmp[distbin].extend([bin_hic_dict[i, j]] * n_insignificant)

    logger.info("finished phase 2: determining background bins")
    del bin_pval_dict

    bg_hic_dict = {}
    for distbin, hic_list in bg_hic_dict_tmp.items():
        hic_array = np.array(hic_list)
        hic_array.sort()
        bg_hic_dict[distbin] = hic_array
    del bg_hic_dict_tmp
    logger.info("finised phase 3: sorting hic values in background bins")

    logger.info("started processing input file")
    with open(args.pred_file) as pred_handle, open(args.out_file, 'w') as out_handle:
        header = pred_handle.readline().split()
        header.extend(['hic_pval', 'bg_count', 'dist_bin'])
        print(*header, sep='\t', file=out_handle)
        pred_handle.seek(0)

        parser = PredFileParser(pred_handle)
        for line_no, rec in enumerate(parser.parse(), start=1):
            tss_bin = rec.gene.tss // resolution
            start_bin = rec.DHS.start // resolution

            i, j = (start_bin, tss_bin) if start_bin <= tss_bin else (tss_bin, start_bin)

            hic_val = bin_hic_dict[(i, j)]
            try:
                distbin = dist_binning.dist_to_bin((j - i) * resolution)
            except ValueError:
                continue
            hic_bg = bg_hic_dict[distbin]
            # adding jitter could be a cleaner solution here.
            ins_pos_left = np.searchsorted(hic_bg, hic_val, side='left')
            ins_pos_right = np.searchsorted(hic_bg, hic_val, side='right')
            ins_pos = np.random.randint(ins_pos_left, ins_pos_right + 1)
            hic_pval = 1 - ins_pos / len(hic_bg)

            bg_count = bg_ind[i, j]

            tokens = parser.lastline.split('\t')
            tokens.extend([hic_pval, bg_count, distbin])
            print(*tokens, sep='\t', file=out_handle)
            if line_no % 50000 == 0:
                logger.info("%s lines processed", line_no)
    logger.info("all done. Exiting.")


if __name__ == '__main__':
    main()
