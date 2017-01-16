import argparse
import logging
import pickle
import sys

import numpy as np
from ngsbiotools.parsers import PredFileParser
from enh_scripts import add_logging_args, init_logger
logger = logging.getLogger()


def main():
    parser = argparse.ArgumentParser('extract_hic')
    parser.add_argument('hic_matrix')
    parser.add_argument('pred_file')
    parser.add_argument('out_file')
    parser.add_argument('--bg_bandwidth', type=int)
    parser.add_argument('--full_bg', action="store_true")
    add_logging_args(parser)
    args = parser.parse_args()
    init_logger(args)

    logger.info("starting...")
    logger.info("invoked script via %r", " ".join(sys.argv))

    # load hic matrix from binary dump
    with open(args.hic_matrix, 'rb') as matrix_dump:
        matrix_obj = pickle.load(matrix_dump)
        assert matrix_obj.matrix_type == "band"
        hic_matrix = matrix_obj.matrix
        m, n = matrix_obj.shape
        resolution = matrix_obj.resolution
    logger.info("successfully loaded the HiC matrix")

    if args.bg_bandwidth:
        bw = args.bg_bandwidth // resolution
    full_bg = args.full_bg
    if full_bg:
        bg_cache = {}

    logger.info("started processing input file")
    with open(args.pred_file) as pred_handle, open(args.out_file, 'w') as out_handle:
        header = pred_handle.readline().split()
        header.extend(['hic_i', 'hic_j', 'hic_strength'])
        print(*header, sep='\t', file=out_handle)
        pred_handle.seek(0)

        parser = PredFileParser(pred_handle)
        hic_cache = {}
        for line_no, rec in enumerate(parser.parse(), start=1):
            tss_bin = rec.gene.tss // resolution
            start_bin = rec.DHS.start // resolution

            i, j = start_bin, tss_bin
            if i > j:
                i, j = j, i

            if full_bg:
                bg_hic = bg_cache.get(j - i, None)
                if bg_hic is None:
                    bg_hic = np.nanmean(hic_matrix[j - i].todense())
                    bg_cache[j - i] = bg_hic
                try:
                    hic = hic_matrix[j - i][0, i] - bg_hic
                except IndexError:
                    logger.warn('interaction in line: %s not covered by HiC matrix. Skipping',
                                line_no)
                    continue
            elif args.bg_bandwidth:
                if (i, j) in hic_cache:
                    hic = hic_cache[i, j]
                else:
                    hic_bg = extract_background(hic_matrix, i, j, m, n, bw)
                    hic = hic_matrix[j - i][0, i] - np.nanmean(hic_bg)
                    hic_cache[i, j] = hic
            else:
                try:
                    hic = hic_matrix[j - i][0, i]
                except IndexError:
                    logger.warn('interaction in line: %s not covered by HiC matrix. Skipping',
                                line_no)
                    continue

            # skip, if we don't have hic information
            if np.isfinite(hic):
                tokens = parser.lastline.split('\t')
                tokens.extend([i, j, hic])
                print(*tokens, sep='\t', file=out_handle)
                if line_no % 50000 == 0:
                    logger.info("%s lines processed", line_no)
    logger.info("all done. Exiting.")


def extract_background(diag_dict, i, j, m, n, bw):
    d = j - i
    high = min(m, n - d)
    hic_vec = diag_dict[d]
    hic_bg = hic_vec[0, max(i - bw, 0):min(i + bw + 1, high)].toarray().ravel()
    return hic_bg

if __name__ == '__main__':
    main()
