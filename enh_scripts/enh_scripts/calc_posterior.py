import argparse
import pickle
import logging
from functools import lru_cache

from ngsbiotools.parsers import PredFileParser
from enh_scripts import init_logger, add_logging_args

logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser(
        prog="calc_posterior",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("in_file", help="input file")
    parser.add_argument("out_file", help="output file")
    parser.add_argument("pval_bf", help="seralized p-value bayes factor")
    parser.add_argument("dist_bf", help="serialized distance bayes factor")
    parser.add_argument("hic_bf", help="serialized hic bayes factor")

    parser.add_argument("--pval_col_name", default="pval_e",
                        help="name of the p-value column")
    parser.add_argument("--dist_col_name", default="distance",
                        help="name of the distance column")
    parser.add_argument("--hicpval_col_name", default="hic_pval",
                        help="name of the hic p-value column")
    add_logging_args(parser)
    return parser


def main():
    arg_parser = create_parser()
    args = arg_parser.parse_args()
    init_logger(args)

    pval_col = args.pval_col_name
    dist_col = args.dist_col_name
    hic_col = args.hicpval_col_name
    post_col = "posterior"
    pval_bf_col = "%s_bf" % pval_col
    dist_bf_col = "%s_bf" % dist_col
    hic_bf_col = "%s_bf" % hic_col

    with open(args.pval_bf, "rb") as pval_pickle:
        pval_bf = pickle.load(pval_pickle)
    with open(args.dist_bf, "rb") as dist_pickle:
        dist_bf = pickle.load(dist_pickle)
    with open(args.hic_bf, "rb") as hic_pickle:
        hic_bf = pickle.load(hic_pickle)
    eta0 = pval_bf.bg_cut

    with open(args.in_file, "r") as in_handle,\
            open(args.out_file, "w") as out_handle:

        header = in_handle.readline().split()
        header.extend([pval_bf_col, dist_bf_col, hic_bf_col, post_col])
        print(*header, sep='\t', file=out_handle)
        in_handle.seek(0)

        @lru_cache(maxsize=None)
        def eval_hicbf(dist, hic_pval):
            return hic_bf(dist, hic_pval)

        @lru_cache(maxsize=None)
        def eval_pvalbf(pval):
            return pval_bf(pval)

        @lru_cache(maxsize=None)
        def eval_distbf(dist):
            return dist_bf(dist)

        parser = PredFileParser(in_handle)
        for line_no, rec in enumerate(parser.parse(), start=1):
            pval = rec.pval_e
            dist = abs(rec.distance)
            hic_pval = rec.hic_pval

            pval_bf_inv = eval_pvalbf(pval)
            dist_bf_inv = eval_distbf(dist)
            try:
                hic_bf_inv = eval_hicbf(dist, hic_pval)
            # distance out of range
            except ValueError:
                continue

            bf_inv = pval_bf_inv * dist_bf_inv * hic_bf_inv * (1 - eta0) / eta0
            posterior = bf_inv / (bf_inv + 1)

            tokens = parser.lastline.split('\t')
            tokens.extend([pval_bf_inv, dist_bf_inv, hic_bf_inv, posterior])
            print(*tokens, sep='\t', file=out_handle)
            if line_no % 50000 == 0:
                logger.info("%s lines processed", line_no)
    logger.info("all done. Exiting")


if __name__ == "__main__":
    main()
