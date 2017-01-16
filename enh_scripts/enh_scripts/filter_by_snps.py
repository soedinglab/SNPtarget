import argparse
from collections import defaultdict
import logging
import sys

from enh_scripts import add_logging_args, init_logger
from ngsbiotools.ivtree import IVTree, Interval
from ngsbiotools.parsers import PredFileParser
logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser("filter_by_snps")
    parser.add_argument("pred_file")
    parser.add_argument("snp_db")
    parser.add_argument("snp_list")
    parser.add_argument("outfile")
    parser.add_argument("--tol_dist", default=5000)
    parser.add_argument("--posterior_cutoff", type=float, default=0)
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info("invoked script via %r", " ".join(sys.argv))

    logger.info("started reading snp_db")
    snp_db = {}
    with open(args.snp_db) as db_handle:
        for line in db_handle:
            snp_id, chrom, pos_str = line.split()
            snp_db[int(snp_id)] = chrom, int(pos_str)

    logger.info("started preparing the search trees")
    # using a interval tree here probably only pays off for huge snp lists
    rel_snps = defaultdict(IVTree)
    with open(args.snp_list) as list_handle:
        for line in list_handle:
            snp_id = int(line.strip())
            chrom, pos = snp_db[snp_id]
            snp_iv = Interval(pos - args.tol_dist, pos + args.tol_dist)
            snp_iv.snp_pos = pos
            snp_iv.snp_id = snp_id
            rel_snps[chrom].insert(snp_iv)

    logger.info("started processing the prediction file")
    with open(args.outfile, "w") as outfile:
        header = ["snp_id", "snp_pos", "gene_symb", "distance",
                  "dhs_snp_dist", "pval_bf", "dist_bf", "hic_bf", "posterior"]
        print(*header, sep="\t", file=outfile)

        with open(args.pred_file) as pred_handle:
            pparser = PredFileParser(pred_handle)
            for line_no, rec in enumerate(pparser.parse(), start=1):
                chrom = rec.DHS.chrom
                dhs_iv = Interval(rec.DHS.start, rec.DHS.end)
                dhs_mid = rec.DHS.start + (rec.DHS.end - rec.DHS.start) // 2
                snp = rel_snps[chrom].query_overlap(dhs_iv)
                if snp and rec.posterior > args.posterior_cutoff:
                    dhs_snp_dist = abs(dhs_mid - snp.snp_pos)
                    print(snp.snp_id, snp.snp_pos, rec.gene_symb,
                          rec.distance, dhs_snp_dist, rec.pval_e_bf,
                          rec.distance_bf, rec.hic_pval_bf, rec.posterior,
                          sep="\t", file=outfile)
                if line_no % 50000 == 0:
                    logger.info("%s lines processed", line_no)
        logger.info("all done. Bye bye")


if __name__ == '__main__':
    main()
