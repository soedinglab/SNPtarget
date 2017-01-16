import argparse
import logging
import os
import sys

import numpy as np
from ngsbiotools.parsers import PredFileParser, DHS_Parser
from enh_scripts import add_logging_args, init_logger
logger = logging.getLogger()


def load_dhs_db(db_path, chrom):
    dhs_access = {}
    dhs_fpath = os.path.join(db_path, chrom + ".dhs.gz")
    with DHS_Parser(dhs_fpath) as dhs_parser:
        for dhs in dhs_parser.read():
            dhs_access[dhs.id] = np.sum(dhs.access_vals)
    return dhs_access


def main():
    parser = argparse.ArgumentParser('extract_access')
    parser.add_argument('dhs_db')
    parser.add_argument('pred_file')
    parser.add_argument('out_file')
    add_logging_args(parser)
    args = parser.parse_args()
    init_logger(args)

    logger.info("starting...")
    logger.info("invoked script via %r", " ".join(sys.argv))

    dhs_access_dict = {}

    logger.info("started processing input file")
    with open(args.pred_file) as pred_handle, open(args.out_file, 'w') as out_handle:
        header = pred_handle.readline().split()
        header.extend(['total_access'])
        print(*header, sep='\t', file=out_handle)
        pred_handle.seek(0)

        parser = PredFileParser(pred_handle)
        for line_no, rec in enumerate(parser.parse(), start=1):
            chrom = rec.DHS.chrom
            dhs_id = rec.DHS_id

            if chrom not in dhs_access_dict:
                dhs_access_dict[chrom] = load_dhs_db(args.dhs_db, chrom)
                logger.info("loaded dhs database for chromosome '%s'", chrom)

            if dhs_id not in dhs_access_dict[chrom]:
                logger.warn("Unknown dhs_id '%s'", dhs_id)
                continue
            total_access = dhs_access_dict[chrom][dhs_id]

            tokens = parser.lastline.split('\t')
            tokens.append(total_access)
            print(*tokens, sep='\t', file=out_handle)

            if line_no % 50000 == 0:
                logger.info("%s lines processed", line_no)

    logger.info("all done. Exiting.")


def extract_background(matrix, i, j, m, n, bw):
    i_ind = np.arange(i-bw, i+bw+1)
    j_ind = np.arange(j-bw, j+bw+1)
    ind = (i_ind >= 0) & (j_ind >= 0) & (i_ind < m) & (j_ind < n)
    hic_bg = matrix[(i_ind[ind], j_ind[ind])].toarray()
    return hic_bg

if __name__ == '__main__':
    main()
