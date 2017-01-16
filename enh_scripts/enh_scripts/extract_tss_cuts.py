import argparse
import pickle
import logging
import sys
import os
from multiprocessing import Pool
from collections import defaultdict

import numpy as np
import pysam

from ngsbiotools.parsers import DNaseDataParser
from enh_scripts import init_logger, add_logging_args

logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser("extract_cuts")
    parser.add_argument("trscr_file")
    parser.add_argument("data_info_file")
    parser.add_argument("output_file")
    parser.add_argument("--bam_db", default=".")
    parser.add_argument("--upstream_bp", type=int, default=150)
    parser.add_argument("--n_proc", type=int, default=1)
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)
    upstream_bp = args.upstream_bp

    logger.info("invoked script via %r", " ".join(sys.argv))

    tss_md = []
    tss_locs = []
    with open(args.trscr_file) as trscr_file:
        # remove header
        # trscr_file.readline()
        for line in trscr_file:
            toks = line.split()
            gene_name = toks[6]
            chrom = toks[0]
            tss = int(toks[4])
            strand = toks[3]
            if strand == '+':
                start = tss - upstream_bp
                end = tss
            else:
                start = tss
                end = tss + upstream_bp

            md = {}
            md["id"] = gene_name
            md["chrom"] = chrom
            md["loc"] = tss
            md["flanking"] = upstream_bp
            tss_md.append(md)
            tss_locs.append((chrom, (start, end)))

    logger.info("successfully read tss database to memory")

    tissue_names = []
    bam_files = []
    lab_list = []
    with open(args.data_info_file) as dnase_handle:
        dnase_parser = DNaseDataParser(dnase_handle)
        for rec in dnase_parser.parse():
            tissue_names.append(rec.tissue)
            bam_files.append(rec.bam_files)
            lab_list.append(rec.batch)

    logger.info("creating cut matrix")
    jobs = []
    access_dict = defaultdict(list)
    with Pool(args.n_proc) as pool:
        for tissue, bam_list in zip(tissue_names, bam_files):
            for rep_bam in bam_list:
                bam_path = os.path.join(args.bam_db, rep_bam)
                params = [bam_path, tss_locs]
                result = pool.apply_async(extract_cuts, params)
                jobs.append((tissue, rep_bam, result))

        for tissue, rep_bam, result in jobs:
            access_dict[tissue].append(result.get())
            logger.info("finished processing %s", rep_bam)

    with open(args.output_file, "wb") as pkl_handle:
        for tissue, data in access_dict.items():
            for vector in data:
                assert len(vector) == len(tss_md)
        output = (access_dict, tissue_names, tss_md, lab_list)
        pickle.dump(output, pkl_handle)
    logger.info("successfully exported data matrix")
    logger.info("all done - good bye")


def extract_cuts(bam_path, extr_regions):
    d = len(extr_regions)
    cuts = np.zeros(d)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for tss_no, (chrom, (start, end)) in enumerate(extr_regions):
            for read in bam.fetch(chrom, start, end):
                # cuts are always at the 5' end
                if read.is_reverse:
                    cut = read.reference_end
                else:
                    cut = read.reference_start
                if start <= cut < end:
                    cuts[tss_no] += 1
    return cuts


if __name__ == '__main__':
    main()
