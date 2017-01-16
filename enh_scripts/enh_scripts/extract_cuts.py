import argparse
import pickle
import logging
import sys
import os
from multiprocessing import Pool
from collections import defaultdict

import numpy as np
import pysam

from ngsbiotools.parsers import SNPdbParser, DNaseDataParser
from enh_scripts import init_logger, add_logging_args

logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser("extract_cuts")
    parser.add_argument("snp_db")
    parser.add_argument("snp_list")
    parser.add_argument("data_info_file")
    parser.add_argument("output_file")
    parser.add_argument("--bam_db", default=".")
    parser.add_argument("--snp_flanking", type=int, default=75)
    parser.add_argument("--n_proc", type=int, default=1)
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info("invoked script via %r", " ".join(sys.argv))

    logger.info("started reading snp database")
    snp_db = {}
    with open(args.snp_db) as snpdb_handle:
        snpdb_parser = SNPdbParser(snpdb_handle)
        for rec in snpdb_parser.parse():
            snp_db[rec.snp_id] = (rec.chrom, rec.location)

    snp_md = []
    snp_locs = []
    with open(args.snp_list) as snplist_handle:
        for line in snplist_handle:
            snp_id = line.strip()
            snp_loc = snp_db.get(snp_id, None)
            if snp_loc is None:
                logger.warn("could not find snp with id %s in the database", snp_id)
            else:
                chrom, loc = snp_loc
                md = {}
                md["id"] = snp_id
                md["chrom"] = chrom
                md["loc"] = loc
                md["flanking"] = args.snp_flanking
                snp_md.append(md)
                snp_locs.append(snp_loc)
    logger.info("successfully read snp database to memory")

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
    flanking = args.snp_flanking
    jobs = []
    access_dict = defaultdict(list)
    with Pool(args.n_proc) as pool:
        for tissue, bam_list in zip(tissue_names, bam_files):
            for rep_bam in bam_list:
                bam_path = os.path.join(args.bam_db, rep_bam)
                params = [bam_path, snp_locs, flanking]
                result = pool.apply_async(extract_cuts, params)
                jobs.append((tissue, rep_bam, result))

        for tissue, rep_bam, result in jobs:
            access_dict[tissue].append(result.get())
            logger.info("finished processing %s", rep_bam)

    with open(args.output_file, "wb") as pkl_handle:
        output = (access_dict, tissue_names, snp_md, lab_list)
        pickle.dump(output, pkl_handle)
    logger.info("successfully exported data matrix")
    logger.info("all done - good bye")


def extract_cuts(bam_path, snp_locs, flanking):
    d = len(snp_locs)
    cuts = np.zeros(d)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for snp_no, (chrom, pos) in enumerate(snp_locs):
            # here we have to shift to zero based indexing
            start = pos - 1 - flanking
            end = pos + flanking
            for read in bam.fetch(chrom, start, end):
                # cuts are always at the 5' end
                if read.is_reverse:
                    cut = read.reference_end
                else:
                    cut = read.reference_start
                if start <= cut < end:
                    cuts[snp_no] += 1
    return cuts


if __name__ == '__main__':
    main()
