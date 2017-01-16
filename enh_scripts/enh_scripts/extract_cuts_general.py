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
    parser.add_argument("region_file")
    parser.add_argument("data_info_file")
    parser.add_argument("output_file")
    parser.add_argument("--bam_db", default=".")
    parser.add_argument("--n_proc", type=int, default=1)
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info("invoked script via %r", " ".join(sys.argv))

    metadata = []
    regions = []
    logger.info("started reading regions")
    with open(args.region_file) as region_handle:
        for i, line in enumerate(region_handle):
            seqid, start_str, end_str, *_ = line.split()
            start = int(start_str)
            end = int(end_str)
            md = {}
            md['chrom'] = seqid
            md['start'] = start
            md['end'] = end
            metadata.append(md)
            regions.append((seqid, start, end))

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
                params = [bam_path, regions]
                result = pool.apply_async(extract_cuts, params)
                jobs.append((tissue, rep_bam, result))

        for tissue, rep_bam, result in jobs:
            access_dict[tissue].append(result.get())
            logger.info("finished processing %s", rep_bam)

    with open(args.output_file, "wb") as pkl_handle:
        output = (access_dict, tissue_names, metadata, lab_list)
        pickle.dump(output, pkl_handle)
    logger.info("successfully exported data matrix")
    logger.info("all done - good bye")


def extract_cuts(bam_path, regions):
    d = len(regions)
    cuts = np.zeros(d)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for region_id, (chrom, start, end) in enumerate(regions):
            # here we have to shift to zero based indexing
            start = start - 1
            for read in bam.fetch(chrom, start, end):
                # cuts are always at the 5' end
                if read.is_reverse:
                    cut = read.reference_end
                else:
                    cut = read.reference_start
                if start <= cut < end:
                    cuts[region_id] += 1
    return cuts


if __name__ == '__main__':
    main()
