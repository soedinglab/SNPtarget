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
    parser.add_argument('chrom')
    parser.add_argument('width', type=int)
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
                params = [bam_path, args.width, args.chrom]
                result = pool.apply_async(extract_cuts, params)
                jobs.append((tissue, rep_bam, result))

        for tissue, rep_bam, result in jobs:
            access_dict[tissue].append(result.get())
            logger.info("finished processing %s", rep_bam)

    max_start = 0
    for cut_list in access_dict.values():
        for cut_dict in cut_list:
            new_max_start = max(cut_dict.keys())
            max_start = max(max_start, new_max_start)

    metadata = []
    for i in range(max_start + 1):
        md = {}
        md['chrom'] = args.chrom
        md['start'] = i * args.width + 1
        md['end'] = (i + 1) * args.width
        metadata.append(md)

    for tissue, cut_list in access_dict.items():
        new_list = []
        for cut_dict in cut_list:
            cut_arr = np.zeros(max_start + 1)
            cut_arr[list(cut_dict.keys())] = list(cut_dict.values())
            new_list.append(cut_arr)
        access_dict[tissue] = new_list

    with open(args.output_file, 'wb') as pkl_handle:
        output = (access_dict, tissue_names, metadata, lab_list)
        pickle.dump(output, pkl_handle)
    logger.info("successfully exported data matrix")
    logger.info("all done - good bye")


def extract_cuts(bam_path, width, chrom):
    cuts = defaultdict(int)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for ali in bam.fetch(chrom):
            if ali.is_reverse:
                cut = ali.reference_end
            else:
                cut = ali.reference_start
            cuts[(cut + 1) // width] += 1
    return cuts


if __name__ == '__main__':
    main()
