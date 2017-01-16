import argparse
import pickle
import logging
import sys
import os
import glob
from collections import defaultdict

import numpy as np

from enh_scripts import init_logger, add_logging_args
from enh_scripts import HUMAN_CHROMS
logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser("normalize_cuts")
    parser.add_argument("cut_db")
    parser.add_argument("access_db")
    parser.add_argument("error_db")
    parser.add_argument("--batch_correction", action="store_true")
    parser.add_argument("--perc_cutoff", type=float, default=99)
    parser.add_argument("--dump_raw_data")
    parser.add_argument("--max_d_dump", type=int)
    parser.add_argument("--dump_proc_data")
    parser.add_argument("--pseudo_counts", type=int, default=0)
    add_logging_args(parser)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info("invoked script via %r", " ".join(sys.argv))

    if not os.path.exists(args.access_db):
        os.makedirs(args.access_db)
    if not os.path.exists(args.error_db):
        os.makedirs(args.error_db)

    def cutdata_generator():
        for chrom in HUMAN_CHROMS:
            cut_glob = os.path.join(args.cut_db, "%s_*" % chrom)
            cut_db_files = glob.glob(cut_glob)
            if len(cut_db_files) != 1:
                logger.warn("no or ambiguous cut data for chromosome %s. "
                            "Skipping.", chrom)
                continue
            cut_file, = cut_db_files
            with open(cut_file, "rb") as pkl_file:
                # access_dict, tissue_names, snp_md, batch_list
                yield (chrom,) + pickle.load(pkl_file)

    cur_d = 0
    chrom_slices = {}
    for chrom, access_dict, tissue_names, snp_md, batch_list in cutdata_generator():
        slice_d = len(snp_md)
        chrom_slices[chrom] = slice(cur_d, cur_d + slice_d)
        cur_d += slice_d

    row_reps_dict = {}
    cur_n = 0
    for tissue_name in tissue_names:
        rep_count = len(access_dict[tissue_name])
        row_reps_dict[tissue_name] = list(range(cur_n, cur_n + rep_count))
        cur_n += rep_count

    logger.info("assembling cut matrix")
    n, d = cur_n, cur_d
    access_matrix = np.empty((n, d))

    snp_md_db = {}
    for chrom, access_dict, tissue_names, snp_md, _ in cutdata_generator():
        chrom_slice = chrom_slices[chrom]
        snp_md_db[chrom] = snp_md
        for tissue_name in tissue_names:
            rep_rows = row_reps_dict[tissue_name]
            reps = access_dict[tissue_name]
            for row, rep in zip(rep_rows, reps):
                access_matrix[row, chrom_slice] = rep
    logger.info("successfully initialized cut matrix")

    if args.dump_raw_data:
        logger.info("writing raw cut data to binary pickle")
        rep_dict = {}
        for tissue_name in tissue_names:
            rep_rows = row_reps_dict[tissue_name]
            if not args.max_d_dump:
                rep_dict[tissue_name] = [access_matrix[r] for r in rep_rows]
            else:
                max_d = args.max_d_dump
                rep_dict[tissue_name] = [access_matrix[r, :max_d] for r in rep_rows]
        with open(args.dump_raw_data, "wb") as handle:
            pkl_data = rep_dict
            pickle.dump(pkl_data, handle)

    access_matrix += args.pseudo_counts

    # start normalization procedure
    logger.info("started normalization procedure")
    n, d = access_matrix.shape

    batch_map = defaultdict(list)
    for row, batch in enumerate(batch_list):
        batch_map[batch].append(row)

    # quantile normalization
    mean_arr = np.zeros(d)
    for row in access_matrix:
        mean_arr += np.sort(row)
    mean_arr /= n

    for i in range(n):
        access_matrix[i, np.argsort(access_matrix[i])] = mean_arr

    # cap at given percentile
    pc = args.perc_cutoff
    cutoff = np.percentile(mean_arr, pc)
    access_matrix[access_matrix > cutoff] = cutoff
    logger.info("finished quantile normalization")

    if args.batch_correction:
        # correct batch effects
        for rows in batch_map.values():
            corrector = np.mean(np.log(access_matrix[rows, :]), axis=0)
            for row in rows:
                access_matrix[row, :] = np.exp(np.log(access_matrix[row, :]) - corrector)
        logger.info("finished batch correction")

    min_cuts = np.min(access_matrix)
    max_cuts = np.max(access_matrix)
    access_matrix -= min_cuts
    access_matrix /= (max_cuts - min_cuts)
    logger.info("finished rescaling between [0, 1]")

    if args.dump_proc_data:
        logger.info("writing data to binary pickle")
        rep_dict = {}
        for tissue_name in tissue_names:
            rep_rows = row_reps_dict[tissue_name]
            if not args.max_d_dump:
                rep_dict[tissue_name] = [access_matrix[r] for r in rep_rows]
            else:
                max_d = args.max_d_dump
                rep_dict[tissue_name] = [access_matrix[r, :max_d] for r in rep_rows]
        with open(args.dump_proc_data, "wb") as handle:
            pkl_data = rep_dict
            pickle.dump(pkl_data, handle)

    logger.info("writing merged batch output files")
    for chrom in HUMAN_CHROMS:
        if chrom not in chrom_slices:
            continue

        chrom_slice = chrom_slices[chrom]
        n, d = len(tissue_names), len(snp_md_db[chrom])
        norm_acc_mat = np.empty((n, d))
        for row, tissue_name in enumerate(tissue_names):
            rep_rows = row_reps_dict[tissue_name]
            rep_block = access_matrix[rep_rows, chrom_slice]
            rep_mean = np.mean(rep_block, axis=0)
            norm_acc_mat[row] = rep_mean

        outfile = os.path.join(args.access_db, "%s.norm_access" % chrom)
        with open(outfile, "w") as out_handle:
            print("id", "chr", "pos", *tissue_names, sep=" ", file=out_handle)
            for snp_no, snp_md in enumerate(snp_md_db[chrom]):
                try:
                    snp_id = snp_md["id"]
                    chrom = snp_md["chrom"]
                    loc = snp_md["loc"]
                    print(snp_id, chrom, loc, *norm_acc_mat[:, snp_no],
                          sep=" ", file=out_handle)
                except KeyError:
                    chrom = snp_md['chrom']
                    start = snp_md['start']
                    end = snp_md['end']
                    mid = start + (end - start) // 2
                    reg_id = '|'.join(str(s) for s in (chrom, start, end))
                    print(reg_id, chrom, mid, *norm_acc_mat[:, snp_no],
                          sep=' ', file=out_handle)

    logger.info("writing estimated error output file")
    for chrom in HUMAN_CHROMS:
        if chrom not in chrom_slices:
            continue

        chrom_slice = chrom_slices[chrom]
        n, d = len(tissue_names), len(snp_md_db[chrom])
        err_acc_mat = np.empty((n, d))
        for row, tissue in enumerate(tissue_names):
            rep_rows = row_reps_dict[tissue]
            k = len(rep_rows)
            norm_rep = access_matrix[rep_rows, chrom_slice]
            acc_avg = np.mean(norm_rep, axis=0)
            err_rep = np.sum((norm_rep - acc_avg)**2, axis=0)
            if k > 1:
                err_rep /= (k - 1)
            err_acc_mat[row] = np.sqrt(err_rep)

        outfile = os.path.join(args.error_db, "%s.error" % chrom)
        with open(outfile, "w") as out_handle:
            print("id", "chr", "pos", *tissue_names, sep=" ", file=out_handle)
            assert len(snp_md_db[chrom]) == d
            for snp_no, snp_md in enumerate(snp_md_db[chrom]):
                # really ugly hack to support two different formats
                # should be fixed at some point
                try:
                    snp_id = snp_md["id"]
                    chrom = snp_md["chrom"]
                    loc = snp_md["loc"]
                    print(snp_id, chrom, loc, *err_acc_mat[:, snp_no],
                          sep=" ", file=out_handle)
                except KeyError:
                    chrom = snp_md['chrom']
                    start = snp_md['start']
                    end = snp_md['end']
                    mid = start + (end - start) // 2
                    reg_id = '|'.join(str(s) for s in (chrom, start, end))
                    print(reg_id, chrom, mid, *err_acc_mat[:, snp_no],
                          sep=' ', file=out_handle)

    logger.info("all done. exiting")


if __name__ == '__main__':
    main()
