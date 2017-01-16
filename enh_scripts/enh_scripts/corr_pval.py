import argparse
import ctypes as ct
import bisect
import sys
from functools import total_ordering
import logging
from multiprocessing import Process, JoinableQueue
from multiprocessing import sharedctypes as sct
import os
import random
import glob
from collections import namedtuple
import numpy as np
from scipy import stats

from ngsbiotools.ivtree import Interval, IVTree
from ngsbiotools.parsers import RegionParser
from ngsbiotools.misc.misc import centered_pearsonr as pearsonr
from enh_scripts import init_logger, add_logging_args


logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser(
        prog="corr_pval.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("region_file", help="region file for correlations")
    parser.add_argument("--seed", "-s", type=int, default=None,
                        help="seed of the random generator")
    parser.add_argument("chrom", help="chromosome for predictions")
    parser.add_argument("--bg_sample_size", "-k", type=int, default=10000,
                        help="sample size for p-value calculation")

    parser.add_argument("--scramble_count", type=int, default=10,
                        help="number of scrambled tss-dhs pairs for each gene")

    help_text = ("minimum distance between two TSS regions"
                 "on the same chromosome")
    parser.add_argument("--min_tss_skip", type=int, default=1e7,
                        help=help_text)

    help_text = "the minimum number of dhs in the background model"
    parser.add_argument("--min_bg_samples", type=int, default=10000,
                        help=help_text)

    help_text = "the p-value threshold in probing KS-Test statistic"
    parser.add_argument("--ks_probe_threshold", type=float, default=0.01,
                        help=help_text)

    help_text = "the p-value threshold in background KS-Test statistic"
    parser.add_argument("--ks_bg_threshold", type=float, default=0.01,
                        help=help_text)

    help_text = ("the minimum count of different tss sampled in the "
                 "background model")
    parser.add_argument("--min_bg_tss", type=int, default=10, help=help_text)

    help_text = "the relative splitting position for the KS-Test"
    parser.add_argument("--ks_split", type=float, default=0.5,
                        help=help_text)

    help_text = "the sample size used for random probing"
    parser.add_argument("--probing_size", type=int, default=100,
                        help=help_text)
    parser.add_argument("--n_threads", type=int, default=3,
                        help="number of worker threads")
    parser.add_argument("--dist_ks_thresh", type=float, default=1e-3)
    parser.add_argument("--high_dist_thresh", type=float, default=1.5e6)
    parser.add_argument("--min_n_high_dist", type=int, default=50)
    parser.add_argument('--tss_min_dist', type=int, default=0)

    parser.add_argument("bg_stat_file")
    parser.add_argument("out_file", help="output file for predictions")
    parser.add_argument("access_db", help="path to the accessiblity database")
    parser.add_argument("expr_db", help="path to expression database")
    add_logging_args(parser)

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()
    init_logger(args)

    logger.info("invoked script via %r", " ".join(sys.argv))
    np.random.seed(args.seed)
    random.seed(args.seed)

    chrom_e = args.chrom
    chrom_g = args.chrom
    tss_min_dist = args.tss_min_dist

    logger.info("preparing sample database")
    all_chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

    acc_db = []
    id_row_mapping = {}
    probing_samples = []
    random_samples = {}
    next_row = 0
    full_rows = []
    probe_rows = []

    for chrom in all_chroms:
        acc_fglob = os.path.join(args.access_db, chrom + ".*")
        expr_fglob = os.path.join(args.expr_db, chrom + ".*")

        try:
            acc_fpath, = glob.glob(acc_fglob)
            expr_fpath, = glob.glob(expr_fglob)
        except ValueError:
            # not exactly one file matching
            logger.warn("cannot find unique database for chromosome %r", chrom)
            continue

        logger.debug("building the accessibility interval tree for chromosome %s",
                     chrom)
        sample_tree = IVTree()
        with open(acc_fpath) as acc_handle:
            # discard header
            acc_handle.readline()

            for acc_line in acc_handle:
                acc_id, _, acc_pos, *access_strs = acc_line.split()
                pos = int(acc_pos)
                access_vals = np.array(access_strs, dtype=float)
                access_vals -= access_vals.mean()

                if not np.isclose(access_vals.std(), 0):
                    iv = Interval(pos, pos)
                    iv.access_vals = access_vals
                    iv.id = acc_id
                    sample_tree.insert(iv)

        if chrom == args.chrom:
            dhs_tree = sample_tree

        logger.debug("building the expression database for chromosome %s",
                     chrom)
        chrom_expr_db = {}
        with open(expr_fpath) as expr_handle:
            # discard header
            expr_handle.readline()

            for expr_line in expr_handle:
                _, _, _, expr_gene, _, *expr_strs = expr_line.split()
                expr_vals = np.array(expr_strs, dtype=float)
                expr_vals -= expr_vals.mean()
                if not np.isclose(expr_vals.std(), 0):
                    expr = Expression(expr_gene, expr_vals)
                    chrom_expr_db[expr_gene] = expr

        if chrom == args.chrom:
            expr_db = chrom_expr_db

        logger.debug("creating the sampling database for chromosome %s",
                     chrom)

        with open(args.region_file) as reg_handle:
            n_probe = args.probing_size
            reg_parser = RegionParser(reg_handle)
            for reg in reg_parser.parse():
                if (reg.chrom != chrom or
                        reg.gene_name not in chrom_expr_db):
                    continue
                iv = Interval(reg.start, reg.end)
                sample_dhs = list(sample_tree.query_all_overlaps(iv))
                sample_dhs = min_tss_dist_filter(sample_dhs, reg, tss_min_dist)
                n_full = len(sample_dhs)

                # for the sampling procedure it is important to have enough sites per
                # tss. Here I demand that the number of dhs should at least be as large
                # as the number of dhs used in the probing procedure
                if n_full < args.probing_size:
                    continue

                full_ids = []
                for dhs in sample_dhs:
                    row = id_row_mapping.get(dhs.id)
                    if row is None:
                        id_row_mapping[dhs.id] = next_row
                        acc_db.append(dhs.access_vals)
                        row = next_row
                        next_row += 1
                    full_ids.append(row)
                full_ids = np.array(full_ids, dtype=np.int)
                probe_ids = np.random.choice(full_ids, n_probe)

                sh_full_ids = sct.RawArray(ct.c_int, len(full_ids))
                sh_full_ids[:] = full_ids
                full_rows.append(sh_full_ids)

                sh_probe_ids = sct.RawArray(ct.c_int, n_probe)
                sh_probe_ids[:] = probe_ids
                probe_rows.append(sh_probe_ids)

                tss = reg.start + (reg.end - reg.start) // 2
                probing_samples.append((reg.gene_name, chrom, tss, sh_probe_ids))
                random_samples[reg.gene_name] = sh_full_ids

        logger.info("finished sampling database for chromosome %s", chrom)

    logger.info("writing sampling database to shared memory representation")
    n_db = len(acc_db)
    d_db = len(acc_db[0])

    sh_acc_db = sct.RawArray(ct.c_double, n_db * d_db)

    cur_pointer = 0
    for i in range(n_db):
        sh_acc_db[cur_pointer:cur_pointer + d_db] = acc_db[i]
        cur_pointer += d_db
    del acc_db

    logger.info("successfully created the complete sampling database")
    logger.info("there are %s tss in the sampling database",
                len(random_samples))

    logger.info("reading the regions")
    with open(args.region_file) as reg_handle:
        reg_queue = JoinableQueue()
        reg_count = 0
        reg_parser = RegionParser(reg_handle)
        for reg in reg_parser.parse():
            if reg.chrom == args.chrom and reg.gene_name in expr_db:
                reg_queue.put({"region": TupleWrapper(reg._asdict())})
                reg_count += 1

    with open(args.out_file, "w") as out_handle,\
            open(args.bg_stat_file, "w") as bg_stat_handle:
        data_queue = JoinableQueue(32)
        out_queue = JoinableQueue(32)

        # start the worker subprocesses
        dworker = DataWorker(reg_queue, data_queue, expr_db, dhs_tree, tss_min_dist)
        del dhs_tree
        del expr_db
        dworker.start()

        logger.info("started data worker")

        # calculating the p-values is highly parallelized
        for i in range(args.n_threads):
            cworker = CorrWorker(
                data_queue, out_queue, random_samples, probing_samples,
                n_db, d_db, sh_acc_db,
                args.min_bg_samples, args.probing_size, args.min_tss_skip,
                args.min_bg_tss, args.ks_probe_threshold, args.ks_bg_threshold,
                args.ks_split,
                args.dist_ks_thresh, args.high_dist_thresh, args.min_n_high_dist,
            )
            cworker.start()
            logger.info("started correlation worker %s", i)

        # create the output file writer process
        del random_samples
        del probing_samples
        del sh_acc_db
        wworker = WriterWorker(
            out_handle, bg_stat_handle, out_queue, chrom_e, chrom_g,
            reg_count
        )
        wworker.start()

        reg_queue.join()
        data_queue.join()
        out_queue.join()

        logger.info("finished predictions. Exiting.")


class DataWorker(Process):

    def __init__(self, reg_queue, data_queue, expr_db, dhs_tree, min_dist):
        Process.__init__(self)
        self.daemon = True

        self.reg_queue = reg_queue
        self.data_queue = data_queue
        self.expr_db = expr_db
        self.dhs_tree = dhs_tree
        self.min_dist = min_dist

    def run(self):
        min_dist = self.min_dist
        while True:
            qmsg = self.reg_queue.get()
            reg = qmsg["region"]
            probe = self.expr_db[reg.gene_name]
            reg_iv = Interval(reg.start, reg.end)
            dhs = list(self.dhs_tree.query_all_overlaps(reg_iv))
            dhs = min_tss_dist_filter(dhs, reg, min_dist)

            dists = np.empty(len(dhs))
            r_mid = reg.start + (reg.end - reg.start) // 2
            for i, d in enumerate(dhs):
                d_mid = d.start + (d.end - d.start) // 2
                dists[i] = d_mid - r_mid
            dists = np.abs(dists)

            qans = {
                "region": reg,
                "probe": probe,
                "dhs": dhs,
                "distances": dists,
            }
            self.data_queue.put(qans)
            self.reg_queue.task_done()


class CorrWorker(Process):

    # TODO this is getting very ugly...
    def __init__(self, data_queue, out_queue, random_samples, probing_samples,
                 n_db, d_db, sh_acc_arr, min_n_samples, min_dhs_per_tss,
                 min_dist, min_n_tss, ks_probe_thresh, ks_bg_thresh,
                 ks_split, dist_ks_thresh, high_dist_thresh, min_n_high_dist):

        Process.__init__(self)
        self.daemon = True

        self.data_queue = data_queue
        self.out_queue = out_queue
        self.min_n_samples = min_n_samples
        self.min_dist = min_dist
        self.min_n_tss = min_n_tss
        self.ks_probe_thresh = ks_probe_thresh
        self.ks_bg_thresh = ks_bg_thresh
        self.ks_split = ks_split
        self.min_dhs_per_tss = min_dhs_per_tss
        self.dist_ks_thresh = dist_ks_thresh
        self.high_dist_thresh = high_dist_thresh
        self.min_n_high_dist = min_n_high_dist

        # convert shared memory types
        acc_db = np.ctypeslib.as_array(sh_acc_arr)
        acc_db.shape = (n_db, d_db)
        self.acc_db = acc_db

        self.probing_samples = probing_samples
        self.random_samples = random_samples

    def run(self):
        min_n_samples = self.min_n_samples
        min_dist = self.min_dist
        min_n_tss = self.min_n_tss
        ks_probe_thresh = self.ks_probe_thresh
        ks_bg_thresh = self.ks_bg_thresh
        ks_split = self.ks_split
        probing_samples = self.probing_samples
        random_samples = self.random_samples
        acc_db = self.acc_db
        min_dhs_per_tss = self.min_dhs_per_tss

        dist_ks_thresh = self.dist_ks_thresh
        high_dist_thresh = self.high_dist_thresh
        min_n_high_dist = self.min_n_high_dist

        gene_names = []
        for gene_name, chrom, tss, probe_rows in probing_samples:
            gene_names.append(gene_name)
        gene_vec = np.array(gene_names)

        while True:
            qmsg = self.data_queue.get()
            skip = qmsg.get("skip", False)
            reg = qmsg["region"]
            logger.debug("start prediction of gene %r", reg.gene_name)
            # prepare auxiliary variables
            probe = qmsg["probe"]
            expr = probe.expr_vals
            dists = qmsg["distances"]
            dhs_list = qmsg["dhs"]
            if len(dhs_list) < min_dhs_per_tss:
                skip = True
            reg_center = reg.start + (reg.end - reg.start) // 2
            reg_strand = reg.strand

            if not skip:
                # calculate correlation distribution for predictions
                score_list = []

                for dhs in dhs_list:
                    score = pearsonr(dhs.access_vals, expr)
                    score_list.append(score)
                scores = np.array(score_list)

                high_dist_scores = scores[dists > high_dist_thresh]
                n_high_dists = len(high_dist_scores)

                scores.sort()
                split_score = np.percentile(scores, q=ks_split * 100)
                split_score_index = np.searchsorted(scores, split_score)
                lower_scores = scores[:split_score_index]

                if n_high_dists < min_n_high_dist:
                    msg = "not enough high distance predictions (%d/%d) for '%s'"
                    logger.debug(msg, n_high_dists, min_n_high_dist, reg.gene_name)
                    # not enough high distance interaction pairs
                    skip = True

            if skip:
                self.out_queue.put({"gene_name": reg.gene_name, "skip": True})
                self.data_queue.task_done()
                continue

            ks_pval_vec = -np.ones(len(probing_samples))
            for probe_no, (gene_name, chrom, tss, probe_rows) in enumerate(probing_samples):

                # do not sample from regions that could contain true positives
                if chrom == reg.chrom and abs(tss - reg_center) < min_dist:
                    continue

                n = len(probe_rows)
                probe_scores = np.empty(n)
                for i, row in enumerate(probe_rows):
                    score = pearsonr(acc_db[row], expr)
                    probe_scores[i] = score
                probe_scores.sort()

                probe_split_index = np.searchsorted(probe_scores, split_score)
                lower_probe_scores = probe_scores[:probe_split_index]
                ks_stat, pval = stats.ks_2samp(lower_scores, lower_probe_scores)
                ks_pval_vec[probe_no] = pval

            sampled = 0
            cur_samples = 0
            bg_model_scores = []
            bg_tss_names = []

            ks_probe_mask = ks_pval_vec >= ks_probe_thresh
            shuffled_genes = np.random.permutation(gene_vec[ks_probe_mask])
            for gene in shuffled_genes:

                if sampled >= min_n_tss:
                    break

                rows = random_samples[gene]
                n = len(rows)

                bg_scores = np.empty(n)
                for i, row in enumerate(rows):
                    score = pearsonr(acc_db[row], expr)
                    bg_scores[i] = score
                bg_scores.sort()

                # check if the background model fits for low correlation p-values
                bg_split_index = np.searchsorted(bg_scores, split_score)
                lower_bg_scores = bg_scores[:bg_split_index]
                ks_stat, pval = stats.ks_2samp(lower_scores, lower_bg_scores)

                if pval < ks_bg_thresh:
                    continue

                # check if the background model fits for high distances
                high_dist_pvals = np.empty(n_high_dists)
                for i, score in enumerate(high_dist_scores):
                    left_pos = bisect.bisect_left(bg_scores, score)
                    right_pos = bisect.bisect_right(bg_scores, score)
                    insert_pos = np.round(left_pos + (right_pos - left_pos) / 2)
                    pval = 1 - insert_pos / n
                    high_dist_pvals[i] = pval

                # since there are virtually no true positives at very high distances, we
                # assume that the p-value distribution is uniform
                dist_ks_pval = stats.kstest(high_dist_pvals, stats.uniform.cdf).pvalue

                if dist_ks_pval < dist_ks_thresh:
                    continue

                bg_model_scores.append(bg_scores)
                bg_tss_names.append(gene)
                cur_samples += len(bg_scores)
                sampled += 1

            # we did not manage to draw a proper background model
            if cur_samples < min_n_samples:
                msg = "insufficient sample size (%s/%s) from %s tss"
                logger.debug(msg, cur_samples, min_n_samples, sampled)
                qans = {
                    "gene_name": reg.gene_name,
                    "skip": True,
                }
                self.out_queue.put(qans)
                self.data_queue.task_done()
                continue

            else:
                bg_model = np.concatenate(bg_model_scores)
                bg_model.sort()
                msg = "successfully collected %s samples from %s tss"
                logger.debug(msg, cur_samples, sampled)

            results = []
            for dhs in dhs_list:
                score = pearsonr(dhs.access_vals, expr)
                left_pos = bisect.bisect_left(bg_model, score)
                right_pos = bisect.bisect_right(bg_model, score)
                insert_pos = np.round(left_pos + (right_pos - left_pos) / 2)
                pval_e = 1 - insert_pos / cur_samples

                res_data = [
                    dhs.start, dhs.end, dhs.id,
                    probe.symbol, reg_strand, reg_center,
                    score, score, pval_e
                ]
                res = Result(*res_data)
                results.append(res)

            # calculate statistics on the background model
            pred_count = len(scores)
            bg_stats = {}
            bg_stats["bg_mean"] = np.mean(bg_model)
            bg_stats["bg_median"] = np.median(bg_model)
            bg_stats["bg_std"] = np.std(bg_model)
            bg_stats["bg_n"] = cur_samples
            bg_stats["bg_tss"] = len(bg_model_scores)
            bg_stats["bg_tss_names"] = bg_tss_names
            bg_stats["pred_mean"] = np.mean(scores)
            bg_stats["pred_median"] = np.median(scores)
            bg_stats["pred_std"] = np.std(scores)
            bg_stats["pred_n"] = pred_count

            pred_distr = scores[:split_score_index]
            bg_distr_split_index = np.searchsorted(bg_model, split_score)
            bg_distr = bg_model[:bg_distr_split_index]
            ks_stat, ks_pval = stats.ks_2samp(pred_distr, bg_distr)
            bg_stats["ks_stat"] = ks_stat
            bg_stats["ks_pval"] = ks_pval

            qans = {
                "gene_name": reg.gene_name,
                "predictions": results,
                "bg_stats": bg_stats,
                "skip": False,
            }
            self.out_queue.put(qans)

            self.data_queue.task_done()
            logger.debug("finished prediction on gene %r", reg.gene_name)


class WriterWorker(Process):

    def __init__(self, fit_handle, bg_stat_handle, out_queue, chrom_e, chrom_g,
                 n_total):
        Process.__init__(self)
        self.fit_handle = fit_handle
        self.bg_stat_handle = bg_stat_handle
        self.out_queue = out_queue
        self.chrom_e = chrom_e
        self.chrom_g = chrom_g
        self.n_total = n_total
        self.daemon = True

    def run(self):

        # prepare variables
        n_finished = 0
        fit_handle = self.fit_handle
        bgstat_handle = self.bg_stat_handle

        # print headers
        print(
            "DHS", "DHS_id", "gene", "gene_symb", "cor", "score", "pval_e",
            sep="\t", file=fit_handle
        )
        print(
            "gene", "bg_median", "bg_mean", "bg_std", "bg_n", "bg_tss",
            "pred_median", "pred_mean", "pred_std", "pred_n",
            "ks_stat", "ks_pval", "sampled_tss",
            sep="\t", file=bgstat_handle
        )

        while True:
            qmsg = self.out_queue.get()
            skip = qmsg.get("skip", False)
            gene = qmsg["gene_name"]

            n_finished += 1
            if not skip:
                results = qmsg["predictions"]
                for result in results:
                    d_start, d_end = result.dhs_start, result.dhs_end,
                    d_id = result.dhs_id
                    p_symb = result.probe_symbol
                    reg_strand, reg_center = result.gene_strand, result.gene_tss
                    cor, score, pval_e = result.cor, result.score, result.pvalue
                    print(
                        "{}|{}|{}".format(self.chrom_e, d_start, d_end),
                        d_id,
                        "{}|{}|{}|{}".format(p_symb, self.chrom_g,
                                             reg_strand, reg_center),
                        p_symb,
                        "{:.5f}".format(cor),
                        "{:.5f}".format(score),
                        pval_e,
                        sep="\t", file=fit_handle
                    )

                sts = qmsg["bg_stats"]

                print(
                    gene,
                    sts["bg_median"], sts["bg_mean"], sts["bg_std"],
                    sts["bg_n"], sts["bg_tss"],
                    sts["pred_median"], sts["pred_mean"], sts["pred_std"],
                    sts["pred_n"],
                    sts["ks_stat"], sts["ks_pval"],
                    "|".join(sts["bg_tss_names"]),
                    sep="\t", file=bgstat_handle
                )
                logger.info("finished gene '%s' (%i/%i)",
                            gene, n_finished, self.n_total)
            else:
                logger.info("skipped gene '%s' (%i/%i)",
                            gene, n_finished, self.n_total)

            fit_handle.flush()
            bgstat_handle.flush()
            self.out_queue.task_done()


@total_ordering
class PvalWrapper:
    def __init__(self, pval, gene):
        self.pval = pval
        self.gene = gene

    def __eq__(self, other):
        return self.pval == other.pval

    def __lt__(self, other):
        return self.pval < other.pval


Expression = namedtuple("Expression", ["symbol", "expr_vals"])
result_fields = [
    "dhs_start", "dhs_end", "dhs_id",
    "probe_symbol", "gene_strand", "gene_tss",
    "cor", "score", "pvalue"
]
Result = namedtuple("Result", result_fields)


class TupleWrapper(dict):

    def __init__(self, nt_dict):
        super().__init__(nt_dict)

    def __getattr__(self, attrib):
        return super().__getitem__(attrib)

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, d):
        return TupleWrapper(d)


def min_tss_dist_filter(dhs_list, region, min_dist):
    reg_center = region.start + (region.end - region.start) // 2
    out_list = []
    for dhs in dhs_list:
        dhs_center = dhs.start + (dhs.end - dhs.start) // 2
        if abs(dhs_center - reg_center) >= min_dist:
            out_list.append(dhs)
    return out_list


if __name__ == "__main__":
    main()
