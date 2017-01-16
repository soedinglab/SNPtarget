import argparse
import ctypes as ct
import bisect
import sys
from functools import total_ordering
import heapq
import logging
from multiprocessing import Process, JoinableQueue
from multiprocessing import sharedctypes as sct
import gc
import os
import random
import glob
from collections import namedtuple
import numpy as np
from scipy import stats

from ngsbiotools.ivtree import Interval, IVTree
from ngsbiotools.parsers import RegionParser
from ngsbiotools.misc import pearson_err
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
    parser.add_argument("--min_random_samples", type=int, default=10000,
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
    parser.add_argument("bg_stat_file")
    parser.add_argument("out_file", help="output file for predictions")
    parser.add_argument("access_db", help="path to the accessiblity database")
    parser.add_argument("access_error_db",
                        help="path to database with accessibility error estimates")
    parser.add_argument("expr_db", help="path to expression database")
    parser.add_argument("expr_error_db",
                        help="path to database with expression error estimates")
    add_logging_args(parser)

    parser.add_argument("--pure_cor", action="store_true")
    parser.add_argument("--no_error", action="store_true")
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

    logger.info("preparing sample database")
    all_chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

    acc_db = []
    err_db = []
    id_row_mapping = {}
    probing_samples = []
    random_samples = {}
    next_row = 0
    full_rows = []
    probe_rows = []

    for chrom in all_chroms:
        acc_fglob = os.path.join(args.access_db, chrom + ".*")
        acc_err_fglob = os.path.join(args.access_error_db, chrom + ".*")
        expr_fglob = os.path.join(args.expr_db, chrom + ".*")
        expr_err_fglob = os.path.join(args.expr_error_db, chrom + ".*")

        try:
            acc_fpath, = glob.glob(acc_fglob)
            acc_err_fpath, = glob.glob(acc_err_fglob)
            expr_fpath, = glob.glob(expr_fglob)
            expr_err_fpath, = glob.glob(expr_err_fglob)
        except ValueError:
            # not exactly one file matching
            logger.warn("cannot find unique database for chromosome %r", chrom)
            continue

        logger.debug("building the accessibility interval tree for chromosome %s",
                     chrom)
        sample_tree = IVTree()
        with open(acc_fpath) as acc_handle, open(acc_err_fpath) as err_handle:
            # discard header
            acc_handle.readline()
            err_handle.readline()

            for acc_line, err_line in zip(acc_handle, err_handle):
                acc_id, _, acc_pos, *access_strs = acc_line.split()
                err_id, _, err_pos, *error_strs = err_line.split()
                assert acc_id == err_id
                assert acc_pos == err_pos
                pos = int(acc_pos)
                access_vals = np.array(access_strs, dtype=float)
                access_vals -= access_vals.mean()
                error_vals = np.array(error_strs, dtype=float)

                if not np.isclose(access_vals.std(), 0):
                    iv = Interval(pos, pos)
                    iv.access_vals = access_vals
                    iv.access_errs = error_vals
                    iv.id = acc_id
                    sample_tree.insert(iv)

        if chrom == args.chrom:
            dhs_tree = sample_tree

        logger.debug("building the expression database for chromosome %s",
                     chrom)
        chrom_expr_db = {}
        with open(expr_fpath) as expr_handle, open(expr_err_fpath) as err_handle:
            # discard header
            expr_handle.readline()
            err_handle.readline()

            for expr_line, err_line in zip(expr_handle, err_handle):
                _, _, _, expr_gene, _, *expr_strs = expr_line.split()
                _, _, _, err_gene, _, *err_strs = err_line.split()
                assert(expr_gene == err_gene)
                expr_vals = np.array(expr_strs, dtype=float)
                expr_vals -= expr_vals.mean()
                err_vals = np.array(err_strs, dtype=float)
                if not np.isclose(expr_vals.std(), 0):
                    expr = Expression(expr_gene, expr_vals, err_vals)
                    chrom_expr_db[expr_gene] = expr

        if chrom == args.chrom:
            expr_db = chrom_expr_db

        logger.debug("creating the sampling database for chromosome %s",
                     chrom)

        with open(args.region_file) as reg_handle:
            probing_bytes = 0
            full_bytes = 0
            n_probe = args.probing_size
            reg_parser = RegionParser(reg_handle)
            for reg in reg_parser.parse():
                if (reg.chrom != chrom or
                        reg.gene_name not in chrom_expr_db):
                    continue
                iv = Interval(reg.start, reg.end)
                sample_dhs = list(sample_tree.query_all_overlaps(iv))
                n_full = len(sample_dhs)
                if n_full == 0:
                    continue

                full_ids = []
                for dhs in sample_dhs:
                    row = id_row_mapping.get(dhs.id)
                    if row is None:
                        id_row_mapping[dhs.id] = next_row
                        acc_db.append(dhs.access_vals)
                        err_db.append(dhs.access_errs)
                        row = next_row
                        next_row += 1
                    full_ids.append(row)
                full_ids = np.array(full_ids, dtype=np.int)
                probe_ids = np.random.choice(full_ids, n_probe)

                sh_full_ids = sct.RawArray(ct.c_int, len(full_ids))
                sh_full_ids[:] = full_ids
                full_bytes += ct.sizeof(sh_full_ids)
                full_rows.append(sh_full_ids)

                sh_probe_ids = sct.RawArray(ct.c_int, n_probe)
                sh_probe_ids[:] = probe_ids
                probing_bytes += ct.sizeof(sh_probe_ids)
                probe_rows.append(sh_probe_ids)

                tss = reg.start + (reg.end - reg.start) // 2
                probing_samples.append((reg.gene_name, chrom, tss, sh_probe_ids))
                random_samples[reg.gene_name] = sh_full_ids

        logger.info("finished sampling database for chromosome %s", chrom)

    logger.info("writing sampling database to shared memory representation")
    n_db = len(acc_db)
    d_db = len(acc_db[0])

    sh_acc_db = sct.RawArray(ct.c_double, n_db * d_db)
    sh_err_db = sct.RawArray(ct.c_double, n_db * d_db)
    acc_bytes = ct.sizeof(sh_acc_db)
    err_bytes = ct.sizeof(sh_err_db)

    cur_pointer = 0
    for i in range(n_db):
        sh_acc_db[cur_pointer:cur_pointer + d_db] = acc_db[i]
        sh_err_db[cur_pointer:cur_pointer + d_db] = err_db[i]
        cur_pointer += d_db
    del acc_db
    del err_db

    logger.info("probing database size: %.2f MiB", probing_bytes / 2**20)
    logger.info("background database size: %.2f MiB", full_bytes / 2**20)
    logger.info("access database size: %.2f MiB", acc_bytes / 2**20)
    logger.info("error database size: %.2f MiB", err_bytes / 2**20)
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
        gc.collect()
        dworker = DataWorker(reg_queue, data_queue, expr_db, dhs_tree)
        del dhs_tree
        del expr_db
        dworker.start()

        logger.info("started data worker")

        # calculating the p-values is highly parallelized
        gc.collect()
        for i in range(args.n_threads):
            cworker = CorrWorker(
                data_queue, out_queue, random_samples, probing_samples,
                n_db, d_db, sh_acc_db, sh_err_db,
                args.min_random_samples, args.min_tss_skip,
                args.min_bg_tss, args.ks_probe_threshold, args.ks_bg_threshold,
                args.ks_split, args.pure_cor, args.no_error
            )
            cworker.start()
            logger.info("started correlation worker %s", i)

        # create the output file writer process
        del random_samples
        del probing_samples
        del sh_acc_db
        del sh_err_db
        gc.collect()
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

    def __init__(self, reg_queue, data_queue, expr_db, dhs_tree):
        Process.__init__(self)
        self.daemon = True

        self.reg_queue = reg_queue
        self.data_queue = data_queue
        self.expr_db = expr_db
        self.dhs_tree = dhs_tree

    def run(self):
        while True:
            qmsg = self.reg_queue.get()
            reg = qmsg["region"]
            probe = self.expr_db[reg.gene_name]
            reg_iv = Interval(reg.start, reg.end)
            dhs = list(self.dhs_tree.query_all_overlaps(reg_iv))

            qans = {
                "region": reg,
                "probe": probe,
                "dhs": dhs,
            }
            self.data_queue.put(qans)
            self.reg_queue.task_done()


class CorrWorker(Process):

    def __init__(self, data_queue, out_queue, random_samples, probing_samples,
                 n_db, d_db, sh_acc_arr, sh_err_arr, min_n_samples, min_dist=1e7, min_n_tss=10,
                 ks_probe_thresh=0.2, ks_bg_thresh=0.1, ks_split=0.5, pure_cor=False, no_error=False):

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
        self.pure_cor = pure_cor
        self.no_error = no_error

        # convert shared memory types
        acc_db = np.ctypeslib.as_array(sh_acc_arr)
        acc_db.shape = (n_db, d_db)
        self.acc_db = acc_db

        err_db = np.ctypeslib.as_array(sh_err_arr)
        err_db.shape = (n_db, d_db)
        self.err_db = err_db

        self.probing_samples = probing_samples
        self.random_samples = random_samples

    def run(self):
        min_n_samples = self.min_n_samples
        min_dist = self.min_dist
        min_n_tss = self.min_n_tss
        ks_probe_thresh = self.ks_probe_thresh
        ks_bg_thresh = self.ks_bg_thresh
        ks_split = self.ks_split
        pure_cor = self.pure_cor
        no_error = self.no_error
        probing_samples = self.probing_samples
        random_samples = self.random_samples
        acc_db = self.acc_db
        err_db = self.err_db

        ks_pval_stack = []
        for i in range(20 * min_n_tss):
            ks_pval_stack.append(PvalWrapper(0, None))
        spare_wrapper = PvalWrapper(0, None)

        while True:
            qmsg = self.data_queue.get()
            skip = qmsg.get("skip", False)
            reg = qmsg["region"]
            logger.debug("start prediction of gene %r", reg.gene_name)
            # prepare auxiliary variables
            probe = qmsg["probe"]
            expr = probe.expr_vals
            expr_err = probe.expr_errs

            dhs_list = qmsg["dhs"]
            if len(dhs_list) == 0:
                skip = True
            reg_center = reg.start + (reg.end - reg.start) // 2
            reg_strand = reg.strand

            if skip:
                self.out_queue.put({"gene_name": reg.gene_name, "skip": True})
                self.data_queue.task_done()
                continue

            # calculate correlation distribution for predictions
            score_list = []
            for dhs in dhs_list:
                cor, r, dr = pearson_err(dhs.access_vals, dhs.access_errs, expr, expr_err)
                if no_error:
                    dr = 1
                score = r / dr if not pure_cor else cor
                score_list.append(score)
            scores = np.array(score_list)
            scores.sort()
            split_score = np.percentile(scores, q=ks_split * 100)
            split_score_index = np.searchsorted(scores, split_score)
            lower_scores = scores[:split_score_index]

            # do a full probe screening
            for wrapper in ks_pval_stack:
                wrapper.pval = 0
                wrapper.gene = None

            for gene_name, chrom, tss, probe_rows in probing_samples:

                # do not sample from regions that could contain true positives
                if chrom == reg.chrom and abs(tss - reg_center) < min_dist:
                    continue

                n = len(probe_rows)
                probe_scores = np.empty(n)
                for i, row in enumerate(probe_rows):
                    cor, r, dr = pearson_err(acc_db[row], err_db[row], expr, expr_err)
                    if no_error:
                        dr = 1
                    score = r / dr if not pure_cor else cor
                    probe_scores[i] = score
                probe_scores.sort()

                probe_split_index = np.searchsorted(probe_scores, split_score)
                lower_probe_scores = probe_scores[:probe_split_index]
                ks_stat, pval = stats.ks_2samp(lower_scores, lower_probe_scores)
                spare_wrapper.pval = pval
                spare_wrapper.gene = gene_name

                # push to a minimum heap which keeps the maximal observed
                # p-values
                spare_wrapper = heapq.heappushpop(ks_pval_stack, spare_wrapper)

            sampled = 0
            cur_samples = 0
            bg_model_scores = []
            bg_tss_names = []

            for wrapper in sorted(ks_pval_stack, reverse=True):

                if wrapper.pval < ks_probe_thresh or sampled >= min_n_tss:
                    break

                rows = random_samples[wrapper.gene]
                n = len(rows)

                bg_scores = np.empty(n)
                for i, row in enumerate(rows):
                    cor, r, dr = pearson_err(acc_db[row], err_db[row], expr, expr_err)
                    if no_error:
                        dr = 1
                    score = r / dr if not pure_cor else cor
                    bg_scores[i] = score
                bg_scores.sort()

                bg_split_index = np.searchsorted(bg_scores, split_score)
                lower_bg_scores = bg_scores[:bg_split_index]
                ks_stat, pval = stats.ks_2samp(lower_scores, lower_bg_scores)

                if pval < ks_bg_thresh:
                    continue

                bg_model_scores.append(bg_scores)
                bg_tss_names.append(wrapper.gene)
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
                cor, r, dr = pearson_err(dhs.access_vals, dhs.access_errs, expr, expr_err)
                if no_error:
                    dr = 1
                score = r / dr if not pure_cor else cor
                left_pos = bisect.bisect_left(bg_model, score)
                right_pos = bisect.bisect_right(bg_model, score)
                insert_pos = np.round(left_pos + (right_pos - left_pos) / 2)
                pval_e = 1 - insert_pos / cur_samples

                res_data = [
                    dhs.start, dhs.end, dhs.id,
                    probe.symbol, reg_strand, reg_center,
                    cor, score, pval_e
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


Expression = namedtuple("Expression", ["symbol", "expr_vals", "expr_errs"])
result_fields = [
    "dhs_start", "dhs_end", "dhs_id",
    "probe_symbol", "gene_strand", "gene_tss",
    "cor", "score", "pvalue"
]
Result = namedtuple("Result", result_fields)


class ResultStruct(ct.Structure):
    _fields_ = [
        ("rho", ct.c_double),
        ("r", ct.c_double),
        ("dr", ct.c_double),
    ]


class TupleWrapper(dict):

    def __init__(self, nt_dict):
        super().__init__(nt_dict)

    def __getattr__(self, attrib):
        return super().__getitem__(attrib)

    def __getstate__(self):
        return dict(self)

    def __setstate__(self, d):
        return TupleWrapper(d)


if __name__ == "__main__":
    main()
