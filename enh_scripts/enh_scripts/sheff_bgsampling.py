import argparse
import bisect
import sys
from functools import total_ordering
import heapq
import logging
from multiprocessing import Process, JoinableQueue
import os
import random

import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt

from ngsbiotools.ivtree import Interval, IVTree
from ngsbiotools.parsers import DHS_Parser, Expr_Parser, RegionParser
from ngsbiotools.misc import pearsonr
from ngsbiotools.parsers.expr_parser import Expr
from enh_scripts import init_logger, add_logging_args

logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser(
        prog="sheffield_bg",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("region_file", help="region file for correlations")
    parser.add_argument("--seed", "-s", type=int, default=None,
                        help="seed of the random generator")
    parser.add_argument("chrom", help="chromosome for predictions")
    parser.add_argument("--bg_sample_size", "-k", type=int, default=10000,
                        help="sample size for p-value calculation")

    parser.add_argument("--random_mode", action="store_true",
                        help="scramble tss and corresponding dhs")
    parser.add_argument("--scramble_count", type=int, default=10,
                        help="number of scrambled tss-dhs pairs for each gene")

    help_text = ("minimum distance between two TSS regions"
                 "on the same chromosome")
    parser.add_argument("--min_tss_skip", type=int, default=1e7,
                        help=help_text)

    help_text = "the minimum number of dhs in the background model"
    parser.add_argument("--min_random_samples", type=int, default=10000,
                        help=help_text)

    help_text = ("the maximum number of sampling draws for the "
                 "background model. Deactivated when non-positive.")
    parser.add_argument("--max_random_draws", type=int, default=100,
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

    parser.add_argument("--cor_method", choices=["pearson", "spearman"],
                        default="pearson", help="correlation implementation")
    parser.add_argument("--n_threads", type=int, default=3,
                        help="number of worker threads")
    parser.add_argument("--bg_stat_file")
    parser.add_argument("out_file", help="output file for predictions")
    parser.add_argument("dhs_db", help="path to DHS database")
    parser.add_argument("expr_db", help="path to expression database")

    help_text = "sample and predict on the same chromosome"
    parser.add_argument("--intrachrom", action="store_true", help=help_text)

    help_text = "the directory for the prediction vs background qq-plots"
    default_dir = os.path.join(os.curdir, "qqplots")
    parser.add_argument("--qq_dir", default=default_dir, help=help_text)

    help_text = "create qq-plots for all predictions"
    parser.add_argument("--qq_plotting", action="store_true", help=help_text)
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

    logger.info("preparing sample database")

    all_chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    if args.intrachrom:
        all_chroms = [args.chrom]

    probing_samples = []
    random_db = []
    random_samples = {}

    for chrom in all_chroms:
        logger.info("building the DHS interval tree for chromosome %s",
                    chrom)
        sample_tree = IVTree()
        sample_fpath = os.path.join(args.dhs_db, chrom + ".dhs.gz")
        with DHS_Parser(sample_fpath) as dhs_parser:
            for dhs in dhs_parser.read():
                if not np.isclose(dhs.access_vals.std(), 0):
                    sample_tree.insert(dhs)

        if chrom == args.chrom:
            dhs_tree = sample_tree

        logger.info("building the expression database for chromosome %s",
                    chrom)
        ran_expr_db = {}
        ran_expr_fpath = os.path.join(args.expr_db, chrom + ".expr")
        with Expr_Parser(ran_expr_fpath) as expr_parser:
            for expr in expr_parser.read():
                if not np.isclose(expr.expr_vals.std(), 0):
                    ran_expr_db[expr.symbol] = expr

        logger.info("creating the sampling database for chromosome %s",
                    chrom)

        with open(args.region_file) as reg_handle:
            reg_parser = RegionParser(reg_handle)
            for reg in reg_parser.parse():
                if (reg.chrom != chrom or
                        reg.gene_name not in ran_expr_db):
                    continue
                iv = Interval(reg.start, reg.end)
                sample_dhs = list(sample_tree.query_overlapping(iv))
                full_sample = [d.access_vals for d in sample_dhs]
                if len(full_sample) == 0:
                    continue
                probe_sample = []
                for i in range(args.probing_size):
                    probe_sample.append(random.choice(full_sample))
                tss = reg.start + (reg.end - reg.start) // 2
                probing_sample = (reg.gene_name, chrom, tss, probe_sample)
                probing_samples.append(probing_sample)
                random_samples[reg.gene_name] = full_sample

                if args.random_mode:
                    if reg.chrom == args.chrom:
                        continue
                    expr_probe = ran_expr_db[reg.gene_name]
                    random_db.append((reg, expr_probe, sample_dhs))

        logger.info("finished sampling database for chromosome %s", chrom)

    logger.info("successfully created the complete sampling database")
    logger.info("there are %s tss in the sampling database",
                len(random_samples))

    logger.info("creating the expression database")
    expr_db = {}
    expr_fpath = os.path.join(args.expr_db, args.chrom + ".expr")
    with Expr_Parser(expr_fpath) as expr_parser:
        for expr in expr_parser.read():
            if not np.isclose(expr.expr_vals.std(), 0):
                expr_db[expr.symbol] = expr

    logger.info("reading the regions")
    with open(args.region_file) as reg_handle:
        reg_queue = JoinableQueue()
        reg_count = 0
        reg_parser = RegionParser(reg_handle)
        for reg in reg_parser.parse():
            if reg.chrom == args.chrom and reg.gene_name in expr_db:
                reg_queue.put({"region": reg})
                reg_count += 1

    with open(args.out_file, "w") as out_handle,\
            open(args.bg_stat_file, "w") as bg_stat_handle:
        data_queue = JoinableQueue(32)
        out_queue = JoinableQueue(32)
        plotting_queue = JoinableQueue(32)

        # start the worker subprocesses
        if not args.random_mode:
            dworker = DataWorker(
                reg_queue, data_queue, expr_db, dhs_tree
            )
        else:
            dworker = RandomDataWorker(
                reg_queue, data_queue, expr_db, random_db,
                args.scramble_count
            )
        dworker.start()
        logger.info("started data worker")

        # calculating the p-values is highly parallelized
        for i in range(args.n_threads):
            cworker = CorrWorker(
                data_queue, out_queue, plotting_queue, probing_samples,
                random_samples, args.min_random_samples, args.min_tss_skip,
                args.min_bg_tss, args.ks_probe_threshold, args.ks_bg_threshold,
                args.ks_split
            )
            cworker.start()
            logger.info("started correlation worker %s", i)

        # create the output file writer process
        wworker = WriterWorker(
            out_handle, bg_stat_handle, out_queue, chrom_e, chrom_g,
            reg_count
        )
        wworker.start()

        # create the qq-plotting process
        plotting_worker = QQPlottingWorker(
            plotting_queue, args.qq_dir, enabled=args.qq_plotting
        )
        plotting_worker.start()

        reg_queue.join()
        data_queue.join()
        out_queue.join()
        plotting_queue.join()

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
            dhs = list(self.dhs_tree.query_overlapping(reg_iv))

            qans = {
                "region": reg,
                "probe": probe,
                "dhs": dhs,
            }
            self.data_queue.put(qans)
            self.reg_queue.task_done()


class RandomDataWorker(Process):

    def __init__(self, reg_queue, data_queue, expr_db, random_db,
                 random_draws):
        Process.__init__(self)
        self.daemon = True

        self.reg_queue = reg_queue
        self.data_queue = data_queue
        self.expr_db = expr_db
        self.random_db = random_db
        self.random_draws = random_draws

    def run(self):
        ran_db = self.random_db
        ran_draws = self.random_draws

        while True:
            qmsg = self.reg_queue.get()
            reg = qmsg["region"]

            probe = self.expr_db[reg.gene_name]

            for region, mock_probe, dhs_list in random.sample(ran_db, ran_draws):
                mock_probe = Expr(
                    mock_probe.chrom,
                    mock_probe.start_loc,
                    mock_probe.end_loc,
                    probe.symbol,
                    mock_probe.probeset_id,
                    probe.expr_vals
                )
                qans = {
                    "region": region,
                    "probe": mock_probe,
                    "dhs": dhs_list,
                }
                self.data_queue.put(qans)

            self.reg_queue.task_done()


class CorrWorker(Process):

    def __init__(self, data_queue, out_queue, plotting_queue, random_probes,
                 random_samples, min_n_samples, min_dist=1e7, min_n_tss=10,
                 ks_probe_thresh=0.2, ks_bg_thresh=0.1, ks_split=0.5):

        Process.__init__(self)
        self.daemon = True

        self.data_queue = data_queue
        self.out_queue = out_queue
        self.plotting_queue = plotting_queue
        self.random_probes = random_probes
        self.random_samples = random_samples
        self.min_n_samples = min_n_samples
        self.min_dist = min_dist
        self.min_n_tss = min_n_tss
        self.ks_probe_thresh = ks_probe_thresh
        self.ks_bg_thresh = ks_bg_thresh
        self.ks_split = ks_split

    def run(self):
        random_probes = self.random_probes
        random_samples = self.random_samples
        min_n_samples = self.min_n_samples
        min_dist = self.min_dist
        min_n_tss = self.min_n_tss
        ks_probe_thresh = self.ks_probe_thresh
        ks_bg_thresh = self.ks_bg_thresh
        ks_split = self.ks_split

        ks_pval_stack = []
        for i in range(20 * min_n_tss):
            ks_pval_stack.append(PvalWrapper(0, None))
        spare_wrapper = PvalWrapper(0, None)

        while True:
            qmsg = self.data_queue.get()
            skip = qmsg.get("skip", False)
            reg = qmsg["region"]
            if skip:
                self.out_queue.put({"gene_name": reg.gene_name, "skip": True})
                self.data_queue.task_done()
                continue
            logger.debug("start prediction of gene %r", reg.gene_name)
            # prepare auxiliary variables
            probe = qmsg["probe"]
            expr = probe.expr_vals
            dhs_list = qmsg["dhs"]
            reg_center = reg.start + (reg.end - reg.start) // 2
            reg_strand = reg.strand

            # calculate correlation distribution for predictions
            cor_list = []
            for dhs in dhs_list:
                cor = pearsonr(dhs.access_vals, expr)
                cor_list.append(cor)
            cors = np.array(cor_list)
            cors.sort()
            split_cor = np.percentile(cors, q=ks_split*100)
            split_cor_index = np.searchsorted(cors, split_cor)
            lower_cors = cors[:split_cor_index]

            # do a full probe screening
            for wrapper in ks_pval_stack:
                wrapper.pval = 0
                wrapper.gene = None

            for probe_id, chrom, tss, probe_sample in random_probes:

                # do not sample from regions that could contain true positives
                if chrom == reg.chrom and abs(tss - reg_center) < min_dist:
                    continue

                sample_size = len(probe_sample)
                probe_cors = np.empty(sample_size)
                for i, acc in enumerate(probe_sample):
                    probe_cors[i] = pearsonr(acc, expr)
                probe_cors.sort()

                probe_split_index = np.searchsorted(probe_cors, split_cor)
                lower_probe_cors = probe_cors[:probe_split_index]
                ks_stat, pval = stats.ks_2samp(lower_cors, lower_probe_cors)

                spare_wrapper.pval = pval
                spare_wrapper.gene = probe_id

                # push to a minimum heap which keeps the maximal observed
                # p-values
                spare_wrapper = heapq.heappushpop(ks_pval_stack, spare_wrapper)

            sampled = 0
            cur_samples = 0
            bg_model_cors = []

            for wrapper in sorted(ks_pval_stack, reverse=True):

                if wrapper.pval < ks_probe_thresh or sampled >= min_n_tss:
                    break

                random_sample = random_samples[wrapper.gene]
                bg_cors = np.empty(len(random_sample))
                for i, acc in enumerate(random_sample):
                    bg_cors[i] = pearsonr(acc, expr)
                bg_cors.sort()

                bg_split_index = np.searchsorted(bg_cors, split_cor)
                lower_bg_cors = bg_cors[:bg_split_index]
                ks_stat, pval = stats.ks_2samp(lower_cors, lower_bg_cors)

                if pval < ks_bg_thresh:
                    continue

                bg_model_cors.append(bg_cors)
                cur_samples += len(bg_cors)
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
                bg_model = np.concatenate(bg_model_cors)
                bg_model.sort()
                msg = "successfully collected %s samples from %s tss"
                logger.debug(msg, cur_samples, sampled)

            results = []
            for dhs in dhs_list:
                cor = pearsonr(dhs.access_vals, expr)
                left_pos = bisect.bisect_left(bg_model, cor)
                right_pos = bisect.bisect_right(bg_model, cor)
                insert_pos = np.round(left_pos + (right_pos - left_pos) / 2)
                pval_e = 1 - insert_pos / cur_samples
                res = (
                    (dhs.start, dhs.end, dhs.id),
                    (probe.symbol, reg_strand, reg_center),
                    (cor, pval_e),
                )
                results.append(res)

            # calculate statistics on the background model
            pred_count = len(cors)
            bg_stats = {}
            bg_stats["bg_mean"] = np.mean(bg_model)
            bg_stats["bg_median"] = np.median(bg_model)
            bg_stats["bg_std"] = np.std(bg_model)
            bg_stats["bg_n"] = cur_samples
            bg_stats["bg_tss"] = len(bg_model_cors)
            bg_stats["pred_mean"] = np.mean(cors)
            bg_stats["pred_median"] = np.median(cors)
            bg_stats["pred_std"] = np.std(cors)
            bg_stats["pred_n"] = pred_count

            pred_distr = cors[:split_cor_index]
            bg_distr_split_index = np.searchsorted(bg_model, split_cor)
            bg_distr = bg_model[:bg_distr_split_index]
            ks_stat, ks_pval = stats.ks_2samp(pred_distr, bg_distr)
            bg_stats["ks_stat"] = ks_stat
            bg_stats["ks_pval"] = ks_pval

            # estimate eta0
            positions = []
            median_cor = np.median(cors)
            for distr in bg_model_cors:
                distr.sort()
                pos = np.searchsorted(distr, median_cor)
                rel_pos = pos / len(distr)
                positions.append(rel_pos)
            eta0 = 0.5 / np.mean(positions)
            bg_stats["eta0"] = eta0

            qans = {
                "gene_name": reg.gene_name,
                "predictions": results,
                "bg_stats": bg_stats,
                "skip": False,
            }
            self.out_queue.put(qans)

            plotting_ans = {
                "gene_name": reg.gene_name,
                "bg_cors": bg_model_cors,
                "pred_cors": cors,
            }
            self.plotting_queue.put(plotting_ans)

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
            "DHS", "DHS_id", "gene", "gene_symb", "cor", "pval_e",
            sep="\t", file=fit_handle
        )
        print(
            "gene", "bg_median", "bg_mean", "bg_std", "bg_n", "bg_tss",
            "pred_median", "pred_mean", "pred_std", "pred_n",
            "ks_stat", "ks_pval", "eta0",
            sep="\t", file=bgstat_handle
        )

        while True:
            qmsg = self.out_queue.get()
            skip = qmsg.get("skip", False)
            gene = qmsg["gene_name"]

            n_finished += 1
            if not skip:
                result = qmsg["predictions"]
                for dhs, probe, statistics in result:
                    d_start, d_end, d_id = dhs
                    p_symb, reg_strand, reg_center = probe
                    cor, pval_e = statistics
                    print(
                        "{}|{}|{}".format(self.chrom_e, d_start, d_end),
                        d_id,
                        "{}|{}|{}|{}".format(p_symb, self.chrom_g,
                                             reg_strand, reg_center),
                        p_symb,
                        "{:.5f}".format(cor),
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
                    sts["ks_stat"], sts["ks_pval"], sts["eta0"],
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


class QQPlottingWorker(Process):

    def __init__(self, plotting_queue, output_dir, quantile_count=500,
                 enabled=True):
        Process.__init__(self)
        self.daemon = True

        self.plotting_queue = plotting_queue
        self.quantile_count = quantile_count
        self.output_dir = output_dir
        self.enabled = enabled

    def run(self):
        percentiles = np.linspace(0, 100, self.quantile_count)
        while True:
            qmsg = self.plotting_queue.get()
            if not self.enabled:
                self.plotting_queue.task_done()
                continue
            skip = qmsg.get("skip", False)
            gene_name = qmsg["gene_name"]
            if skip:
                self.plotting_queue.task_done()
                logger.debug("skipped plotting of gene %r", gene_name)
                continue

            bg_cors = qmsg["bg_cors"]
            pred_cors = qmsg["pred_cors"]
            bg_quant = np.percentile(bg_cors, q=percentiles)
            pred_quant = np.percentile(pred_cors, q=percentiles)

            min_quant = min(bg_quant[0], pred_quant[0])
            max_quant = max(bg_quant[-1], bg_quant[-1])

            fig, ax = plt.subplots(figsize=(8, 8))
            fig.suptitle(gene_name, fontsize=16)
            ax.set_xlabel("Background Correlations", fontsize=14)
            ax.set_ylabel("Prediction Correlations", fontsize=14)

            # draw reference line
            ax.plot([min_quant, max_quant], [min_quant, max_quant],
                    color="r", lw=2)
            # scatter qq data
            ax.scatter(bg_quant, pred_quant, color="b", alpha=0.3, s=30)

            outfile = os.path.join(self.output_dir, gene_name + ".png")
            try:
                fig.savefig(outfile)
                logger.debug("successfully created qq-plot for gene %r",
                             gene_name)
            except Exception as e:
                logger.warn(e)
                logger.debug("could not create qq-plot for gene %r", gene_name)
            finally:
                plt.close(fig)

            self.plotting_queue.task_done()


@total_ordering
class PvalWrapper:
    def __init__(self, pval, gene):
        self.pval = pval
        self.gene = gene

    def __eq__(self, other):
        return self.pval == other.pval

    def __lt__(self, other):
        return self.pval < other.pval


if __name__ == "__main__":
    main()
