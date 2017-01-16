import argparse
import numpy as np
import os
import glob

from collections import namedtuple
from ngsbiotools.ivtree import Interval, IVTree
from ngsbiotools.parsers import RegionParser


def create_parser():
    parser = argparse.ArgumentParser(
        prog="list_enh_gene_pairs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("region_file", help="region file for correlations")
    parser.add_argument("access_db", help="path to the accessiblity database")
    parser.add_argument("expr_db", help="path to expression database")
    parser.add_argument("chrom")

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    chrom = args.chrom

    acc_fglob = os.path.join(args.access_db, chrom + ".*")
    acc_fpath, = glob.glob(acc_fglob)

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

    expr_fglob = os.path.join(args.expr_db, chrom + ".*")
    expr_fpath, = glob.glob(expr_fglob)

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

    print('DHS', 'gene', sep='\t')
    with open(args.region_file) as reg_handle:
        reg_parser = RegionParser(reg_handle)
        for reg in reg_parser.parse():
            if (reg.chrom != chrom or
                    reg.gene_name not in chrom_expr_db):
                continue
            iv = Interval(reg.start, reg.end)

            reg_center = reg.start + (reg.end - reg.start) // 2
            for dhs in sample_tree.query_all_overlaps(iv):
                print(
                    "{}|{}|{}".format(chrom, dhs.start, dhs.end),
                    "{}|{}|{}|{}".format(reg.gene_name, chrom,
                                         reg.strand, reg_center),
                    sep='\t'
                )


Expression = namedtuple("Expression", ["symbol", "expr_vals"])


if __name__ == "__main__":
    main()
