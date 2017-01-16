import argparse
from collections import defaultdict
import numpy as np
import os
import glob


def create_parser():
    parser = argparse.ArgumentParser(
        prog="select_tss",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("trscr_file")
    parser.add_argument("tss_access_db")
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    gene_db = defaultdict(int)
    gene_tss = {}
    acc_fglob = os.path.join(args.tss_access_db, "*")
    for access_path in glob.glob(acc_fglob):
        with open(access_path) as acc_handle:
            acc_handle.readline()
            for line in acc_handle:
                gene, chrom, tss, *values = line.split()
                access = np.array(values, dtype=float).sum()
                if gene_db[gene] <= access:
                    gene_db[gene] = access
                    gene_tss[gene] = tss

    with open(args.trscr_file) as trscr_handle:
        header = trscr_handle.readline()
        print(header.strip())
        for line in trscr_handle:
            toks = line.split()
            tss = toks[4]
            gene_name = toks[6]

            try:
                best_tss = gene_tss[gene_name]
            except KeyError:
                continue

            if best_tss == tss:
                print(line.strip())


if __name__ == "__main__":
    main()
