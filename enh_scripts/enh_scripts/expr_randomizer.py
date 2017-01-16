import argparse
from enh_scripts import HUMAN_CHROMS
from collections import defaultdict
import random
import numpy as np


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('expr_db')
    parser.add_argument('--shuffle', action='store_true')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    if not args.shuffle:
        sample_db = defaultdict(list)
        with open(args.expr_db) as expr_db:
            expr_db.readline()
            for line in expr_db:
                chrom, _, _, _, _, *expr_data = line.split()
                for sample_chrom in HUMAN_CHROMS:
                    if chrom != sample_chrom:
                        sample_db[sample_chrom].append(expr_data)

    with open(args.expr_db) as expr_db:
        header = expr_db.readline()
        print(header.strip())
        for line in expr_db:
            chrom, start, end, gene, probeid, *values = line.split()
            if args.shuffle:
                new_values = np.random.permutation(values)
            else:
                new_values = random.choice(sample_db[chrom])
            print(chrom, start, end, gene, probeid, *new_values, sep=' ')


if __name__ == '__main__':
    main()
