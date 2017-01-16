import argparse
from collections import defaultdict

from ngsbiotools.ivtree import IVTree, Interval


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('region_file')
    parser.add_argument('pred_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    reg_db = defaultdict(IVTree)
    with open(args.region_file) as reg_file:
        for line in reg_file:
            chrom, start_str, end_str, *_ = line.split()
            reg_iv = Interval(int(start_str), int(end_str))
            reg_db[chrom].insert(reg_iv)

    with open(args.pred_file) as pred_file:
        header = pred_file.readline()
        print(header.strip())
        for line in pred_file:
            toks = line.split()
            chrom, dhs_start, dhs_end = toks[0].split('|')
            dhs_iv = Interval(int(dhs_start), int(dhs_end))
            if reg_db[chrom].has_overlap(dhs_iv):
                print(line.strip())


if __name__ == '__main__':
    main()
