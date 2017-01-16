import argparse
from ngsbiotools.ivtree import Interval, IVTree


def create_parser():
    parser = argparse.ArgumentParser('filter_promoters')
    parser.add_argument('access_file')
    parser.add_argument('trscr_file')
    parser.add_argument('chrom')
    parser.add_argument('--updown_bp', type=int, default=500)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    tss_tree = IVTree()
    with open(args.trscr_file) as trscr_file:
        line = trscr_file.readline()
        for line in trscr_file:
            toks = line.split()
            if toks[0] != args.chrom:
                continue
            tss = int(toks[4])
            tss_iv = Interval(tss - args.updown_bp, tss + args.updown_bp)
            tss_tree.insert(tss_iv)

    with open(args.access_file) as access_file:
        header = access_file.readline()
        print(header.strip())
        for line in access_file:
            snp_id, chrom, pos_str, *_ = line.split()
            snp_pos = int(pos_str)
            snp_iv = Interval(snp_pos, snp_pos)
            if tss_tree.has_overlap(snp_iv):
                continue
            print(line.strip())


if __name__ == '__main__':
    main()
