import argparse
import numpy as np


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('snp_access_file')
    parser.add_argument('cell_line')
    parser.add_argument('access_thresh', type=float)
    return parser


def run(snp_access_file, cell_line, acc_thresh):
    with open(snp_access_file) as acc:
        header = acc.readline().split()
        cell_lines = header[3:]
        cell_index = cell_lines.index(cell_line)

        for line in acc:
            toks = line.split()
            acc_vals = np.array(toks[3:], dtype=float)

            if acc_vals[cell_index] > acc_thresh:
                print(toks[0])


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args.snp_access_file, args.cell_line, args.access_thresh)


if __name__ == '__main__':
    main()
