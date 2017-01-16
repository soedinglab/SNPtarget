import argparse
import numpy as np


def create_parser():
    parser = argparse.ArgumentParser("snp_access_mask")
    parser.add_argument("--min_access", default=0.1, type=float)
    parser.add_argument("--remove_non_informative", action="store_true")
    parser.add_argument("infile")
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    min_access = args.min_access
    sparse = args.remove_non_informative

    with open(args.infile) as handle:
        header = handle.readline().split()
        print(*header, sep="\t")
        for line in handle:
            id_str, chrom, pos, *vals = line.split()
            access_arr = np.array(vals, dtype=float)
            access_arr[access_arr < min_access] = min_access

            if sparse:
                sd = access_arr.std()
                if np.allclose(sd, 0):
                    continue

            print(id_str, chrom, pos, *access_arr, sep="\t")



if __name__ == '__main__':
    main()
