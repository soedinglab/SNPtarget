import argparse
import numpy as np


def create_parser():
    parser = argparse.ArgumentParser("expr_error_estim")
    parser.add_argument("expr_file")
    return parser


def error_fun(x, a, b, c):
        return a * np.exp(-b * x * x) + c


def main():
    parser = create_parser()
    args = parser.parse_args()

    # the fitting happens in expr_err_estim.ipynb
    params = 0.662, 0.086, 0.042

    with open(args.expr_file) as expr_handle:
        header = expr_handle.readline().split()
        print(*header, sep=" ")

        for line in expr_handle:
            chrom, start, end, gene, ps_id, *vals = line.split()
            expr_arr = np.array(vals, dtype=float)
            err = np.sqrt(error_fun(expr_arr, *params))
            print(chrom, start, end, gene, ps_id, *err, sep=" ")


if __name__ == '__main__':
    main()
