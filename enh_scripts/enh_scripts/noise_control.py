import argparse
import gzip
from itertools import chain
import os

import numpy as np


def adv_open(fname):
    name, ext = os.path.splitext(fname)
    if ext == ".gz":
        return gzip.open(fname, "rt")
    else:
        return open(fname, "r")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", "-m", choices=["log_smooth", "min_value"],
                        default="min_value")
    parser.add_argument("--type", "-t", choices=["dhs", "expr"],
                        default="dhs")
    parser.add_argument("--value", "-v", type=float, default=0.1)
    parser.add_argument("--min_maxvalue", type=float, default=0)
    parser.add_argument("file")
    parser.add_argument("--rounding_precision", type=int, default=4)
    parser.add_argument("--header", choices=["y", "n"], default="y")
    parser.add_argument("--standardize", action="store_true")
    parser.add_argument("--std_epsilon", type=float, default=0.001)
    args = parser.parse_args()

    value_start = (3 if args.type == "dhs" else 5)
    with adv_open(args.file) as f:
        if args.header == "y":
            header = f.readline()
            head_toks = header.split()
            print(*head_toks, sep="\t")
        for line in f:
            toks = line.split()
            vals = np.array(toks[value_start:], dtype=float)
            max_val = np.max(vals)
            if max_val < args.min_maxvalue:
                continue
            if args.mode == "log_smooth":
                vals = np.log2(2**vals + 2**args.value)
            else:
                vals[vals < args.value] = args.value
            if args.standardize:
                mean = np.mean(vals)
                std = np.std(vals)
                vals = (vals - mean) / (std + args.std_epsilon)
            vals = np.round(vals, 4)
            print(*chain(toks[:value_start], vals), sep="\t")

if __name__ == "__main__":
    main()
