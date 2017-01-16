import argparse
import glob
import os

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def create_parser():
    parser = argparse.ArgumentParser("plot_pval_score")
    parser.add_argument("input_dir")
    parser.add_argument("--ext", default=".pred")
    parser.add_argument("output_dir")
    parser.add_argument("--normed_pval", action="store_true")
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    fglob = os.path.join(args.input_dir, "*" + args.ext)
    for fname in glob.glob(fglob):
        basename = os.path.basename(fname)
        df = pd.read_table(fname)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        ax1.hist(df.pval_e, bins=20, normed=args.normed_pval)
        ax2.hist(df.score, bins=40, log=True)
        out_fname = os.path.join(args.output_dir, basename + ".png")
        fig.savefig(out_fname)
        plt.close()

if __name__ == '__main__':
    main()
