import argparse

import numpy as np
import pandas as pd


def create_parser():
    parser = argparse.ArgumentParser('subsample_dist_bf')
    parser.add_argument('all_ia')
    parser.add_argument('pred_ia')
    parser.add_argument('out_file')
    parser.add_argument('--n_bins', type=int, default=40)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    all_df = pd.read_table(args.all_ia)
    min_dist = all_df.distance.min()
    max_dist = all_df.distance.max()

    bins = np.linspace(min_dist, max_dist, args.n_bins + 1)
    bin_labels = np.arange(args.n_bins)

    all_cut = pd.cut(all_df.distance, bins, labels=bin_labels)
    bin_probs = all_cut.value_counts() / len(all_cut)

    del all_df

    pred_df = pd.read_table(args.pred_ia)
    n_pred = len(pred_df)
    pred_cut = pd.cut(pred_df.distance, bins, labels=bin_labels)

    indices = []
    for cat in pred_cut.cat.categories:
        n_sample = n_pred * bin_probs[cat]
        ind = pred_cut[pred_cut == cat].index
        indices.append(np.random.choice(ind, size=int(n_sample)))
    pred_list = [pred_df.loc[ind] for ind in indices]
    joined_df = pd.concat(pred_list)
    joined_df.to_csv(args.out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
