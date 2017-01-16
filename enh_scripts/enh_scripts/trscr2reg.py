import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gc_trscr_file", help="GENCODE transcript file")
    parser.add_argument("out_file", help="output file")
    help_text = "+/- range from TSS"
    parser.add_argument("--range", type=int, default=100000, help=help_text)
    args = parser.parse_args()

    gc_tr = pd.read_table(args.gc_trscr_file)

    # filter levels
    tr = gc_tr[gc_tr.level < 3]

    grp = tr.groupby(["chrom", "strand", "gene_name"])

    def choose_tss(df, int_range):
        all_tss = df.tss.copy()
        all_tss.sort()
        med_tss = all_tss.iloc[len(all_tss)//2]
        return pd.DataFrame({
            "chrom": [df.chrom.iloc[0]],
            "gene_name": [df.gene_name.iloc[0]],
            "strand": [df.strand.iloc[0]],
            "start": [med_tss - int_range],
            "end": [med_tss + int_range]
        }, columns=["chrom", "start", "end", "gene_name", "strand"])

    result = grp.apply(choose_tss, int_range=args.range)
    result.reset_index(drop=True, inplace=True)

    result.to_csv(args.out_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
