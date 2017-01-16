import argparse
import csv

import numpy as np


def main():
    parser = argparse.ArgumentParser(
            prog="posterior_pval",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("pred_file", help="file with predictions")
    parser.add_argument("shuffled_file", help="file with shuffled predictions")
    parser.add_argument("output_file", help="output file with p-values")
    parser.add_argument("--gene_col_name", default="gene_symb",
                        help="name of the gene column")
    parser.add_argument("--posterior_col_name", default="posterior",
                        help="name of the posterior column")
    parser.add_argument("--pvalue_col_name", default="post_pvalue",
                        help="name of the newly create p-value column")
    parser.add_argument("--mode", choices=["left", "right", "mid", "rand"],
                        default="rand", help="insertion strategy")
    parser.add_argument("--seed", type=int, default=None,
                        help="seed of the random generator")
    args = parser.parse_args()

    gene_col = args.gene_col_name
    post_col = args.posterior_col_name
    pval_col = args.pvalue_col_name

    with open(args.shuffled_file, "r") as shuff_handle:
        dialect = csv.Sniffer().sniff(shuff_handle.read(2024))
        shuff_handle.seek(0)
        reader = csv.DictReader(shuff_handle, dialect=dialect)

        unread_reader = RecordUnreadReader(reader)
        random_db = {}
        cur_gene = None
        cur_post = []
        for row in unread_reader:
            if row[gene_col] != cur_gene:
                if len(cur_post) > 0:
                    post_array = np.array(cur_post)
                    post_array.sort()
                    random_db[cur_gene] = post_array
                    cur_post.clear()
                cur_gene = row[gene_col]
                unread_reader.unread()
            else:
                posterior = float(row[post_col])
                cur_post.append(posterior)

        if len(cur_post) > 0:
            post_array = np.array(cur_post)
            post_array.sort()
            random_db[cur_gene] = post_array

    with open(args.pred_file, "r") as in_handle,\
            open(args.output_file, "w") as out_handle:

        dialect = csv.Sniffer().sniff(in_handle.read(2056))
        in_handle.seek(0)
        reader = csv.DictReader(in_handle, dialect=dialect)
        header = reader.fieldnames[:]
        header.append(pval_col)
        writer = csv.DictWriter(out_handle, dialect=dialect, fieldnames=header)
        writer.writeheader()

        for row in reader:
            gene_symb = row[gene_col]
            bg_post = random_db[gene_symb]
            posterior = float(row[post_col])

            left_insert = np.searchsorted(bg_post, posterior, side="left")
            right_insert = np.searchsorted(bg_post, posterior, side="right")

            if args.mode == "left":
                insert_pos = left_insert
            elif args.mode == "right":
                insert_pos = right_insert
            elif args.mode == "mid":
                insert_pos = left_insert + (right_insert - left_insert) // 2
            elif args.mode == "rand":
                insert_pos = np.random.randint(left_insert, right_insert + 1)

            pvalue = 1 - insert_pos / len(bg_post)
            row[pval_col] = pvalue
            writer.writerow(row)


class RecordUnreadReader():

    def __init__(self, reader):
        self._iter = iter(reader)
        self._unread_flag = False
        self._buffer = None

    def unread(self):
        self._unread_flag = True

    def __iter__(self):
        return self

    def __next__(self):
        if self._buffer is not None and self._unread_flag:
            self._unread_flag = False
            return self._buffer
        else:
            next_item = next(self._iter)
            self._unread_flag = False
            self._buffer = next_item
            return next_item


if __name__ == "__main__":
    main()
