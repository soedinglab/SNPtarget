from __future__ import print_function

import argparse
import os
import glob
import sys

__version__ = 1.0


def main():
    parser = argparse.ArgumentParser('snp_finder')
    parser.add_argument('posterior_threshold', type=float)
    parser.add_argument('--snp_ids', nargs='+', default=[])
    parser.add_argument('--genes', nargs='+', default=[])
    parser.add_argument('--database', default='.')
    parser.add_argument('--chroms', nargs='+', default=[])
    parser.add_argument('--max_hits', type=int, default=10000)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()

    search_files = []
    if len(args.chroms) > 0:
        for chrom in args.chroms:
            db_file = os.path.join(args.database, chrom + '.post_pred')
            search_files.extend(glob.glob(db_file))
    else:
        glob_pattern = os.path.join(args.database, '*.post_pred')
        search_files = glob.glob(glob_pattern)

    if len(search_files) == 0:
        msg = (
            'There are no files in the database that match the given criteria. '
            'Please check if the --chroms argument is valid and that --database '
            'points to the directory containing the *.post_pred files'
        )
        print(msg, file=sys.stdout)
        sys.exit(1)

    snp_ids = set()
    for snp_id in args.snp_ids:
        snp_ids.add(snp_id)

    genes = set()
    for gene in args.genes:
        genes.add(gene)

    thresh = args.posterior_threshold
    max_hits = args.max_hits

    hits = []
    for pred_file in search_files:
        with open(pred_file) as infile:
            header = infile.readline().split()
            for line in infile:

                if len(hits) >= max_hits:
                    break

                toks = line.split()
                snp_id = toks[1]
                gene = toks[3]
                posterior = float(toks[-1])

                if len(snp_ids) > 0:
                    if snp_id not in snp_ids:
                        continue
                if len(genes) > 0:
                    if gene not in genes:
                        continue
                if posterior < thresh:
                    continue
                hits.append(toks)

    def posterior_cmp(toks):
        return float(toks[-1])

    sorted_hits = sorted(hits, key=posterior_cmp, reverse=True)

    print(*header, sep='\t')
    for hit in sorted_hits:
        print(*hit, sep='\t')


if __name__ == '__main__':
    main()
