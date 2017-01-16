from __future__ import print_function

from pyivtree import IVTree
from pyivtree import GenomicInterval as Interval
from collections import defaultdict
import argparse
import os

__version__ = 1.0


def main():
    parser = argparse.ArgumentParser('snp_finder')
    parser.add_argument('posterior_threshold', type=float)
    parser.add_argument('regions')
    parser.add_argument('--database', default='.')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()
    thresh = args.posterior_threshold
    search_ivs = defaultdict(list)
    with open(args.regions) as regfile:
        for line in regfile:
            chrom, start, end, seq_id, *_ = line.split()
            iv = Interval(int(start), int(end))
            iv.seq_id = seq_id
            search_ivs[chrom].append(iv)

    hits = []
    header = None
    for chrom, iv_list in search_ivs.items():
        db_file = os.path.join(args.database, chrom + '.post_pred')
        db_tree = IVTree()
        with open(db_file) as db_file:
            header = db_file.readline().split()
            for line in db_file:
                toks = line.split()
                chrom, start, end = toks[0].split('|')
                db_iv = Interval(int(start), int(end))
                db_iv.toks = toks
                db_tree.insert(db_iv)
        for iv in iv_list:
            for hit_iv in db_tree.query_all_overlaps(iv):
                posterior = float(hit_iv.toks[-1])
                if posterior < thresh:
                    continue
                hits.append(hit_iv.toks)

    def posterior_cmp(toks):
        return float(toks[-1])

    sorted_hits = sorted(hits, key=posterior_cmp, reverse=True)

    if header is not None:
        print(*header, sep='\t')
        for hit in sorted_hits:
            print(*hit, sep='\t')


if __name__ == '__main__':
    main()
