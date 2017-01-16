import argparse
from collections import defaultdict

from ngsbiotools.ivtree import IVTree, Interval
from ngsbiotools.parsers import TabFileParser


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('annotation_file')
    parser.add_argument('prediction_file')
    parser.add_argument('--snp_db')
    return parser


def run(annotation_file, prediction_file, snp_db_file):

    if snp_db_file:
        snp_db = set()
        with open(snp_db_file) as snp_handle:
            for snp_id in snp_handle:
                snp_db.add(snp_id.strip())

    ia_db = defaultdict(IVTree)
    with open(annotation_file) as annot:
        annot.readline()
        for line in annot:
            chrom1, start1, end1, chrom2, start2, end2, *_ = line.split()
            # ignore interchromosomal interactions
            if chrom1 != chrom2:
                continue
            anchor_left = Interval(int(start1), int(end1))
            anchor_right = Interval(int(start2), int(end2))
            anchor_left.partner = anchor_right
            ia_db[chrom1].insert(anchor_left)

    with open(prediction_file) as pred_file:
        pred_parser = TabFileParser(pred_file)
        header = pred_parser.col_names
        header.append('ia_annot')
        print(*header, sep='\t')

        for rec in pred_parser.parse():

            if snp_db_file:
                if rec.DHS_id not in snp_db:
                    continue

            dhs_chrom, dhs_start_str, dhs_end_str = rec.DHS.split('|')
            _, gene_chrom, _, gene_tss_str = rec.gene.split('|')
            assert dhs_chrom == gene_chrom

            dhs_start = int(dhs_start_str)
            dhs_end = int(dhs_end_str)
            dhs_iv = Interval(dhs_start, dhs_end)

            gene_tss = int(gene_tss_str)
            gene_iv = Interval(gene_tss, gene_tss)

            if dhs_start < gene_tss:
                left_anchor = dhs_iv
                right_anchor = gene_iv
            else:
                left_anchor = gene_iv
                right_anchor = dhs_iv

            overlap = False
            for ovl in ia_db[dhs_chrom].query_all_overlaps(left_anchor):
                if ovl.partner.overlaps(right_anchor):
                    overlap = True
                    break
            all_toks = pred_parser.lastline.split('\t')
            all_toks.append(int(overlap))
            print(*all_toks, sep='\t')


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args.annotation_file, args.prediction_file, args.snp_db)


if __name__ == '__main__':
    main()
