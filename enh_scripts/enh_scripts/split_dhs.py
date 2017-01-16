import argparse
from collections import defaultdict
import gzip
from itertools import chain
import os

from ngsbiotools.parsers import TranscriptParser
from pyivtree import GenomicInterval as Interval
from pyivtree import IVTree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dhs_file")
    parser.add_argument("--tss_upstr", "-u", type=int, default=500)
    parser.add_argument("--tss_downstr", "-d", type=int, default=500)
    parser.add_argument("--trscr_file")
    parser.add_argument("--output_dir", "-o", default=".")
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        raise ValueError("The output directory has to be a directory.")
    if args.tss_upstr < 0:
        raise ValueError("TSS distance has to be a positive integer.")

    tss_trees = defaultdict(IVTree)
    if args.trscr_file:
        with open(args.trscr_file, "r") as trscr_file:
            trscr_parser = TranscriptParser(trscr_file)
            for trscr in trscr_parser.parse():
                chrom = trscr.chrom
                strand = trscr.strand
                tss = trscr.tss
                gene_name = trscr.gene_name
                if strand == "+":
                    tss_iv = Interval(tss - args.tss_upstr,
                                      tss + args.tss_downstr)
                else:
                    tss_iv = Interval(tss - args.tss_downstr,
                                      tss + args.tss_upstr)
                tss_iv.gene_name = gene_name
                tss_trees[chrom].insert(tss_iv)

    chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    try:
        chrom_files = {}
        for chrom in chroms:
            filename = chrom + ".dhs.gz"
            path = os.path.join(args.output_dir, filename)
            chrom_files[chrom] = gzip.open(path, "wt")
        if args.trscr_file:
            ppath = os.path.join(args.output_dir, "promotor.dhs.gz")
            prom_file = gzip.open(ppath, "wt")

        with gzip.open(args.dhs_file, "rt") as dhs_file:
            header = dhs_file.readline()
            head_toks = header.split()
            head_toks.insert(3, "dhs_id")

            # print headers
            for cfile in chrom_files.values():
                print(*head_toks, sep="\t", file=cfile)
            if args.trscr_file:
                print(*chain(head_toks, ["gene_name"]), sep="\t",
                      file=prom_file)

            dhs_id = 1
            dhs_id_pat = "dhs_{id:08d}"
            for line in dhs_file:
                tok = line.split()
                tok.insert(3, dhs_id_pat.format(id=dhs_id))
                chrom = tok[0]

                # check if the dhs overlaps a tss
                if args.trscr_file:
                    start = int(tok[1])
                    end = int(tok[2])
                    dhs_iv = Interval(start, end)

                    overlap = False
                    for tss in tss_trees[chrom].query_all_overlaps(dhs_iv):
                        print(*chain(tok, [tss.gene_name]), sep="\t",
                              file=prom_file)
                        overlap = True
                    if not overlap:
                        print(*tok, sep="\t", file=chrom_files[chrom])
                else:
                    print(*tok, sep="\t", file=chrom_files[chrom])
                dhs_id += 1
    finally:
        for cfile in chrom_files.values():
            cfile.close()
        if args.trscr_file:
            prom_file.close()

if __name__ == '__main__':
    main()
