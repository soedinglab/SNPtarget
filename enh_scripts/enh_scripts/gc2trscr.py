import argparse
import gzip
import re
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gc_file", help="the GENCODE annotation file")
    args = parser.parse_args()

    header = [
        "chrom", "start", "end", "strand", "tss",
        "trscr_id", "gene_name", "gene_id", "level"
    ]

    print(*header, sep="\t")

    gid_pat = re.compile("(?P<gid>ENSGR*[0-9]+)(\.[0-9]+)*")
    with gzip.open(args.gc_file, "rt") as gc_file:
        for line in gc_file:
            line = line.strip()
            # ignore header lines
            if line.startswith("#"):
                continue
            tok = line.split("\t")

            # only use transcripts
            if tok[2] != "transcript":
                continue

            tr_chrom = tok[0]
            tr_start = tok[3]
            tr_end = tok[4]
            tr_strand = tok[6]
            tr_tss = tr_start if tr_strand == "+" else tr_end

            info_toks = tok[8].strip(";").split(";")
            info = {}
            for entry in info_toks:
                entry = entry.strip()
                try:
                    key, value = entry.split(" ")
                except ValueError:
                    print("Cannot split \"{}\"".format(entry), file=sys.stderr)
                value = value.strip('"')
                info[key] = value

            gid_str = info["gene_id"]
            gid_mat = gid_pat.match(gid_str)
            if not gid_mat:
                print("Cannot extract gene_id from '{}'".format(gid_str),
                      file=sys.stderr)
                continue
            gene_id = gid_mat.group("gid")
            print(
                tr_chrom, tr_start, tr_end, tr_strand,
                tr_tss, info["transcript_id"], info["gene_name"],
                gene_id, info["level"],
                sep="\t"
            )

if __name__ == "__main__":
    main()
