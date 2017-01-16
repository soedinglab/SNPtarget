import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("expr_file")
    parser.add_argument("--output_dir", "-o", default=".")
    args = parser.parse_args()

    if not os.path.isdir(args.output_dir):
        raise ValueError("Output has to be a directory.")

    chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
    try:
        chrom_files = {}
        for chrom in chroms:
            filename = chrom + ".expr"
            path = os.path.join(args.output_dir, filename)
            chrom_files[chrom] = open(path, "w")

        with open(args.expr_file, "rt") as expr_file:
            header = expr_file.readline()
            toks = header.split()

            # print headers
            for cfile in chrom_files.values():
                print(*toks, sep="\t", file=cfile)

            for line in expr_file:
                tok = line.split()
                chrom = tok[0]
                print(*tok, sep="\t", file=chrom_files[chrom])

    finally:
        for cfile in chrom_files.values():
            cfile.close()

if __name__ == "__main__":
    main()
