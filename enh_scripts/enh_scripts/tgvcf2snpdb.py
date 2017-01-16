import argparse
from . import open_or_stdin


def create_parser():
    parser = argparse.ArgumentParser("tgvcf2snpdb")
    parser.add_argument("vcf_file")
    parser.add_argument("--minimum_maf", type=float)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    min_maf = args.minimum_maf
    with open_or_stdin(args.vcf_file) as vcf_handle:
        print("snp_id", "chrom", "location", sep="\t")
        for line in vcf_handle:
            if line.startswith("#"):
                continue
            chrom, loc, ids, *tokens = line.split()

            annot = {}
            for info_str in tokens[4].split(";"):
                try:
                    key, value = info_str.split("=")
                    annot[key] = value
                except:
                    continue

            if min_maf is None:
                for snp_id in ids.split(";"):
                    print(snp_id, "chr%s" % chrom, loc, sep="\t")
            else:
                for af_str in annot["AF"].split(","):
                    af = float(af_str)
                    if af >= min_maf:
                        for snp_id in ids.split(";"):
                            print(snp_id, "chr%s" % chrom, loc, sep="\t")
                        continue


if __name__ == "__main__":
    main()
