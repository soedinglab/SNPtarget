import argparse
import os


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('def_file')
    parser.add_argument('db')
    parser.add_argument('thresh', type=float)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    with open(args.def_file) as def_file:
        for line in def_file:
            _, locus, lead_snp_id = line.split()

            loc, loc_id, _ = locus.split('.')
            snp_file = '.'.join((loc, loc_id, 'summary'))
            ld_file = '.'.join((loc, loc_id, 'LD', 'p1'))

            snp_list = []
            snp_path = os.path.join(args.db, snp_file)
            with open(snp_path) as snps:
                for line in snps:
                    snp_id, *_ = line.split()
                    snp_list.append(snp_id)
            lead_pos = snp_list.index(lead_snp_id)

            ld_path = os.path.join(args.db, ld_file)
            with open(ld_path) as ld_snps:
                for i, line in enumerate(ld_snps):
                    if i == lead_pos:
                        break
                for i, lead_cor in enumerate(line.split()):
                    if float(lead_cor) > args.thresh:
                        print(snp_list[i])


if __name__ == '__main__':
    main()
