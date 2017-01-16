import argparse


def create_parser():
    parser = argparse.ArgumentParser("filter_by_snps")
    parser.add_argument("pred_file")
    parser.add_argument("snp_list")
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    snp_db = set()
    with open(args.snp_list) as db_handle:
        for line in db_handle:
            snp_id = line.strip()
            snp_db.add(snp_id)

        with open(args.pred_file) as pred_handle:
            header = pred_handle.readline().split()
            print(*header, sep='\t')
            for line in pred_handle:
                toks = line.split()
                if toks[1] in snp_db:
                    print(*toks, sep='\t')


if __name__ == '__main__':
    main()
