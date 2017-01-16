import argparse
import re


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('encode_chiapet_file')
    return parser


def main(chiapet_file):
    annot_regex = re.compile(
        '(chr[0-9XY]+):([0-9]+)..([0-9]+)-'  # left anchor
        '(chr[0-9XY]+):([0-9]+)..([0-9]+),'  # right anchor
        '([0-9]+)'                           # number of tags

    )
    with open(chiapet_file) as infile:
        header = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'tags']
        print(*header, sep='\t')

        for line in infile:
            toks = line.split()
            annot_col = toks[3]
            match = annot_regex.match(annot_col)
            assert match is not None
            hits = [match.group(i + 1) for i in range(7)]
            print(*hits, sep='\t')


if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    main(args.encode_chiapet_file)
