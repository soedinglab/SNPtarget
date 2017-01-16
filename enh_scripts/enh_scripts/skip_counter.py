import argparse
import re


def main():
    parser = argparse.ArgumentParser("skip_counter")
    parser.add_argument("log_file")
    args = parser.parse_args()

    regex = re.compile("(?P<state>finished|skipped) gene '(?P<gene>.*)'")
    print("gene", "skip", sep="\t")
    with open(args.log_file) as handle:
        for line in handle:
            match = regex.search(line)
            if match:
                gene = match.group("gene")
                skip = 1 if match.group("state") == "skipped" else 0
                print(gene, skip, sep="\t")


if __name__ == '__main__':
    main()
