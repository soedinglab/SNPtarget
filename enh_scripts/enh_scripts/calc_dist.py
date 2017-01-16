import argparse
import csv


def main():
    parser = argparse.ArgumentParser(
        prog="calc_dist",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("in_file", help="input file")
    parser.add_argument("out_file", help="output file")
    parser.add_argument("--gene_col_name", default="gene",
                        help="name of the gene column")
    parser.add_argument("--dhs_col_name", default="DHS",
                        help="name of the dhs column")
    parser.add_argument("--dist_col_name", default="distance",
                        help="name of the newly created distance column")
    args = parser.parse_args()

    dist_col = args.dist_col_name
    gene_col = args.gene_col_name
    dhs_col = args.dhs_col_name

    with open(args.in_file, "r") as in_handle,\
            open(args.out_file, "w") as out_handle:

        dialect = csv.Sniffer().sniff(in_handle.read(2056))
        in_handle.seek(0)
        reader = csv.DictReader(in_handle, dialect=dialect)
        header = reader.fieldnames[:]
        header.append(dist_col)
        writer = csv.DictWriter(out_handle, dialect=dialect, fieldnames=header)
        writer.writeheader()

        for row in reader:
            e_chrom, e_start, e_end = row[dhs_col].split("|")
            g_symb, g_chrom, g_strand, g_center = row[gene_col].split("|")
            e_start, e_end = int(e_start), int(e_end)
            g_center = int(g_center)

            e_center = e_start + (e_end - e_start) // 2
            eg_dist = -(g_center - e_center)

            # upstream-downstream sensitive distance
            if g_strand == "-":
                eg_dist *= -1

            row[dist_col] = eg_dist
            writer.writerow(row)


if __name__ == "__main__":
    main()
