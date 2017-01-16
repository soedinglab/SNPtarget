import argparse
from xml.etree.ElementTree import iterparse
import sys


def create_parser():
    parser = argparse.ArgumentParser("dbsnp2loc")
    return parser


def main():
    parser = create_parser()
    parser.parse_args()

    print("snp_id", "chrom", "location", sep="\t")
    xmlns = "{http://www.ncbi.nlm.nih.gov/SNP/docsum}"
    rs_tag = "%sRs" % xmlns
    comp_tag = "%sAssembly/%sComponent" % (xmlns, xmlns)
    loc_tag = "%sMapLoc" % xmlns
    for event, rs in iterparse(sys.stdin, events=["end"]):
        if rs.tag == rs_tag:
            rs_id = rs.get("rsId")
            component = rs.find(comp_tag)
            chrom = component.get("chromosome")
            start_str = component.get("start")
            maploc = component.find(loc_tag)
            pos_str = maploc.get("asnFrom")
            if start_str is None or pos_str is None:
                print("cannot find mapping for snp %s" % rs_id, file=sys.stderr)
                continue
            ref_pos = int(start_str) + int(pos_str) + 1
            print(rs_id, "chr%s" % chrom, ref_pos, sep="\t")

            rs.clear()

if __name__ == "__main__":
    main()
