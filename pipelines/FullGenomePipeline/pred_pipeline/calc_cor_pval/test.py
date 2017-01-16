from ngsbiotools.parsers import RegionParser

with open('regions.reg') as reg:
    parser = RegionParser(reg)
    for rec in parser.parse():
        print(rec)
        break
