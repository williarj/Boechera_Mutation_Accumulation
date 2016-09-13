import argparse, sys
sys.path.append("/data/robert.williamson/bin/")
import summary
from collections import defaultdict

def main():
    reader = summary.Reader(open(args.summary, 'rb'))

    site_counts = defaultdict(int)
    snp_counts = defaultdict(int)
    het_counts = defaultdict(int)
    alt_counts = defaultdict(int)
    for record in reader:
        for ind in record.Genos:
            if ind == "N":
                continue
            het = False
            if ind == "H":
                het = True
            alt = False
            if ind == "A":
                alt = True
    
            for type in record.Types:
                if record.ALT != ".":
                    snp_counts[type] += 1
                site_counts[type] += 1
                if het:
                    het_counts[type] += 1
                if alt:
                    alt_counts[type] += 1

    outputhets(site_counts, snp_counts, het_counts, alt_counts)

def parseArgs():
    parser = argparse.ArgumentParser(description="gets the total heterozygosity for all site types in a summary.")
    parser.add_argument("summary", type=str, help="the summary to process")

    return parser.parse_args()

def outputhets(sites, snp, hets, alt):
    print("Type Sites SNPs Hets Alt %")
    for t in sites.keys():
        print("%s %s %s %s %s %s" % (t, sites[t], snp[t], hets[t], alt[t], 1/float(sites[t])/hets[t]))

if __name__=="__main__":
    args=parseArgs()
    main()
