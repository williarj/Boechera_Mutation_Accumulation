import sys
sys.path.append("/data/robert.williamson/bin")
import argparse
import summary
from collections import defaultdict

def main():
    exclude = defaultdict(lambda : defaultdict(int))
    if args.exclude:
        for line in open(args.exclude):
            sline = line.split()
            #esclude[CHROM][POS] = 1
            exclude[sline[0]][int(sline[1])] = 1

    include = defaultdict(lambda : defaultdict(int))
    if args.include:
        for line in open(args.include):
            sline = line.split()
            include[sline[0]][int(sline[1])] = 1

    regions = defaultdict(list)
    if args.regions:
        for line in open(args.regions):
            sline = line.split()
            #regions[CHROM].append((START, END))
            regions[sline[0]].append(map(int, (sline[1], sline[2])))

    reader = summary.Reader(open(args.summary,'rb'))
    #w = []# (scaf, pos, homo, het, called)
    w = defaultdict(list)
    print("SAMP SCAF MIDPOINT HOMO HET CALLED")
    for record in reader:
        if w and record.CHROM != w["all"][0][0]:
            outputAll(w)
            w = defaultdict(list)
        elif w and len(w["all"]) >= args.window:
            outputAll(w)
            w = defaultdict(list)

        if args.type and args.type not in record.Types:
            #print(record.Types)
            continue

        if args.exclude and exclude[record.CHROM].has_key(record.POS):
            #in exlude list
            continue

        if args.include and not include[record.CHROM].has_key(record.POS):
            #not in include list
            continue

        while args.regions and regions[record.CHROM] and regions[record.CHROM][0][1] < record.POS:
            regions[record.CHROM].pop(0)

        if args.regions and recions[record.CHROM] and record.POS >= regions[record.CHROM][0][0] and record.POS <= regions[record.CHROM][0][1]:
            #site is in an excluded region
            continue

        homo = 0
        het = 0
        called = 0
        for s, g in record.Genotypes.items():
            if g != "N":
                called += 1
            
            if g == "A":
                homo += 1
                if args.samples:
                    w[s].append((record.CHROM, record.POS, 1, 0, 1))
            elif g == "H":
                het += 1
                if args.samples:
                    w[s].append((record.CHROM, record.POS, 0, 1, 1))
            elif g == "R" and args.samples:
                w[s].append((record.CHROM, record.POS, 0, 0, 1))
        site = (record.CHROM, record.POS, homo, het, called)
        w["all"].append(site)
    outputAll(w)
    return

def outputAll(window):
    for k in window.keys():
        if args.samples and k == "all":
            continue
        output(window[k], k)

def output(window, samp):
    scaf = window[0][0]
    midpoint = (window[0][1]+window[-1][1])/2
    homo = 0
    het = 0
    called = 0
    for s, p, h, e, c in window:
        homo += h
        het += e
        called += c
    #homo = sum([h for s, p, h, e, c in window])
    #het = sum([e for s, p, h, e, c in window])
    #called = sum([c for s, p, h, e, c in window])
    print("%s %s %s %s %s %s" % (samp, scaf, midpoint, homo, het, called))

def parseArgs():
    parser = argparse.ArgumentParser(description="takes in a summary and in windows counts the number of hets the number of homo ref and hoo alt and the number of called sites")
    parser.add_argument("summary", help="The summary to read in.")
    parser.add_argument("-w", "--window", default=10000, type=int, help="The window size")
    parser.add_argument("-t", "--type", default=False, help="If this option is used only sites with this type are counted.")
    parser.add_argument("-s", "--samples", action="store_true", help="If this flag is used separate counts are kept for each sample.")
    parser.add_argument("-e", "--exclude", type=str, help="This option can be used to provide a list of sites that should be excluded from the analysis (e.g. fixed derived sites)")
    parser.add_argument("-r", "--regions", type=str, help="This option can be used to exclude whole regions from the analysis.")
    parser.add_argument("-i", "--include", type=str, help="This option can be used to specify a list of sites to include, all other sites will be excluded.")

    return parser.parse_args()

if __name__=="__main__":
    args = parseArgs()
    sys.stderr.write("%s\n"%args)
    main()
