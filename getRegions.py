import sys
sys.path.append("/data/robert.williamson/bin/")
import summary
import argparse
from collections import namedtuple
import random

def main():
    summaryReader = summary.Reader(open(args.summary))

    if args.sample not in summaryReader.Samples:
        sys.stderr.write("No sample called %s in the summary. Possible samples: %s\n Did you mispell something?\n" % (args.sample, summaryReader.Samples))
        sys.exit()

    regionGetter = getNextRegion()
    nextRegion = next(regionGetter)
    try:
        for record in summaryReader:
            if nextRegion.scaf != record.CHROM:
                continue

            if nextRegion.end < record.POS:
                outputRegion(nextRegion)
                nextRegion = next(regionGetter)

            if nextRegion.start <= record.POS and nextRegion.end >= record.POS:
                nextRegion.reference.append(record.REF)
                try:
                    if len(record.REF) != len(record.ALT):
                        sys.stderr.write("INDEL at %s %s. skipping this position.\n" % (record.CHROM, record.POS))
                        continue
                    elif record.Genotypes[args.sample] == "R":
                        nextRegion.new_hap1.append(record.REF)
                        nextRegion.new_hap2.append(record.REF)
                    elif record.Genotypes[args.sample] == "A":
                        nextRegion.new_hap1.append(record.ALT)
                        nextRegion.new_hap2.append(record.ALT)
                    elif record.Genotypes[args.sample] == "H":
                        if random.random() <= 0.5:
                            nextRegion.new_hap1.append(record.ALT)
                            nextRegion.new_hap2.append(record.REF)
                        else:
                            nextRegion.new_hap1.append(record.REF)
                            nextRegion.new_hap2.append(record.ALT)
                    else:
                        nextRegion.new_hap1.append("N")
                        nextRegion.new_hap2.append("N")
                except KeyError:
                    sys.stderr.write("Missing genotype for %s at %s %s. Adding Ns\n" % (args.sample, record.CHROM, record.POS))
                    nextRegion.new_hap1.append("N")
                    nextRegion.new_hap2.append("N")
            if nextRegion.end == record.POS:
                outputRegion(nextRegion)
                nextRegion = next(regionGetter)

    except NoMoreRegions:
        pass

class NoMoreRegions(Exception):
    pass
        
RegionData = namedtuple('RegionData', ['scaf','start','end','samp1','hap1','hap2','reference','new_hap1', 'new_hap2'])

def outputAllRegions(region):
    if args.haplotypes:
        return
    if len(region.new_hap1) != len(region.new_hap2):
        sys.stderr.write("WARNING: region at %s %s has 2 different haplotype lengths.\n" % (region.scaf, region.start))
    print(">%s_haplotype1 %s %s %s"%(region.samp1, region.scaf, region.start, region.end))
    print("%s"% ("".join(region.new_hap1)))
    print(">%s_haplotype2 %s %s %s"%(region.samp1, region.scaf, region.start, region.end))
    print("%s"%("".join(region.new_hap2)))

def outputRegion(region):
    if not args.haplotypes:
        outputAllRegions(region)
        return
    newfile = "%s.%s.%s_%s.fa" % (args.prefix, region.scaf, region.start, region.end)
    if len(region.hap1) != len(region.hap2) or len(region.hap1) != len(region.new_hap1):
        sys.stderr.write("Haplotypes for %s not the same length, skipping.\n" % (newfile))
        #sys.stderr.write("%s\n%s\n"% (region.hap1, "".join(region.new_hap1)))
        return

    openFile = open(newfile, "w")

    openFile.write(">reference %s %s %s\n" % (region.scaf, region.start, region.end))
    openFile.write("%s\n" % ("".join(region.reference)))
    openFile.write(">%s_haplotype1 %s %s %s\n" % (region.samp1, region.scaf, region.start, region.end))
    openFile.write("%s\n" % (region.hap1))
    openFile.write(">%s_haplotype2 %s %s %s\n" % (region.samp1, region.scaf, region.start, region.end))
    openFile.write("%s\n" % (region.hap2))
    openFile.write(">%s_haplotype1 %s %s %s\n" % (args.sample, region.scaf, region.start, region.end))
    openFile.write("%s\n" % ("".join(region.new_hap1)))
    openFile.write(">%s_haplotype2 %s %s %s\n" % (args.sample, region.scaf, region.start, region.end))
    openFile.write("%s\n" % ("".join(region.new_hap2)))
    openFile.close()
    sys.stderr.write("Output %s.\n" % (newfile))

def getNextRegion():
    region = ""
    samp = ""
    scaf = ""
    start = 0
    end = 0
    seq1 = ""
    seq2 = ""
    if args.haplotypes:
        haplo = open(args.haplotypes, 'r')
        for line in haplo:
            line = line.rstrip()
            if not line:                    
                continue
        
            if line.startswith(">"):
                if seq1 and seq2:
                    yield RegionData(scaf=scaf, start=start, end=end, samp1=samp, hap1=seq1, hap2=seq2, reference=[], new_hap1=[], new_hap2=[])
                    seq1 = ""
                if not seq1:
                    sline = line.split()
                    samp = sline[0][1:]
                    region = str(sline[1:])
                    scaf = sline[1]#region.split(":")[0]
                    start,end = int(sline[2]), int(sline[3])#map(int, region.split(":")[1].split("-"))
                    seq1 = ""
                    seq2 = ""
                elif not seq2:
                    continue
            else:
                if not seq1:
                    seq1 = line
                else:
                    seq2 = line
        yield RegionData(scaf=scaf, start=start, end=end, samp1=samp, hap1=seq1, hap2=seq2, reference=[], new_hap1=[], new_hap2=[])
    else: #reading from a regions file
        regions = open(args.regions, "r")
        for line in regions:
            scaf, start, end = line.split()
            yield RegionData(scaf=scaf, start=int(start), end=int(end), samp1=args.sample, hap1="", hap2="", reference=[], new_hap1=[], new_hap2=[])
    raise NoMoreRegions

args = None
def parseArgs():
    parser = argparse.ArgumentParser(description="takes a fasta of haplotypes or a list of regions, a summary, and a sample name. Adds the sequences from the given sample and teh reference for each haplotype in the original file (along with those haplotypes) to a new fasta.\nThe fasta containing the haplotypes MUST be sorted in the same way as the summary.")
    parser.add_argument("-a", "--haplotypes", default=None, type=str, help="the fasta of haplotypes, names need to contain the positions of each haplotype.")
    parser.add_argument("-r", "--regions", default=None, type=str, help="a file listing regions to make haplotypes for")
    parser.add_argument("summary", help="the summary file")
    parser.add_argument("sample", help="the sample name to add to the fasta")
    parser.add_argument("prefix", help="the prefix for the output fasta filenames")
    return parser.parse_args()
if __name__=='__main__':
    global args
    args = parseArgs()
    main()
