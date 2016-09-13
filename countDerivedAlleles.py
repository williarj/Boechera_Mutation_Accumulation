import sys
sys.path.append("/data/robert.williamson/bin")
import summary
import argparse

def main():
    # read in samples
    samps = []
    for line in open(args.samples):
        samps.append(line.rstrip())

    #open the summary
    reader = summary.Reader(open(args.summary, 'rb'))
    geno_codes = {v:k for k,v in reader.Genotypes.items()}
    
    #setup the bariables
    snps = 0
    number_het = 0
    unique_snps = 0
    unique_hets = 0
    window = 0
    start = 0
    end = 0
    pChrom = None

    if args.window:
        print("CHROM POS Unique_SNPs Hets Unique_Hets")

    for record in reader:
        samps_AF_alt = 0
        samps_AF_ref = 0
        samps_hets = 0
        if record.REF_NUM == len(reader.Samples) or record.REF_NUM == 0:
            #site is not actually a SNP
            continue
        try:
            for s in samps:
                if record.Genotypes[s] == geno_codes["homozygote alternate"]:
                    samps_AF_alt += 2
                elif record.Genotypes[s] == geno_codes["homozygote reference"]:
                    samps_AF_ref += 2
                elif record.Genotypes[s] == geno_codes["heterozygote"]:
                    samps_AF_alt += 1
                    samps_AF_ref += 1
                    samps_hets += 1
        except KeyError:
            sys.stderr.write("Missing sample %s from like %s %s.\n" % (s, record.CHROM, record.POS))
            continue
        
        if not pChrom:
            pChrom = record.CHROM
         
        if not start:
            start = record.POS

        if args.window and pChrom != record.CHROM:
            unique_snps = 0
            number_het = 0
            unique_hets = 0
            window = 0
            start = record.POS
            pChrom = record.CHROM

        snps += 1
        if samps_AF_alt == samps_hets:
            #all samples are heterozygous
            number_het += 1
        if samps_AF_alt == record.ALT_NUM or samps_AF_ref == record.REF_NUM:
            #either all ref or all alt alleles are in our subsample
            unique_snps += 1
            if samps_AF_alt == samps_hets:
            #and everyone is a heterozygote
                unique_hets += 1
        window += 1
        end = record.POS
        if args.window and window == args.window:
            #scaf midpoint unique hets unique_hets
            print("%s %s %s %s %s" % (pChrom, (start+end)*0.5, unique_snps, number_het, unique_hets))
            unique_snps = 0
            number_het = 0
            unique_hets = 0
            window = 0
            start = 0

    if not args.window:
        print("KIND COUNT PERCENT\nSNPs %s %s\nUnique_SNPs %s %s\nFixed_hets %s %s\nUnique_fixed_hets %s %s" \
            % (snps, snps/float(snps), unique_snps, unique_snps/float(snps), number_het, number_het/float(snps), unique_hets, unique_hets/float(snps)))
    

def parseArgs():
    parser = argparse.ArgumentParser(description="Takes a Summary file with **only** SNPs and a list of samples. Counts the number of SNPs total, and the number of SNPs that have a unique allele in the list of samples and the number of fixed heterozygotes in those samples.")
    parser.add_argument("summary", type=str, help="the summary file to get data from")
    parser.add_argument("samples", type=str, help="the file listing the samples to count hets in.")
    parser.add_argument("-w", "--window", type=int, default = 0, help="if this option is used then sliding windows of the counts are reported.")

    return parser.parse_args()

args = parseArgs()

if __name__=='__main__':
    main()
