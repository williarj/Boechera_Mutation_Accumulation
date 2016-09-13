import sys
sys.path.append("/data/robert.williamson/bin")
import summary
import argparse
from collections import defaultdict

def main():
    #open the summary
    reader = summary.Reader(open(args.summary, 'rb'))
    
    #setup the bariables
    end = 0
    window = defaultdict(lambda : defaultdict(int))
    window_calls = defaultdict(lambda :defaultdict(int))
    start = 0
    pChrom = None
    window_len = 0

    if args.window:
        print("CHROM POS SAMP1 SAMP2 DIV CALLED")

    for record in reader:
        samps_hets = 0
        if record.REF_NUM == len(reader.Samples) or record.REF_NUM == 0:
            #site is not actually a SNP
            continue
        
        for i in range(0, len(record.Genos)):
            for j in range(i+1, len(record.Genos)):
                s1 = reader.Samples[i]
                s2 = reader.Samples[j]
                g = (record.Genos[i], record.Genos[j])
                if "N" not in g:
                    window_calls[s1][s2] += 1
                if (("H", "H") == g or g[0] != g[1]) and "N" not in g:
                    if (g == ("R", "A")) or (g == ("A", "R")):      
                        window[s1][s2] += 1
                    else:
                        window[s1][s2] += 0.5 
        if not pChrom:
            pChrom = record.CHROM
         
        if not start:
            start = record.POS

        window_len += 1
    
        if pChrom != record.CHROM or window_len >= args.window:
            for s in window.keys():
                for s2 in window[s].keys():
                    print("%s %s %s %s %s %s" % (pChrom, (start+end)/2.0, s, s2, window[s][s2], window_calls[s][s2])) 


            window = defaultdict(lambda : defaultdict(int)) 
            window_calls = defaultdict(lambda :defaultdict(int))
            start = record.POS
            pChrom = record.CHROM
            window_len = 0

        end = record.POS
    
    for s in window.keys():
        for s2 in window[s].keys():
            print("%s %s %s %s %s" % (pChrom, (start+end)/2.0, s, s2, window[s][s2], window_calls[s][s2])) 

 

def parseArgs():
    parser = argparse.ArgumentParser(description="Takes a Summary file and counts the number of pairwise differences between each pair of samples..")
    parser.add_argument("summary", type=str, help="the summary file to get data from")
    parser.add_argument("-w", "--window", type=int, default = 20000, help="the size of the sliding windows of the counts are reported.")

    return parser.parse_args()

args = parseArgs()

if __name__=='__main__':
    main()
