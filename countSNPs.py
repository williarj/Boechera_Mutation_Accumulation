import sys, argparse
from collections import defaultdict, namedtuple
from countTrees import parseTree, sortTree

args = None
def main():
    #for each pair of tree/fasta
    #read in the tree, parse it
    #print it pretty to the new file
    #read in the sequences and types
    #count the SNPs
    #close
    homoplasious = 0
    for line in open(args.inputs):
        sline = line.split()
        ID_index = sline.index("ID")
        id = sline[ID_index+1]
        fasta = sline[ID_index-1]

        output = open(args.count_file_name+id, 'w')
        
        tree_filename = args.tree_name.replace("ID", id)
        try:
            tree = parseTree(open(tree_filename).readline())
        except IOError:
            sys.stderr.write("No tree file for %s tried opening %s.\n" % (fasta, tree_filename))
            continue
        output.write(str(tree).replace(" ","").replace("'","")+"\n")
        output.write(prettyTree(sortTree(tree))+"\n")

        seqs = {}
        head = ""
        seq = ""
        for line in open(fasta):
            if line.startswith(">"):
                if seq:
                    seqs[head] = seq
                seq = ""
                head = line[1:].split()[0]
            else:
                seq += line.rstrip()
        if seq:
            seqs[head] = seq


        if "types" not in seqs.keys():
            sys.stderr.write("No types in file %s.\n" % fasta)
            continue
        else:
            if "," in seqs["types"]:
                seqs["types"] = seqs["types"].split(",")
        types = set(seqs["types"])
        snps = countSNPs(tree, seqs, {"homoplasious":[]})
        #sys.stderr.write("%s\n" % snps) 
        sites = []
        for s in snps.values():
            sites += s
        sites = list(set(sites))

        output.write("SNP sites %s %s\n" % (len(sites), sites))
        output.write("Homoplasious sites %s %s\n" % (len(snps["homoplasious"]), snps["homoplasious"]))
        output.write("BRANCH TYPE COUNT TOTAL_NUM\n")
        for branch, s in sorted(snps.items()):
            if branch == "homoplasious" or branch == "reference":
                continue
            counts = {t:0 for t in types}
            for pos in s:
                if pos not in snps["homoplasious"]:
                    counts[seqs["types"][pos]] += 1
            for t in types:
                output.write("%s %s %s %s\n" % (branch, t, counts[t], seqs["types"].count(t)))

        output.close()

def parseArgs():
    parser = argparse.ArgumentParser(description="takes a list of fastas with site types on tope, and trees and then counts SNPs on each branch of each tree")
    parser.add_argument("inputs", help="A file that lists the fastas and ID numbers of the trees from RAxML associated with them.")
    parser.add_argument("tree_name", help="The format and path of the tree files, the ID number should be replaced with 'ID' (which also obviously cant be anywhere else in the filename")
    parser.add_argument("count_file_name", default="SNP_counts_", help="the prefix of the output filenames, the ID number is appended to this")
    parser.add_argument("-t", "--test", action = "store_true", help="tests the function in this script")
    parser.add_argument("-v", "--verbose", action="store_true", help="turns on verbose mode")
    parser.add_argument("-a", "--asexual", type = str, help="the sample name of the asexual in these trees, if provided it replaces the names with 'ASEX'")
    parser.add_argument("-s", "--sexual", type = str, help="the sample name of the sexual in these trees, if provided it replaces the names with 'SEX'")

    return parser.parse_args()

def prettyTree(tree):
    sorted = sortTree(tree)    
    s = str(sorted).replace("[","(").replace("]",")").replace("'", "") 
    s = s.replace("_haplotype1","").replace("_haplotype2","")
    if args.asexual:
        s = s.replace(args.asexual, "ASEX")
    if args.sexual:
        s = s.replace(args.sexual, "SEX")
    return s

"""
seqs needs to have a sequence for every leaf and one for 'referece' and 'types'
one 
"""
def countSNPs(tree, seqs, acc):
    '''
    ~~~~ Test 1
    >>> seqs = {"reference":"AAAAAAAAA", "types":"111111111", "A":"TAAAAAAAA", "B":"TTAAAATAA", "C":"TTTAAAAAA", "D":"TTTTAAAAA"}
    >>> tree = ["reference", ["A", ["B", ["C", "D"]]]]
    >>> list(sorted(countSNPs(tree, seqs, {"homoplasious":[]}).items()))
    [('A', []), ('B', [6]), ('C', []), ('D', [3]), ('[A,[B,[C,D]]]', [0]), ('[B,[C,D]]', [1]), ('[C,D]', [2]), ('homoplasious', []), ('reference', [])]

    ~~~~~ Test 2 homoplasious, seq A site #2 has a back mutation
    >>> tree = ["reference", [["A","B"], ["C","D"]]]
    >>> list(sorted(countSNPs(tree, seqs, {"homoplasious":[]}).items()))
    [('A', []), ('B', [1, 6]), ('C', []), ('D', [3]), ('[A,B]', []), ('[C,D]', [2]), ('[[A,B],[C,D]]', [0]), ('homoplasious', [1]), ('reference', [])]

    ~~~~~ Test 3 homoplasious, seq A site # back mutation again
    >>> tree = ["reference", ["D", ["C", ["B", "A"]]]]
    >>> list(sorted(countSNPs(tree, seqs, {"homoplasious":[]}).items()))
    [('A', []), ('B', [1, 6]), ('C', [2]), ('D', [3]), ('[B,A]', []), ('[C,[B,A]]', []), ('[D,[C,[B,A]]]', [0]), ('homoplasious', [1, 2]), ('reference', [])]

    ~~~~~ Test 4 
    >>> tree = ["reference", [["A","B"], ["C","D"]]]
    >>> seqs["A"] = "TTAATAAAA"
    >>> seqs["B"] = "TTAATATAA"
    >>> list(sorted(countSNPs(tree, seqs, {"homoplasious":[]}).items()))
    [('A', []), ('B', [6]), ('C', []), ('D', [3]), ('[A,B]', [4]), ('[C,D]', [2]), ('[[A,B],[C,D]]', [0, 1]), ('homoplasious', []), ('reference', [])]

    '''

   # sys.stderr.write(str(seqs)+"\n")
    left = tree[0]
    right = tree[1]
    if type(left) is list:
        left = countSNPs(left, seqs, acc)
    else:
        left = seqs[left]

    if type(right) is list:
        right = countSNPs(right, seqs, acc)
    else:
        right = seqs[right]

    newSeq = ""
    left_mutations = []
    right_mutations = []
    homoplasious_mutations = []

    for i, bases in enumerate(zip(left, right)):
        l, r = bases
        if l != r and "N" not in bases:
            if l == seqs["reference"][i]:
                right_mutations.append(i)
            else:
                left_mutations.append(i)
            newSeq += "N"
        elif "N" in bases and seqs["reference"][i] not in bases:
            homoplasious_mutations.append(i)
            newSeq += "N"
        else:
            newSeq += l

    acc[str(tree[0]).replace("'","").replace(" ","")] = left_mutations
    acc[str(tree[1]).replace("'","").replace(" ","")] = right_mutations
    acc["homoplasious"] = list(set(acc["homoplasious"]+homoplasious_mutations))

    if 'reference' == tree[0] or 'reference' == tree[1]:
        return acc
    return newSeq

'''
returns strings of the counts for the branch above a given node (or leaf):
NODE TYPE COUNT
S1 1 10
S1 2 4
[S2, S3] 1 0
[S2, S3] 2 2
'''
def printCounts(counts, node):
    out = []
    for t in counts.keys():
        out.append("%s %s %s" % (str(node).replace("'","").replace(" ",""), t, counts[t]))
    return out

class HomoplasyException(Exception):
    pass

if __name__ == "__main__":
    global args
    args = parseArgs()
    if args.test:
        import doctest
        doctest.testmod(verbose=args.verbose)
    else:    
        main()
