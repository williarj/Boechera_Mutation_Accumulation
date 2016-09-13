import sys, argparse
from collections import defaultdict

def main():
    trees = defaultdict(int)
    trees_seq = defaultdict(int)
    seqs = False

    lengths = {}

    if args.files:
        names = open(args.files)
    else:
        names = args.input

    for f in names:
        #assumes the files are only one line long
        sf = f.split()
        fname = sf[0]
        treeFile = open(args.path+fname.replace("\n",""))
        treeLine = treeFile.readline()
        treeFile.close()

        tree = sortTree(parseTree(treeLine))
        ptree, names = prettyTree(tree)
        trees[ptree] += 1
        if args.lengths:
            tLen = getLengths(treeLine)
            tLen = cleanDictNames(tLen, names)
            if ptree not in lengths.keys():
                lengths[ptree] = tLen
            else:
                for branch in tLen.keys():
                    lengths[ptree][branch] += tLen[branch]
        if len(sf) > 1:
            #print(sf[1])
            trees_seq[ptree] += int(sf[1])
            seqs = True

    for tree in trees.keys():
        if not args.lengths and seqs == False:
            print("%s\t%s" % (tree, trees[tree]))
        elif not args.lengths and seqs == True:
            print("%s\t%s\t%s" % (tree, trees[tree], trees_seq[tree]))
        elif seqs == True:
            for branch in lengths[tree].keys():
                lengths[tree][branch] = lengths[tree][branch] / trees[tree]
            tLengths = str(lengths[tree])[1:-1]
            print("%s\t%s\t%s\t%s" % (tree, trees[tree], tLengths, trees_seq[tree]))
        else:
            for branch in lengths[tree].keys():
                lengths[tree][branch] = lengths[tree][branch] / trees[tree]
            tLengths = str(lengths[tree])[1:-1]
            print("%s\t%s\t%s" % (tree, trees[tree], tLengths))
        

def cleanDictNames(d, names):
    newd = {}
    for k,v in d.items():
        newd[names[k]] = v
    return newd

def parseArgs():
    parser = argparse.ArgumentParser(description="This reads in tree files from RAxML and outputs the counts of every tree topology")
    parser.add_argument("input", nargs="*")
    parser.add_argument("-f", "--files", type = str, help="a file name where the list of tree files can be found, if used any tree files passed in on the command line are ignored")
    parser.add_argument("-a", "--asexual", type = str, help="the sample name of the asexual in these trees, if provided it replaces the names with 'ASEX'")
    parser.add_argument("-s", "--sexual", type = str, help="the sample name of the sexual in these trees, if provided it replaces the names with 'SEX'")
    parser.add_argument("-p", "--path", type = str, default="", help="for use with -f, appends this path to the front of every file name in the list of files")
    parser.add_argument("-l", "--lengths", action="store_true", help="if used the average branch length is also output")

    return parser.parse_args()

args = None

#this makes sure that trees rotated about a node count as the same tree
def sortTree(tree):
    for i, node in enumerate(tree):
        if type(node) is list:
            tree[i] = sorted(sortTree(node))
    tree.sort()
    return tree

def prettyTree(tree):
    sorted = sortTree(tree)    
    s = str(sorted).replace("[","(").replace("]",")").replace("'", "")
    s = s.replace("_haplotype1_haplotype2","_haplotype2").replace("_haplotype1_haplotype1","_haplotype1")
    key = getSimpleNames(s)

    for k,v in key.items():
        s = s.replace(k, v)

#    if args.asexual:
#        s = s.replace(args.asexual, "ASEX")
#    if args.sexual:
#        s = s.replace(args.sexual, "SEX")

    return s, key

def getSimpleNames(s):
    k = {}
    if (s.index(args.asexual+"_haplotype1") < s.index(args.asexual+"_haplotype2")):
        k[args.asexual+"_haplotype1"] = "ASEX1"
        k[args.asexual+"_haplotype2"] = "ASEX2"
    else:
        k[args.asexual+"_haplotype1"] = "ASEX2"       
        k[args.asexual+"_haplotype2"] = "ASEX1"

    if (s.index(args.sexual+"_haplotype1") < s.index(args.sexual+"_haplotype2")):
        k[args.sexual+"_haplotype1"] = "SEX1"       
        k[args.sexual+"_haplotype2"] = "SEX2"       
    else:
        k[args.sexual+"_haplotype1"] = "SEX2"
        k[args.sexual+"_haplotype2"] = "SEX1"         

    return k

def getLengths(s):
    lengths = {}
    lengths.update(getLengthsSub(s, "_haplotype1:"))
    lengths.update(getLengthsSub(s, "_haplotype2:"))
    return lengths

def getLengthsSub(s, hap):
    lengths = {}
    h1 = s.split(hap)           
    for i in range(0, len(h1)-1):
        if h1[i].startswith("0"):
            continue
        name = getName(h1[i]).replace("_haplotype1","")
        nxt = h1[i+1].split(",")
        if ")" in nxt[0]:
            nxt = h1[i+1].split(")")
            h1[i+1] = ")".join(nxt[1:])
        else:
            h1[i+1] = ",".join(nxt[1:])
        lengths[name+hap[:-1]] = float(nxt[0])
    return lengths

 
def getName(s):
    s = s.split("(")[-1]
    s = s.split(",")[-1]
    return s

def parseTree(s):
    node = []
    items = ""
    i = 0
    while i < len(s):
        #print(s[i], node)
        n = s[i]
        if n == ")":
            sitems = items.split(",")
            for taxon in sitems:
                #print(taxon)
                names = taxon.split(":")
                if names[0]:
                    node.append(names[0])
                    
            return sorted(node), s[i+1:]
        elif n == "(":
            res, rest = parseTree(s[i+1:])
            node.append(res)
            s = rest
            i = 0
            continue
        else:
            items += n
        i += 1
    return sorted(node[0])

if __name__ == "__main__":
    global args
    args = parseArgs()
    main()
