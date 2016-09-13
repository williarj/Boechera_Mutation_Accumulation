A repository of programs used to analyze apomictic and sexual boechera genomes.

Most all of these scripts have help info, if run with a -h option will give details on inputs needed.

vcfSummarizer.py - takes a VCF file and a filter python script and creates a filtered and annotated SNP matrix
filterBoechera.py - the filter file for use with vcfSummarizer filters SNPs based on genotype quality and proximity to called indels
summary.py - a parser used to process snap matricies
getRegions.py - used to create the 'pseudohaplotypes' takes a list of regions and a SNP matrix and produces haplotypes
getHeterozygosity.py - takes a SNP matrix file and outputs windows of counts of called and heterozygous sites
countTrees.py - counts the topologies of all the trees given in a list of trees
countPairwiseDiversity.py - calculates pairwise differences between all samples in a SNP matrix
countDerivedAlleles.py - counts, in windows, the number of derived alleles in a subsample of a SNP matrix.
countSNPs.py - given a list of trees and their associated fasts this script counts the number of derived mutations on each branch of the trees.
