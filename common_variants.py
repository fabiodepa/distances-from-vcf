#!/usr/bin/env python3
# common_variants.py
# AUTHOR: Fabio De Pascale
# LAST REVISED: 18/10/2018
#
# Department of Biology, Genomics Unit
# University of Padua

import sys, numpy, re, argparse

parser=argparse.ArgumentParser(description="")
parser.add_argument("-f", "--file", dest="infile", required=True, help="multisample vcf file")
parser.add_argument("-s", action='store_true', help="Design a symmetric matrix")
args=parser.parse_args()
infile=args.infile

args=parser.parse_args()
infile=args.infile

# the function CompareLineLoop scroll through the elements of a list in search of
### 1 where founded it searches other 1 in the list and add 1 to the position in
### the matrix that correspond to the comparison of that two samples.
def CompareLineLoop(line, n, matrix):
    i=0
    while i < (n-1):
        if line[i]==1:
            y=i+1
            while y < n:
                if line[i]==line[y]:
                    matrix[y,i]+=1
                    if args.s:
                        matrix[i,y]+=1
                y+=1
        i+=1
    return

# the fucntion CleanSampleInfotakes the genotype info fields and for each sample,
### each item, verifies the presence of that variant in that sample by looking
### at the genotype: "./." indicate that the given variant is absent in that
### sample; presence or absence are coded by 0=absence and 1=presence. A list of
### 0 and 1 is thus the return object of this function
def CleanSampleInfo(line):
    y=0
    array=[0]*len(line)
    for i in line:
        i=i.split(':')
        i=i[0]
        if i=="./.":
            i=0
        else:
            i=1
        array[y]=i
        y+=1
    return array

def main():
    f=open(infile, 'r')
    for line in f:
        if re.search('#CHROM', line):
            line=line.rstrip()
            line=line.split('\t')
            infofield=line.index('INFO')
            firstsample=infofield+2
            samplenames=line[(firstsample):]
            samplenum=len(samplenames)
            # create matrix
            matrix=numpy.zeros((samplenum,samplenum), dtype=numpy.int)
        elif line[:2]=='##':
            pass
        else:
            line=line.rstrip()
            line=line.split('\t')
            samplefields=line[(firstsample):]
            samplevars=CleanSampleInfo(samplefields)
            CompareLineLoop(samplevars, samplenum, matrix)

    y=0
    print("" + '\t' + '\t'.join(samplenames))
    for i in matrix:
        print(samplenames[y] + '\t' + '\t'.join(map(str, i)))
        y+=1
    return

if __name__ == '__main__':
    main()