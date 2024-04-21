import sys
from itertools import groupby
from collections import defaultdict


#vcfFile = open(sys.argv[1])
#outVCF= open(sys.argv[2],'w')

def joinGT(gt):
    if "1" in gt:
        return "1"
    else:
        return "0"

for l in sys.stdin:
    if l[:2] =="##":
        print(l.strip())
    elif l[:1] == "#":
        l=l.strip().split("\t")
        newHeader=l[:9]
        for i in range(9,len(l)-1,2):
            newHeader.append("%s-%s"%(l[i],l[i+1]))
        print("\t".join(newHeader))
    else:
        l=l.strip().split("\t")
        newLine=l[:9]
        for i in range(9,len(l)-1,2):
            newLine.append("%s|%s"%(joinGT(l[i]),joinGT(l[i+1])))
        print("\t".join(newLine))
            

