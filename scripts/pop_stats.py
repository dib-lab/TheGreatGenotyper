import sys
from itertools import groupby
from pysam import VariantFile



vcfFile =VariantFile(sys.argv[1])

variant_list=set()
if len(sys.argv) > 2:
    for l in open(sys.argv[2]):
        variant_list.add(l.strip())
print(list(variant_list)[:10])        
#minQual=int(sys.argv[2])
#out=VariantFile(sys.argv[3],'w',header=vcfFile.header)


for rec in vcfFile.fetch():
    stats={}
    if len(variant_list) != 0 and rec.info["ID"][0] not in variant_list:
        continue
    for sample in rec.samples.keys():
        gt=rec.samples[sample]["GT"]
        if gt not in stats:
            stats[gt]={}
        gq=str(rec.samples[sample]["GQ"])    
        if gq not in stats[gt]:
            stats[gt][gq]=0
        stats[gt][gq]+=1
    result=""
    for k,v in stats.items():
        result+=str(k[0])+"/"+str(k[1])+": "
        freqs=list(v.items())
        freqs=sorted(freqs,key=lambda x: -x[1] )
        for k2,v2 in freqs:
            result+=k2+":"+str(v2)+", "
        result+="\n"
    print(rec.info["ID"][0])
    print(result)
        
    print("=============")
