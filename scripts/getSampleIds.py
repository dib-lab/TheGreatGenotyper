import sys
sampleNames=[  "HG00731_30_TheGreatGenotyper",  "HG00731_20_TheGreatGenotyper",  "HG00731_10_TheGreatGenotyper",  "HGDP00232",  "HG00731_5_TheGreatGenotyper",  "HGDP00003",  "HGDP00899",  "HGDP00102",  "HGDP00949",  "HGDP01303",  "HGDP00100",  "HGDP00110",  "HGDP01298",  "HGDP00860",  "HGDP01300",  "HGDP00444",  "HGDP00104",  "HGDP00346",  "HGDP00115",  "NA17386",  "HGDP00258",  "HGDP00195",  "Tuba19",  "Y4",  "HGDP00108",  "HG02790",  "K1",  "I3",  "HGDP00311",  "HGDP00433"]

index=0
result=[]
for l in open(sys.argv[1]):
    l=l.strip().split("\t")
   # l[0]=l[0].split("_")[1]
      
    if l[0] in sampleNames:
        result.append(str(index))
        print(l)
    index+=1

print("{ " + ", ".join(result) + "}")
