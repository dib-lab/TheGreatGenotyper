import sys


v_count=0
for l in sys.stdin:
    if l[0] == "#":
        print(l.strip())
    else:
       l=l.strip().split("\t")
       l[2]=f"{l[2]}_{v_count}"
       v_count+=1
       print("\t".join(l))

