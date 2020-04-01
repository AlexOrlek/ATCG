import sys

outputpath=sys.argv[1]
blastdbs=sys.argv[2]

includedsamples=[]
with open('%s/output/distancestats.tsv'%outputpath) as f:
    for indx,line in enumerate(f):
        if indx==0:
            continue
        data=line.strip().split()
        sample1=data[0]
        sample2=data[1]
        includedsamples.extend([sample1,sample2])
includedsamples=set(includedsamples)

originalsamples=[]
with open('%s'%blastdbs) as f:
    for line in f:
        data=line.strip().split('\t')
        sample=data[0]
        originalsamples.append(sample)
originalsamples=set(originalsamples)

excludedsamples=sorted(list(originalsamples.difference(includedsamples)))

f2=open('%s/included.txt'%outputpath,'w')
for sample in sorted(list(includedsamples)):
    f2.write('%s\n'%sample)
f2.close()

f2=open('%s/excluded.txt'%outputpath,'w')
for sample in excludedsamples:
    f2.write('%s\n'%sample)
f2.close()
