import sys

f2=open('%s/excluded.txt'%sys.argv[1],'w')

sampleswithhits=[]
with open('%s/included.txt'%sys.argv[1]) as f:
    for line in f:
        sample=line.strip()
        sampleswithhits.append(sample)
sampleswithhits=set(sampleswithhits)

originalsamples=[]
with open('%s'%sys.argv[2]) as f:
    for line in f:
        data=line.strip().split('\t')
        sample=data[0]
        originalsamples.append(sample)
originalsamples=set(originalsamples)

excludedsamples=sorted(list(originalsamples.difference(sampleswithhits)))

if len(excludedsamples)>0:
    for sample in excludedsamples:
        f2.write('%s\n'%sample)
f2.close()

#flag if there are no blast hits (terminate pipeline in distance.py script, if no hits)
if len(sampleswithhits)==0:
    print('no blast hits')
