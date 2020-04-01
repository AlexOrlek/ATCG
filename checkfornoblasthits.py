import sys

sampleswithhits=[]
with open('%s/included.txt'%sys.argv[1]) as f:
    for line in f:
        sample=line.strip()
        sampleswithhits.append(sample)
sampleswithhits=set(sampleswithhits)


#flag if there are no blast hits (terminate pipeline in distance.py script, if no hits)
if len(sampleswithhits)==0:
    print('no blast hits')
