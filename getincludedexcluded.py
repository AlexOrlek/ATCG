import sys,os
from itertools import combinations

#get included/excluded samples + excluded comparisons

#to get excluded comparisons...
#for -p flag, use filepath1/2 tsv
#for -s flag, use filepathinfo.tsv col 1
#for -1 -2 flags, use filepathinfo1/2.tsv

outputpath=sys.argv[1]
blasttype=sys.argv[2]  #sequencepairs, pairwise, allvallpairwise


#get included samples and included comparisons
includedsamples=[]
includedcomparisons=set()
with open('%s/output/comparisonstats.tsv'%outputpath) as f:  #N.B genome1/genome2 in comparisonstats.tsv are sorted in alphabetical order
    for indx,line in enumerate(f):
        if indx==0:
            continue
        data=line.strip().split()
        sample1=data[0]
        sample2=data[1]
        includedsamples.extend([sample1,sample2])
        includedcomparisons.add((sample1,sample2))
includedsamples=set(includedsamples)

#get excluded samples
originalsamples=[]
with open('%s/filepathinfo.tsv'%outputpath) as f:
    for line in f:
        data=line.strip().split('\t')
        sample=data[0]
        originalsamples.append(sample)
originalsamples=set(originalsamples)

excludedsamples=sorted(list(originalsamples.difference(includedsamples)))

#write included/excluded samples to file
f2=open('%s/included.txt'%outputpath,'w')
for sample in sorted(list(includedsamples)):
    f2.write('%s\n'%sample)
f2.close()

f2=open('%s/excluded.txt'%outputpath,'w')
for sample in excludedsamples:
    f2.write('%s\n'%sample)
f2.close()

#get excluded comparisons

#for blasttype==sequencepairs, compare includedcomparisons with comparisons in sequencepairs file
f2=open('%s/excludedcomparisons.tsv'%outputpath,'w')
excludedcomparisonsamples=set()
if blasttype=='sequencepairs':
    originalcomparisons=set()
    #get samples1 and samples2
    samples1=[]
    samples2=[]
    with open('%s/filepathinfo1.tsv'%outputpath) as f:
        for line in f:
            sample1=line.strip().split('\t')[0]
            samples1.append(sample1)
    with open('%s/filepathinfo2.tsv'%outputpath) as f:
        for line in f:
            sample2=line.strip().split('\t')[0]
            samples2.append(sample2)
    #add sequencepair tuples to originalcomparisons
    for sample1,sample2 in zip(samples1,samples2): 
        samplepair=tuple(sorted([sample1,sample2]))
        originalcomparisons.add(samplepair)
    #get excluded comparisons
    excludedcomparisons=sorted(list(originalcomparisons.difference(includedcomparisons)))
    for comparison in excludedcomparisons:
        sample1=comparison[0]
        sample2=comparison[1]
        excludedcomparisonsamples.add(sample1)
        excludedcomparisonsamples.add(sample2)
        f2.write('%s\t%s\n'%(sample1,sample2))


if blasttype=='allvallpairwise':
    #get samples
    samples=[]
    with open('%s/filepathinfo.tsv'%outputpath) as f:
        for line in f:
            sample=line.strip().split('\t')[0]
            samples.append(sample)
    #get all v all comparisons among samples and add as tuples to originalcomparisons
    samples=sorted(samples)
    originalcomparisons=combinations(samples,2)
    originalcomparisons=set(list(originalcomparisons))
    #get excluded comparisons
    excludedcomparisons=sorted(list(originalcomparisons.difference(includedcomparisons)))
    for comparison in excludedcomparisons:
        sample1=comparison[0]
        sample2=comparison[1]
        excludedcomparisonsamples.add(sample1)
        excludedcomparisonsamples.add(sample2)
        f2.write('%s\t%s\n'%(sample1,sample2))

if blasttype=='pairwise':
    originalcomparisons=set()
    #get samples1 and samples2
    samples1=[]
    samples2=[]
    with open('%s/filepathinfo1.tsv'%outputpath) as f:
        for line in f:
            sample1=line.strip().split('\t')[0]
            samples1.append(sample1)
    with open('%s/filepathinfo2.tsv'%outputpath) as f:
        for line in f:
            sample2=line.strip().split('\t')[0]
            samples2.append(sample2)
    #get comparisons
    for qsample in samples1:
        for ssample in samples2:
            samplepair=tuple(sorted([qsample,ssample]))
            originalcomparisons.add(samplepair)
    #get excluded comparisons
    excludedcomparisons=sorted(list(originalcomparisons.difference(includedcomparisons)))
    for comparison in excludedcomparisons:
        sample1=comparison[0]
        sample2=comparison[1]
        excludedcomparisonsamples.add(sample1)
        excludedcomparisonsamples.add(sample2)
        f2.write('%s\t%s\n'%(sample1,sample2)) 

f2.close()


#flag if there are excluded comparisons involving included samples
if len(excludedcomparisonsamples.intersection(includedsamples))>0:
    print('missingcomparison')