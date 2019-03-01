import sys

outdir=sys.argv[1]


#using stdin (concatenated gff files), or a pre-existing file, get a dictionary of genome:gffdata
genomenamedict={}

gfffile=sys.stdin.read().split('\n') #concatenated gff file from concatenategff.sh script
for line in gfffile:
    line=line.strip()
    if line.startswith('##gff-version') or line.startswith('##sequence-region'):
        continue
    data=line.split('\t')
    genomename=data[0].split('|')[0]
    if genomename in genomenamedict:
        genomenamedict[genomename].append('%s\n'%'\t'.join(data))
    else:
        genomenamedict[genomename]=[]
        genomenamedict[genomename].append('%s\n'%'\t'.join(data))

            
#write gff dictionary to file
genomenames=sorted(genomenamedict.keys())
for genomename in genomenames:
    f2=open('%s/%s.gff'%(outdir,genomename),'w')
    values=genomenamedict[genomename]
    for indx,value in enumerate(values):
        if indx==0:
            f2.write('##gff-version 3\n')
        f2.write(value)
    f2.close()
