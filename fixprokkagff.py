import sys,re, os

#N.B don't need to read stdin first -can read in other files / do other stuff
#assume all inputs are gff3 format


seqnamesfile=sys.argv[1] #a file with original genome names in column 1
outdir=sys.argv[2]
inputtype=sys.argv[3]
if inputtype=='file':
    gfffile=sys.argv[4] #a multigff file


#get original sequence names 
seqnames=[]
with open(seqnamesfile) as f:
    for line in f:
        seqname=line.strip().split('\t')[0]
        seqnames.append(seqname)

seqnames=sorted(list(set(seqnames)))

#from seqnames get prokkaseqnames (if seqname contains '|', substitute '_' and append to prokkaseqnames); get genome names in corresponding order
prokkaseqnames=[]
for seqname in seqnames:
    if '|' in seqname:
        prokkaseqname=re.sub(r"\|","_",seqname)
        prokkaseqnames.append(prokkaseqname)

#using stdin (concatenated gff files), or a pre-existing file, get a dictionary of genome:gffdata
genomenamedict={}
if inputtype=='file':
    with open(gfffile) as f:
        for line in f:
            line=line.strip()
            if line=='##FASTA':
                break
            if line.startswith('##gff-version') or line.startswith('##sequence-region'):
                continue
            data=line.split('\t')
            #replacementindx=[indx for indx,j in enumerate(prokkaseqnames) if j in data[0]]
            replacementindx=[indx for indx,j in enumerate(prokkaseqnames) if j==data[0]]  #data[0] is gff genome name
            if len(replacementindx)==1: #if a genome has been altered by prokka ("|" substituted for "_") ...
                data[0]=seqnames[replacementindx[0]]
            genomename=data[0].split('|')[0]
            if genomename in genomenamedict:
                genomenamedict[genomename].append('%s\n'%'\t'.join(data))
            else:
                genomenamedict[genomename]=[]
                genomenamedict[genomename].append('%s\n'%'\t'.join(data))
else:
    gfffile=sys.stdin.read().split('\n') #concatenated gff file from concatenategff.sh script
    for line in gfffile:
        line=line.strip()
        if line.startswith('##gff-version') or line.startswith('##sequence-region'):
            continue
        data=line.split('\t')
        replacementindx=[indx for indx,j in enumerate(prokkaseqnames) if j==data[0]]  #data[0] is gff genome name
        if len(replacementindx)==1: #if a genome has been altered by prokka ("|" substituted for "_") ...
            data[0]=seqnames[replacementindx[0]]
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

    

#OLD CODE

            #data[0]=re.sub(prokkaseqnames[replacementindx[0]],seqnames[replacementindx[0]],data[0])    #this is from when I was trying to accommodate the sequence region lines at the top so regular expression substitution was required - but just skip lines and replace prokka name with corresponding seqname
