import sys,os,pickle,re,gzip
from Bio import SeqIO

pairsfile=sys.argv[1]
outputpath=sys.argv[2]
bidirectionalblast=str(sys.argv[3])

filepaths1=set()  #used for makeblastdb
filepaths2=set()
filepathpairs=set() #array of sorted filepath pairs (to check for directional duplicates)
comparisonsdict={}  #used for runblast


with open(pairsfile) as f:
    for line in f:
        data=line.strip().split('\t')
        assert len(data)==2,'Error: --sequencepairs file does not have 2 columns of data'
        filepath1=data[0]
        filepath2=data[1]
        #check both are valid filepaths
        if not os.path.isfile(filepath1):
            sys.exit('Error: %s is not a valid filepath'%filepath1)
        if not os.path.isfile(filepath2):
            sys.exit('Error: %s is not a valid filepath'%filepath2)
        #filepath1=os.path.abspath(filepath1)
        #filepath2=os.path.abspath(filepath2)
        #add paths to arrays of unique query/subject filepaths (for makeblastdb) and filepath sorted pairs (to check for duplicate comparisons)
        filepaths1.add(filepath1)
        filepaths2.add(filepath2)
        filepathpair=tuple(sorted([filepath1,filepath2]))
        if filepathpair in filepathpairs:
            sys.exit('Error: the sequence pair %s - %s is given more than once in --sequencepairs file'%(filepath1,filepath2))
        filepathpairs.add(filepathpair)
        #make dictionary of query-subject comparisons
        #if filepath2 not in comparisonsdict:
        #    comparisonsdict[filepath2]=[filepath1]
        #else:
        #    comparisonsdict[filepath2].append(filepath1)


#save filepaths.tsv with sample \t filepath; also write includedsubjects to file
#N.B can then use makeblastdb_dirinput.sh

f3=open('%s/allsubjects.txt'%outputpath,'w')
f2=open('%s/filepathinfo2.tsv'%outputpath,'w')
for filepath in filepaths2:
    if filepath.endswith('.gz'):
        gunzipfilepath=re.sub(r'\.gz$','',filepath)
        extension=os.path.splitext(gunzipfilepath)[1]
        sample=os.path.splitext(os.path.basename(gunzipfilepath))[0]
    else:
        extension=os.path.splitext(filepath)[1]
        sample=os.path.splitext(os.path.basename(filepath))[0]
    assert extension in {'.fa','.fasta','.fna'},'Error: unrecognised file extension: %s fasta files must have extension .fa .fasta or .fna'%extension
    blastdbpath='%s/blastdbs2/%s_db'%(outputpath,sample)
    f2.write('%s\t%s\t%s\n'%(sample,filepath,blastdbpath))
    f3.write('%s\n'%sample)
f2.close()


f2=open('%s/filepathinfo1.tsv'%outputpath,'w')
for filepath in filepaths1:
    if filepath.endswith('.gz'):
        gunzipfilepath=re.sub(r'\.gz$','',filepath)
        extension=os.path.splitext(gunzipfilepath)[1]
        sample=os.path.splitext(os.path.basename(gunzipfilepath))[0]
    else:
        extension=os.path.splitext(filepath)[1]
        sample=os.path.splitext(os.path.basename(filepath))[0]
    assert extension in {'.fa','.fasta','.fna'},'Error: unrecognised file extension: %s fasta files must have extension .fa .fasta or .fna'%extension
    blastdbpath='%s/blastdbs1/%s_db'%(outputpath,sample)
    f2.write('%s\t%s\t%s\n'%(sample,filepath,blastdbpath))
    if bidirectionalblast=='True':
        f3.write('%s\n'%sample)
f2.close()

f3.close()

#pickle comparison dict
#f=open('%s/comparisonsdict.pickle'%outputpath,'wb')
#pickle.dump(comparisonsdict,f)
#f.close()

#get seqlengths and filepathinfo
f2=open('%s/seqlengths.tsv'%outputpath,'w')
f3=open('%s/filepathinfo.tsv'%outputpath,'w')
fastafilepaths=sorted(list(filepaths1|filepaths2))
for fastafilepath in fastafilepaths:
    #open fastafilepath as "infile" file handle
    if fastafilepath.endswith('.gz'):
        infile=gzip.open(fastafilepath,'rt')
    else:
        infile=open(fastafilepath,'r')
    #get sample name (=filename without extension)
    if fastafilepath.endswith('.gz'):
        gunzipfilepath=re.sub(r'\.gz$','',fastafilepath)
        sample=os.path.splitext(os.path.basename(gunzipfilepath))[0]
    else:
        sample=os.path.splitext(os.path.basename(fastafilepath))[0]
    blastdbpath='%s/blastdbs1/%s_db'%(outputpath,sample)
    f3.write('%s\t%s\t%s\n'%(sample,fastafilepath,blastdbpath))
    #parse fasta
    recordlist=[] #list of id,seq tuples
    for recordid,recordseq in SeqIO.FastaIO.SimpleFastaParser(infile):
        #remove description from fasta id if present
        newfastaheader=re.sub(r'(\S+)(?: .*)?',r'\1',recordid)
        #replace | with _ to ensure unique contig names
        newfastaheader=re.sub(r'\|','_',newfastaheader.strip())
        #prefix with sample name to create recordid
        recordid='%s|%s'%(sample,newfastaheader)
        recordlist.append([recordid,recordseq])
    infile.close()
    #print fasta to stdout
    for recordid,recordseq in recordlist:
        f2.write('%s\t%s\n'%(str(recordid),str(len(recordseq))))

f2.close()
f3.close()
