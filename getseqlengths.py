import sys, gzip, os, re
from Bio import SeqIO

seqlenfile=sys.argv[1]
filepathinfo=sys.argv[2]

f2=open(seqlenfile,'w')
with open(filepathinfo) as f:
    for line in f:
        filepaths=line.strip().split('\t')
        fastafilepath=filepaths[1]
        #open fastafilepath as "infile" file handle
        if fastafilepath.endswith('.gz'):
            infile=gzip.open(fastafilepath,'rt')
        else:
            infile=open(fastafilepath,'r')
        #get sample name (=filename without extension)
        if fastafilepath.endswith('.gz'):
            gunzipfilepath=re.sub(r'\.gz$','',fastafilepath)
            samplename=os.path.basename(os.path.splitext(gunzipfilepath)[-2])
        else:
            samplename=os.path.splitext(os.path.basename(fastafilepath))[0]
        #parse fasta
        recordlist=[] #list of id,seq tuples
        for recordid,recordseq in SeqIO.FastaIO.SimpleFastaParser(infile):
            #remove description from fasta id if present
            newfastaheader=re.sub(r'(\S+)(?: .*)?',r'\1',recordid)
            #replace | with _ to ensure unique contig names
            newfastaheader=re.sub(r'\|','_',newfastaheader.strip())
            #prefix with samplename to create recordid
            recordid='%s|%s'%(samplename,newfastaheader)
            recordlist.append([recordid,recordseq])
        infile.close()
        #print fasta to stdout
        for recordid,recordseq in recordlist:
            f2.write('%s\t%s\n'%(str(recordid),str(len(recordseq))))
            
f2.close()


