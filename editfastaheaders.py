import sys, gzip, os, re
from Bio import SeqIO

inputtype=str(sys.argv[1])
if inputtype=='arg':
    filepath=str(sys.argv[2])
elif inputtype=='stdin':
    filepaths=sys.stdin
elif inputtype=='stdin_headersasis':
    filepaths=sys.stdin
else:
    sys.exit('unknown inputtype argument provided to editfastaheaders.py')

def editfasta(filepath):
    """parses fasta file from filepath stdin or filepath arg and return edited fasta with sample| prefixed headers"""
    #open filepath as "infile" file handle
    if filepath.endswith('.gz'):
        infile=gzip.open(filepath,'rt')
    else:
        infile=open(filepath,'r')
    #get sample name (=filename without extension)
    if filepath.endswith('.gz'):
        gunzipfilepath=re.sub(r'\.gz$','',filepath)
        sample=os.path.splitext(os.path.basename(gunzipfilepath))[0]
    else:
        sample=os.path.splitext(os.path.basename(filepath))[0]
    #parse fasta
    recordlist=[] #list of id,seq tuples
    for recordid,recordseq in SeqIO.FastaIO.SimpleFastaParser(infile):
        #remove description from fasta id if present
        newfastaheader=re.sub(r'(\S+)(?: .*)?',r'\1',recordid)
        #replace | with _ to ensure unique contig names
        newfastaheader=re.sub(r'\|','_',newfastaheader.strip())
        #prefix with sample to create recordid
        recordid='%s|%s'%(sample,newfastaheader)
        recordlist.append([recordid,recordseq])
    infile.close()
    #print fasta to stdout
    for recordid,recordseq in recordlist:
        print('>{}\n{}'.format(str(recordid),str(recordseq)))

if inputtype=='arg':
    editfasta(filepath)
if inputtype=='stdin':
    for line in filepaths:
        filepath=line.strip()
        editfasta(filepath)
if inputtype=='stdin_headersasis':
    for line in filepaths:
        filepath=line.strip()
         #open filepath as "infile" file handle
        if filepath.endswith('.gz'):
            infile=gzip.open(filepath,'rt')
        else:
            infile=open(filepath,'r')
        #print fasta to stdout
        for recordid,recordseq in SeqIO.FastaIO.SimpleFastaParser(infile):  #iterate through contigs
            print('>{}\n{}'.format(str(recordid),str(recordseq)))
        infile.close()



